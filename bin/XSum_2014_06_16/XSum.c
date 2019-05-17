// XSum.c
// XSUM - Sum with error compensation [MEX]
// The accuracy of the sum of floating point numbers is limited by the
// truncation error. E.g. SUM([1e16, 1, -1e16]) replies 0 instead of 1. The
// error grows with the length of the input, e.g. the error of SUM(RANDN(N, 1))
// is about EPS*(N / 10).
// Kahan, Knuth, Dekker, Ogita and Rump (and others) have derived methods to
// reduce the influence of rounding errors. Some of them are implemented here as
// fast C-Mex.
//
// Y = XSum(X, N, Method)
// INPUT:
//   X: Double array of any size.
//   N: Dimension to operate on. Optional, default: [], which uses the first
//      non-singelton dimension.
//   Method: String: 'Double', 'Long', 'Kahan', 'Knuth', 'KnuthLong', 'Knuth2'.
//      Optional, default (and recommended): 'Knuth'.
//      Not case sensitive. A description of the methods follows.
//   Calling XSum without inputs displays the compilation date and availability
//   of long double methods.
//
// OUTPUT:
//   Y: Double array, equivalent to SUM, but with compensated error depending
//      on the Method. The high-precision result is rounded to double precision.
//      Y has the same size as X except for the N'th dimension, which is 1.
//      If X is empty, the empty matrix [] is replied.
//
// METHODS:
// The speed is compared to a single-threaded SUM for 1E6 elements. SUM is
// remarkably faster with multi-threading and less than 1E4 elements, so only
// the relations between the methods have an absolute meaning.
// The accuracy depends on the input, so the stated values are rules of thumb!
// Run uTest_XSum to test validity, compare accuracy and speed.
// Double: A thread-safe implementation of Matlab's SUM. At least in Matlab
//         2008a to 2009b the results of the multi-threaded SUM can differ
//         from call to call.
//         Speed:    120%
//         Accuracy: If the compilers supports accumulation in a 80 bit register
//                   (e.g. Open Watcom 1.8), 3 additional digits are gained.
//                   Otherwise the result is equivalent to SUM.
// Long:   The sum is accumulated in a 80 bit long double, if the compiler
//         supports this (e.g. LCC v3.8, Intel compilers).
//         Speed:    120% of SUM.
//         Accuracy: 3-4 more valid digits compared to SUM.
// Kahan:  The local error is subtracted from the next element. See [4].
//         Speed:    490%
//         Accuracy: 1 to 3 more valid digits than SUM.
// Knuth:  Calculate the sum as if it is accumulated in a 128 bit float.
//         See [1], [2], [3].
//         Speed:    250%
//         Accuracy: About 15 more valid digits than SUM.
// Knuth2: Calculate the sum as if it is accumulated in a 192 bit float. This is
//         the Knuth method with a further level of error compensation. See [2]
//         and INTLAB.
//         Speed:    450%
//         Accuracy: 30 more valid digits than SUM.
// KnuthLong: As Knuth, but the temporary variables are long doubles, if this
//         is supported by the compiler. If not available, Knuth2 is called.
//         Speed:    350%
//         Accuracy: 21 more valid digits than SUM.
//
// COMMENTS:
//   The accuracy of 'Double', 'Long' and 'KnuthLong' depends on the compiler
//   and the SSE2 and floating point environment flags! In opposite to that,
//   'Kahan', 'Knuth', 'Knuth2' are using standard precision IEEE arithmetics
//   only and reply equivalent values on different platforms.
//
//   'Double', 'Long', 'Kahan' and 'KnuthLong' are included for accademic
//   reasons. For real world problems 'Knuth' is recommended, because it doubles
//   the accuracy without wasting too much time. The higher accuracy of 'Knuth2'
//   is rarely needed (I've used it to check the results of the Knuth method).
//
//   Sorting the input in ascending order increases the accuracy *sometimes*,
//   but it takes a lot of time: e.g. 30 times longer for 1E7 elements. Consider
//   that sorting can decrease the accuracy, when the values have different
//   signs!
//
//   I've built a version using 384 bit QFLOATs with the free LCC v3.8
//   compiler (not the v2.4 shipped with Matlab): 104 digits precision and 50
//   times slower than SUM. Please mail me, if you want this version.
//
//   The subfunction DoubleDim1/N/V are a nice example for the operating on
//   subvectors in arbitrary dimensions. This is more memory-efficient and in
//   consequence much faster than using PERMUTE (e.g. in GRADIENT).
//
//   On demand I publish an XSum, which calculates the sum of the subvectors in
//   different threads. But I will not use threads for the calculation of
//   partial sums for large vectors (as in Matlab's SUM): this would decrease
//   the numerical stability!
//
// COMPILE:
// This C-file is compiled automatically when the M-function is called the
// first time the function InstallMex is found.
// You can compile XSum.c manually also:
//   mex -O XSum.c
// Linux: Consider C99 comments:
//   mex -O CFLAGS="\$CFLAGS -std=c99" XSum.c
// Pre-compiled MEX files:
//   http://www.n-simon.de/mex
// Run the unit-test uTest_XSum after the compilation.
//
// REFERENCES:
// [1] D.E. Knuth: "The Art of Computer Programming: Seminumerical
//   Algorithms", volume 2. Addison Wesley, Reading, Massachusetts,
//   second edition, 1981.
//
// [2] Takeshi Ogita and Siegfried M. Rump and Shin'ichi Oishi:
//   "Accurate Sum and Dot Product with Applications"
//   Proceedings of 2004 IEEE International Symposium on
//   Computer Aided Control Systems Design, Taipei, pages 152155, 2004
//   NOTE: Error in PDF "OgRuOi04a.pdf" (ask Google), TwoSum algorithm:
//     y = fl((a-(x-z))+(b + z))  ==>  y = fl((a-(x-z))+(b - z))
//
// [3] T.J. Dekker: "A Floating-Point Technique for Extending the
//   Available Precision". Numerische Mathematik, 18:224-242, 1971.
//
// [4] Linnainmaa, S.: Analysis of Some Known Methods of Improving the Accuracy
//   of Floating-point Sums, BIT 14, (1974), 167-202.
//
// Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
//         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008/2010
// Assumed Compatibility: higher Matlab versions, Mac, Linux
// Author: Jan Simon, Heidelberg, (C) 2010-2014 matlab.THISYEAR(a)nMINUSsimon.de
//
// See also: SUM.
// FEX: Alain Barraud: SUMBENCH (#23198), FORWARD_STABLE_SOLVER (#10668),
//      SUMDOTPACK (#8765)
// INTLAB: Prof. Dr. Siegfried M. Rump, http://www.ti3.tu-harburg.de/rump/intlab
// DASUM of XBLAS (see www.netlib.org)

/*
% $JRev: R-E V:030 Sum:fjgiT/MGaDqp Date:16-Jun-2014 00:12:37 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_XSum $
% $File: Tools\Mex\Source\XSum.c $
% History:
% 001: 12-Feb-2010 21:51, First version: Sum with error compensations.
%      I built a version with QFloats of LCC v3.8.
%      XDot is under construction...
% 015: 28-Mar-2010 00:54, Help section and comments improved.
% 023: 30-Nov-2010 11:43, No _control87 under Linux and MacOS.
% 029: 19-Jan-2014 00:01, Brushed up the help section and error handling.
*/

#include "mex.h"
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

// Assume 32 bit addressing for Matlab 6.5:
// See MEX option "compatibleArrayDims" for MEX in Matlab >= 7.7.
#ifndef MWSIZE_MAX
#define mwSize  int32_T           // Defined in tmwtypes.h
#define mwIndex int32_T
#define MWSIZE_MAX MAX_int32_T
#endif

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_ID   "JSimon:XSum:"
#define ERR_HEAD "*** XSum[mex]: "
#define ERROR(id,msg) mexErrMsgIdAndTxt(ERR_ID id, ERR_HEAD msg);

// Under Windows some compilers use 80 bit for LONG DOUBLE. Then setting the
// floating point environment on Intel processors influences the results:
#if defined(__WINDOWS__) || defined(WIN32) || defined(_WIN32) || defined(_WIN64)
#define SET_PRECISION(C) SetPrecision(C);
#define DO_SET_PRECISION 1
#else
#define SET_PRECISION(C)  // No operation under Linux and Mac
#define DO_SET_PRECISION 0
#endif


// Case-insensitive string comparison depends on the compiler: -----------------
// strncmpi, strnicmp, _strnicmp, ...
#define STRCMP_NI _strnicmp


// Prototypes: -----------------------------------------------------------------
void hello(void);
mwSize GetStep(const mwSize *XDim, const mwSize N);
mwSize FirstNonSingeltonDim(const mwSize Xndim, const mwSize *Xdim);

void KnuthDim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void KnuthDimN(double *X, const mwSize *Xdim, const mwSize nX,
               const mwSize N, double *Y);
double KnuthV(double *X, double *XEnd, mwSize Step);

void Knuth2Dim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void Knuth2DimN(double *X, const mwSize *Xdim, const mwSize nX,
                const mwSize N, double *Y);
double Knuth2V(double *X, double *XEnd, mwSize Step);

void KnuthLDim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void KnuthLDimN(double *X, const mwSize *Xdim, const mwSize nX,
                const mwSize N, double *Y);
double KnuthLV(double *X, double *XEnd, mwSize Step);

void KahanDim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void KahanDimN(double *X, const mwSize *Xdim, const mwSize nX,
               const mwSize N, double *Y);
double KahanV(double *X, double *XEnd, mwSize Step);

void DoubleDim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void DoubleDimN(double *X, const mwSize *Xdim, const mwSize nX,
                const mwSize N, double *Y);

void LongDim1(double *X, const mwSize nX, const mwSize rX, double *Y);
void LongDimN(double *X, const mwSize *Xdim, const mwSize nX,
              const mwSize N, double *Y);

#if DO_SET_PRECISION
void SetPrecision(int Command);
#endif

// List of methods:
typedef enum {Unknown_SUM, Double_SUM, Long_SUM, Kahan_SUM,
              Knuth_SUM, KnuthLong_SUM, Knuth2_SUM} Method_t;

// Length of the method names (longest + 1):
#define Method_LEN 10

// Flag to show the warning for missing LONG support once only:
static int NoLong_ShowWarning = 1;

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mwSize *Xdim;
  mwSize       nX, *Ydim, Xndim, N;
  int          dim1Mode;
  Method_t     Method = Knuth_SUM;
  char         Method_In[Method_LEN];
  double       Nd;
  
  // Show general information:
  if (nrhs == 0) {
     hello();
     return;
  }
  
  // Check number of inputs and outputs:
  if (nrhs > 3) {
     ERROR("BadNInput", "1 to 3 inputs required.");
  }
  if (nlhs > 1) {
     ERROR("BadNOutput", "1 output allowed.");
  }
  
  // Input must be a double array - is it worth to implement a sum with higher
  // accuracy for SINGLEs ?!
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
     ERROR("BadDataType", "X must be a real double array.");
  }
  
  // Get first argument:
  nX    = mxGetNumberOfElements(prhs[0]);
  Xdim  = mxGetDimensions(prhs[0]);
  Xndim = mxGetNumberOfDimensions(prhs[0]);
  if (nX == 0) {   // Return the empty matrix if the input is empty:
     plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
     return;
  }
  
  // Get 2nd input [Dim] - default is the first non-singelton dimension:
  if (nrhs >= 2) {
     switch (mxGetNumberOfElements(prhs[1])) {
        case 0:  // Input [Dim] is empty - first non-singelton dimension:
           N        = FirstNonSingeltonDim(Xndim, Xdim);
           dim1Mode = 1;
           break;
           
        case 1:
           // mxGetScalar needs a numeric input:
           if (!mxIsNumeric(prhs[1])) {
              ERROR("BadDimType",
                    "2nd input [Dim] must be a numeric scalar or empty.");
           }
           
           Nd = mxGetScalar(prhs[1]);
           N  = (mwSize) Nd - 1;  // Zero based index
           if (N >= Xndim || Nd <= 0.0 || Nd != floor(Nd)) {
              ERROR("BadDim", "Dimension N out of range.");
           }
           
           // The fast algorithm for the 1st dimension works also for [1 x M]
           // row vectors. [1 x 1 x M] arrays are not caught.
           dim1Mode = (N == 0 || (Xndim == 2 && N == 1 && *Xdim == 1));
           break;
           
        default:
           ERROR("BadDimType",
                 "2nd input [Dim] must be a numeric scalar or empty.");
     }
     
  } else {  // nrhs == 1: Operate on 1st non-singelton dimension:
     N        = FirstNonSingeltonDim(Xndim, Xdim);
     dim1Mode = 1;
  }
  
  // The sum over a scalar is a scalar:
  if (Xdim[N] == 1) {
     plhs[0] = mxDuplicateArray(prhs[0]);
     return;
  }
  
  // Create output:
  if ((Ydim = (mwSize *) mxMalloc(Xndim * sizeof(mwSize))) == NULL) {
     ERROR("OutOfMemory", "Cannot get memory for ndim.");
  }
  memcpy(Ydim, Xdim, Xndim * sizeof(mwSize));
  Ydim[N] = 1;
  plhs[0] = mxCreateNumericArray(Xndim, Ydim, mxDOUBLE_CLASS, mxREAL);
  mxFree(Ydim);
  
  // Parse 3rd input, default is KNUTH:
  if (nrhs == 3) {
     if (!mxIsChar(prhs[2])) {
        ERROR("BadTypeInput3", "3rd input must be a string.");
     }
     
     mxGetString(prhs[2], Method_In, Method_LEN);
     if (STRCMP_NI(Method_In, "Knuth", Method_LEN) == 0) {
        Method = Knuth_SUM;   // The default
     } else if (STRCMP_NI(Method_In, "Knuth2", Method_LEN) == 0) {
        Method = Knuth2_SUM;
     } else if (STRCMP_NI(Method_In, "KnuthLong", Method_LEN) == 0) {
        if (sizeof(long double) > sizeof(double)) {
           Method = KnuthLong_SUM;
        } else {              // Warn once if long doubles are doubles:
           if (NoLong_ShowWarning == 1) {
              mexWarnMsgIdAndTxt("JSim:XSum:NoLongDouble",
                "Using Knuth2 instead of Knuth_Long: "
                "XSum compiled without long double support!");
              NoLong_ShowWarning = 0;
           }
           Method = Knuth2_SUM;
        }
     } else if (STRCMP_NI(Method_In, "Double", Method_LEN) == 0) {
        Method = Double_SUM;
     } else if (STRCMP_NI(Method_In, "Long", Method_LEN) == 0) {
        Method = Long_SUM;
     } else if (STRCMP_NI(Method_In, "Kahan", Method_LEN) == 0) {
        Method = Kahan_SUM;
     } else {             // Bad name, caught later
        Method = Unknown_SUM;
     }
  }
  
  // Call the core calculator for the specific method:
  switch (Method) {
     case Knuth_SUM:
        if (dim1Mode) {
           KnuthDim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           KnuthDimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;
        
     case Knuth2_SUM:
        if (dim1Mode) {
           Knuth2Dim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           Knuth2DimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;
        
     case KnuthLong_SUM:
        if (dim1Mode) {
           KnuthLDim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           KnuthLDimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;

     case Kahan_SUM:
        if (dim1Mode) {
           KahanDim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           KahanDimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;
        
     case Double_SUM:
        if (dim1Mode) {
           DoubleDim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           DoubleDimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;

     case Long_SUM:
        if (dim1Mode) {
           LongDim1(mxGetPr(prhs[0]), nX, Xdim[N], mxGetPr(plhs[0]));
        } else {
           LongDimN(mxGetPr(prhs[0]), Xdim, nX, N, mxGetPr(plhs[0]));
        }
        break;
        
     default:
        ERROR("BadValueInput3", "3rd input [Method] not recognized.");
  }
  
  return;
}

// =============================================================================
mwSize FirstNonSingeltonDim(const mwSize Xndim, const mwSize *Xdim)
{
  // Get first non-singelton dimension - zero based.
  mwSize N;
  
  for (N = 0; N < Xndim; N++) {
     if (Xdim[N] != 1) {
        return (N);
     }
  }
  
  return (0);  // Use the first dimension if all dims are 1
}

// =============================================================================
mwSize GetStep(const mwSize *Xdim, const mwSize N)
{
  // Get step size between elements of a subvector in the N'th dimension.
  const mwSize *XdimEnd, *XdimP;
  mwSize       Step;
  
  Step    = 1;
  XdimEnd = Xdim + N;
  for (XdimP = Xdim; XdimP < XdimEnd; Step *= *XdimP++) ; // empty loop
  
  return (Step);
}

// =============================================================================
void KnuthDim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate Knuth's sum over vectors in the 1st dimension.
  // Inputs/outputs: See DoubleDim1 and DoubleDimN
  
  // METHOD:
  //   % Sum(X): See reference [2]
  //     n = length(X);
  //     for i = 2:n
  //       [X(i), X(i-1)] = TwoSum(X(i), X(i-1))
  //     end
  //     Sum = sum(X(1:n-1)) + X(n);
  //
  //   % 6 ops, better optimization without branch, [1]:  (This is used!)
  //   function1 [x, y] = TwoSum(a, b)
  //     x = fl(a + b)
  //     z = fl(x - a)
  //     y = fl((a - (x - z)) + (b - z))
  //
  //   % 4 ops, but branching impedes pipelining in the processor, [3]:
  //   function2 [x, y] = TwoSum(a, b)
  //     x = fl(a + b)
  //     if abs(a) >= abs(b)
  //       y = fl(b - (x - a))
  //     else
  //       y = fl(a - (x - b))
  
  double *XR, *XEnd, h, q, x, z;
  
  XEnd = X + nX;
  while (X < XEnd) {
     XR = X + rX;
     h  = *X++;
     q  = 0.0;
     while (X < XR) {
        x  = h + *X;                        // Predictor
        z  = x - h;                         // Corrector
        q += ((h - (x - z)) + (*X++ - z));  // Accumulate the local error
        h  = x;
     }
     *Y++ = h + q;                          // Add accumulated local errors
  }
  
  return;
}

// -----------------------------------------------------------------------------
void KnuthDimN(double *X, const mwSize *Xdim, const mwSize nX,
               const mwSize N, double *Y)
{
  // Calculate Knuth's sum over vectors in N'th dimension.
  // Inputs/outputs: See DoubleDim1 and DoubleDimN

  double *XEnd, *SliceEnd;
  mwSize Step, SumStep;
  
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];
  
  // Loop over all dimension:
  XEnd = X + nX;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // X goes through the elements of the subvector in the N'th dimension:
        *Y++ = KnuthV(X, X + SumStep, Step);
     }
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  return;
}

// -----------------------------------------------------------------------------
double KnuthV(double *X, double *XEnd, mwSize Step)
{
  // Calculate Knuth's sum for a subvector.
  // [Step] is the distance between elements of the subvector.
  double q, h, x, z;

  h  = *X;
  X += Step;
  q  = 0.0;
  while (X < XEnd) {
     x  = h + *X;
     z  = x - h;
     q += ((h - (x - z)) + (*X - z));
     h  = x;
     X += Step;
  }

  return (h + q);
}

// =============================================================================
void Knuth2Dim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate Knuth's sum over vectors in the 1st dimension, accumulate the
  // local errors with one extra stage of correction. Finally the sum is
  // obtained as if it was accumulated in a 196 bit double.
  // Inputs/outputs: See DoubleDim1 and DoubleDimN

  double *XR, *XEnd, h, q1, q2, x, z, x2, z2, e1;
  
  XEnd = X + nX;
  while (X < XEnd) {
     XR = X + rX;
     h  = *X++;
     q1  = 0.0;
     q2  = 0.0;
     while (X < XR) {
        x = h + *X;
        z = x - h;
        
        // Accumulate the local error with a further level of correction:
        e1  = ((h - (x - z)) + (*X++ - z));
        x2  = q1 + e1;
        z2  = x2 - q1;
        q2 += ((q1 - (x2 - z2)) + (e1 - z2));
        q1  = x2;
        
        h = x;
     }
     *Y++ = h + q1 + q2;
  }
  
  return;
}

// -----------------------------------------------------------------------------
void Knuth2DimN(double *X, const mwSize *Xdim, const mwSize nX,
                const mwSize N, double *Y)
{
  // Calculate Knuth's sum over vectors in N'th dimension.
  // Inputs/outputs: See DoubleDim1 and DoubleDimN
  
  double *XEnd, *SliceEnd;
  mwSize Step, SumStep;
  
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];
  
  // Loop over all dimension:
  XEnd = X + nX;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // X goes through the elements of the subvector in the N'th dimension:
        *Y++ = Knuth2V(X, X + SumStep, Step);
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  return;
}

// -----------------------------------------------------------------------------
double Knuth2V(double *X, double *XEnd, mwSize Step)
{
  // Calculate Knuth's sum over vectors in the 1st dimension, accumulate the
  // local errors with one extra stage of correction. Finally the sum is
  // obtained as if it was accumulated in a 196 bit double.
  // [Step] is the distance between elements of the subvector.
  double h, q1, q2, x, z, x2, z2, e1;
  
  h   = *X;
  X  += Step;
  q1  = 0.0;
  q2  = 0.0;
  while (X < XEnd) {
     x = h + *X;
     z = x - h;
     
     // Accumulate the local error with a further level of correction:
     e1  = ((h - (x - z)) + (*X - z));
     x2  = q1 + e1;
     z2  = x2 - q1;
     q2 += ((q1 - (x2 - z2)) + (e1 - z2));
     q1  = x2;
     
     h  = x;
     X += Step;
  }
  
  return (h + (q1 + q2));
}

// =============================================================================
void KnuthLDim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate Knuth's sum over vectors in the 1st dimension with intermediate
  // values stored in long doubles.
  // If the compiler uses doubles for long doubles, the error correction fails
  // completely due to setting the precision to 64 bit mantissa. Then the
  // result are comparable with the uncompensated SUM!
  // Inputs/outputs: See DoubleDim1 and DoubleDimN

  double *XR, *XEnd;
  long double h, q, x, z;
  
  SET_PRECISION(64);
  
  XEnd = X + nX;
  while (X < XEnd) {
     XR = X + rX;
     h  = (long double) *X++;
     q  = 0.0L;
     while (X < XR) {
        x  = h + *X;
        z  = x - h;
        q += ((h - (x - z)) + (*X++ - z));
        h  = x;
     }
     *Y++ = (double) (h + q);
  }
  
  SET_PRECISION(0);
  return;
}

// -----------------------------------------------------------------------------
void KnuthLDimN(double *X, const mwSize *Xdim, const mwSize nX,
                const mwSize N, double *Y)
{
  // Calculate Knuth's sum over vectors in N'th dimension with long doubles.
  // This is identical to KnuthDimN, just the subfunction for processing the
  // subvectors differs.
  // Inputs/outputs: See DoubleDim1 and DoubleDimN

  double *XEnd, *SliceEnd;
  mwSize Step, SumStep;
  
  SET_PRECISION(64);
  
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];
  
  // Loop over all dimension:
  XEnd = X + nX;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // X goes through the elements of the subvector in the N'th dimension:
        *Y++ = KnuthLV(X, X + SumStep, Step);
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  SET_PRECISION(0);
  return;
}

// -----------------------------------------------------------------------------
double KnuthLV(double *X, double *XEnd, mwSize Step)
{
  // Calculate Knuth's sum for a subvector with long doubles.
  // [Step] is the distance between elements of the subvector.
  // Caller set precision to 64 bits mantissa.
  
  long double q, h, x, z;
  
  h  = (long double) *X;
  X += Step;
  q  = 0.0L;
  while (X < XEnd) {
     x  = h + *X;
     z  = x - h;
     q += (h - (x - z)) + (*X - z);
     h  = x;
     X += Step;
  }
  
  return ((double) (h + q));
}

// =============================================================================
void KahanDim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate Kahan's sum over vectors of the 1st dimension.
  // Inputs/outputs: See DoubleDim1 and DoubleDimN

  double *XR, *XEnd, c, t, x, s;
  
  // Kahan sum does not work with 64 bit precision, because this conceals the
  // cancellation error!!! 53 bit precision is the default in Matlab, but the
  // user could have changed this by "system_dependent('setprecision', 64)".
  SET_PRECISION(53);
  
  XEnd = X + nX;
  while (X < XEnd) {
     XR = X + rX;
     s  = *X++;
     c  = 0.0;
     while (X < XR) {
       x = *X++ - c;
       t = s + x;
       c = (t - s) - x;
       s = t;
     }
     
     *Y++ = s - c;
  }
  
  SET_PRECISION(0);
  return;
}

// -----------------------------------------------------------------------------
void KahanDimN(double *X, const mwSize *Xdim, const mwSize nX,
               const mwSize N, double *Y)
{
  // Calculate Kahan's sum over vectors in N'th dimension.
  // Inputs/outputs: See DoubleDim1 and DoubleDimN

  double *XEnd, *SliceEnd;
  mwSize Step, SumStep;
  
  SET_PRECISION(53);  // Not working with 64!
    
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];

  // Loop over all dimension:
  XEnd = X + nX;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // X goes through the elements of the subvector in the N'th dimension:
        *Y++ = KahanV(X, X + SumStep, Step);
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  SET_PRECISION(0);
  return;
}

// -----------------------------------------------------------------------------
double KahanV(double *X, double *XEnd, mwSize Step)
{
  // Calculate Kahan's sum for a subvector.
  // [Step] is the distance between elements of the subvector.
  double c = 0.0, s, t, x;
  
  // Cheap first step:
  s  = *X;
  X += Step;
  
  // Loop over elements 2 to end:
  while(X < XEnd) {
    x  = *X - c;
    t  = s + x;
    c  = (t - s) - x;
    s  = t;
    X += Step;
  }
  
  return (s);
}

// =============================================================================
void DoubleDim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate sum over vectors of the 1st dimension in double precision.
  // X:  Pointer to input array
  // nX: Number of elements of X
  // rX: Length or first (non-singleton) dimension
  // Y:  Pointer to output array
  
  double *XR, *XEnd;
  double Sum;
  
  // With enabled usage of 80 bit registers (just some compilers) the result
  // can be more accurate than SUM:
  SET_PRECISION(64);

  XEnd = X + nX;
  while (X < XEnd) {
     Sum = 0.0;
     for (XR = X + rX; X < XR; Sum += *X++) ;  // empty loop
     *Y++ = Sum;
  }
  
  SET_PRECISION(0);
  return;
}

// -----------------------------------------------------------------------------
void DoubleDimN(double *X, const mwSize *Xdim, const mwSize nX,
                const mwSize N, double *Y)
{
  // Calculate sum over vectors in N'th dimension with long double precision.
  // X:  Pointer to input array
  // Xdim: Vector of dimensions of X
  // nX: Number of elements of X
  // N:  Dimension to operate on, zero base.
  // Y:  Pointer to output array
  
  double *XEnd, *SliceEnd, *V, *VEnd;
  mwSize Step, SumStep;
  double Sum;
  
  SET_PRECISION(64);
  
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];
  
  // Loop over all dimension:
  XEnd = X + nX;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // V goes through the elements of the subvector in the N'th dimension:
        Sum  = 0.0L;
        VEnd = X + SumStep;
        for (V = X; V < VEnd; V += Step) {
           Sum += *V;
        }
        *Y++ = Sum;
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  SET_PRECISION(0);
  return;
}

// =============================================================================
void LongDim1(double *X, const mwSize nX, const mwSize rX, double *Y)
{
  // Calculate sum over vectors of the 1st dimension in long double precision.
  //
  // This algorithm can benefit from the 64 bit precision even if the compiler
  // uses doubles as long doubles: [Sum] is accumulated in a 80 bit register!
  // This happens e.g. with Open Watcom 1.8, but neither with LCC v2.4, which is
  // shipped with Matlab, nor with MSVC++ 2008.
  // Inputs/outputs: See DoubleDim1 and DoubleDimN
    
  double *XR, *XEnd;
  long double Sum;
  
  SET_PRECISION(64);

  XEnd = X + nX;
  while (X < XEnd) {
     Sum = 0.0L;
     for (XR = X + rX; X < XR; Sum += *X++) ;  // empty loop
     *Y++ = (double) Sum;
  }
  
  SET_PRECISION(0);
  return;
}

// -----------------------------------------------------------------------------
void LongDimN(double *X, const mwSize *Xdim, const mwSize nX,
              const mwSize N, double *Y)
{
  // Calculate sum over vectors in N'th dimension with long double precision.
  // For the effects of long doubles see LongDim1.
  // Inputs/outputs: See DoubleDim1 and DoubleDimN
   
  double *XEnd, *SliceEnd, *V, *VEnd;
  mwSize Step, SumStep;
  long double Sum;
  
  SET_PRECISION(64);
  
  // Get step size between elements of the subvector and its extent:
  Step    = GetStep(Xdim, N);
  SumStep = Step * Xdim[N];
  
  // Loop over all dimension:
  XEnd = X + nX;
  while (X < XEnd) {
     for (SliceEnd = X + Step; X < SliceEnd; X++) {
        // V goes through the elements of the subvector in the N'th dimension:
        Sum  = 0.0L;
        VEnd = X + SumStep;
        for (V = X; V < VEnd; V += Step) {
           Sum += *V;
        }
        *Y++ = (double) Sum;
     }  // end for
     
     // Move pointer to the next slice:
     X += SumStep - Step;
  }
 
  SET_PRECISION(0);
  return;
}

// =============================================================================
// I've heared that Linux machines always drive in 80 bit mode?!
// At least _control87 is not available under Linux and MacOS. But even under
// MSVC++ 2008 the 80 bit precision does not change the results. Therefore I
// think it is ok to disable the precision control under Linux.
#if DO_SET_PRECISION

void SetPrecision(int Command)
{
  // Set floating point environment for Intel compatible processors.
  // The long double calculations work with 80 bit numbers only, if the
  // precision is set to 64 bits. With 53 bits, the calculations have double
  // precision only.
  // The Kahan sum fails when processed with 64 bit precision, because this
  // steals the necessary bits from the correction of the local error.
  // The Knuth and Knuth2 methods do not profit from setting the precision.
  //
  // It seems that Matlab resets the floating point processor flags at the
  // standard exit of a Mex function. But at least Matlab 5.3 and 6.5 leave the
  // the flags untouched on exits through mexErrMsgTxt. So care for resetting
  // the precision in case of errors!
  //
  // ONLY TESTED FOR: BCC5.5, LCC2.4/3.8, Opwn Watcom 1.8, MSVC 2008
  // I have no idea how this works on processors, which are not supported by
  // the function _control87!
  // If SetPrecision causes troubles on a non-Intel processor, simply delete it!
  // The recommended methods Knuth and Knuth2 do not even use it.
  //
  // INPUT:
  //   Command: int, precision of floating point calculations.
  //            Allowed values: 24 (single), 53 (double), 64 (long double),
  //                            0 (reset to previous value)
  
#if defined(__LCC__)           // Leading underscores needed for macros:
#define MCW_PC _MCW_PC
#define PC_24  _PC_24
#define PC_53  _PC_53
#define PC_64  _PC_64
#endif

// I've experimented with this - without success here:
// MSVC++ 2008: #pragma fenv_access (on)
  
  static unsigned int PC_Back = 0;
  
  if (Command == 0) {                   // Restore original settings
     _control87(PC_Back, MCW_PC);
  } else {
     PC_Back = _control87(MCW_PC, 0);   // Backup original settings
     if (Command == 24) {
        _control87(PC_24, MCW_PC);      // Single precision
     } else if (Command == 53) {
        _control87(PC_53, MCW_PC);      // Double precision
     } else if (Command == 64) {
        _control87(PC_64, MCW_PC);      // Full precision
     } else {                           // Programming error:
        ERROR("BadPrecCommand", "Bad command in SetPrecision.");
     }
  }

  return;
}
#endif  // DO_SET_PRECISION

// =============================================================================
void hello(void)
{
  // Show some information.
  // Most important: availability of LONG DOUBLE
  
  char *hasLong, *yes = "supported", *no = "not supported";
  
  if (sizeof(long double) > sizeof(double)) {
     hasLong = yes;
  } else {
     hasLong = no;
  }
  
  // Actually the compiler would be important also, but it is not available
  // as macro.
  mexPrintf("XSum: Sum with compensated errors\n"
            "Author: Jan Simon\n"
#if defined(__DATE__) && defined(__TIME__)
            "Compiled: " __DATE__ " " __TIME__ "\n"
#endif
            "Long double: %s\n", hasLong);
  
  return;
}
