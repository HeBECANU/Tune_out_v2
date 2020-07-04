% XSUM - Sum with error compensation [MEX]
% The accuracy of the sum of floating point numbers is limited by the
% truncation error. E.g. SUM([1e16, 1, -1e16]) replies 0 instead of 1. The
% error grows with the length of the input, e.g. the error of SUM(RANDN(N, 1))
% is about EPS*(N / 10).
% Kahan, Knuth, Dekker, Ogita and Rump (and others) have derived methods to
% reduce the influence of rounding errors. Some of them are implemented here as
% fast C-Mex.
%
% Y = XSum(X, N, Method)
% INPUT:
%   X: Double array of any size.
%   N: Dimension to operate on. Optional, default: [], which uses the first
%      non-singelton dimension.
%   Method: String: 'Double', 'Long', 'Kahan', 'Knuth', 'KnuthLong', 'Knuth2'.
%      Optional, default (and recommended): 'Knuth'.
%      Not case sensitive. A description of the methods follows.
%   Calling XSum without inputs displays the compilation date and availability
%   of long double methods.
%
% OUTPUT:
%   Y: Double array, equivalent to SUM, but with compensated error depending
%      on the Method. The high-precision result is rounded to double precision.
%      Y has the same size as X except for the N'th dimension, which is 1.
%      If X is empty, the empty matrix [] is replied.
%
% METHODS:
% The speed is compared to a single-threaded SUM for 1E6 elements. SUM is
% remarkably faster with multi-threading and less than 1E4 elements, so only
% the relations between the methods have an absolute meaning.
% The accuracy depends on the input, so the stated values are rules of thumb!
% Run uTest_XSum to test validity, compare accuracy and speed.
% Double: A thread-safe implementation of Matlab's SUM. At least in Matlab
%         2008a to 2009b the results of the multi-threaded SUM can differ
%         from call to call.
%         Speed:    120%
%         Accuracy: If the compilers supports accumulation in a 80 bit register
%                   (e.g. Open Watcom 1.8), 3 additional digits are gained.
%                   Otherwise the result is equivalent to SUM.
% Long:   The sum is accumulated in a 80 bit long double, if the compiler
%         supports this (e.g. LCC v3.8, Intel compilers).
%         Speed:    120% of SUM.
%         Accuracy: 3-4 more valid digits compared to SUM.
% Kahan:  The local error is subtracted from the next element. See [4].
%         Speed:    490%
%         Accuracy: 1 to 3 more valid digits than SUM.
% Knuth:  Calculate the sum as if it is accumulated in a 128 bit float.
%         See [1], [2], [3].
%         Speed:    250%
%         Accuracy: About 15 more valid digits than SUM.
% Knuth2: Calculate the sum as if it is accumulated in a 192 bit float. This is
%         the Knuth method with a further level of error compensation. See [2]
%         and INTLAB.
%         Speed:    450%
%         Accuracy: 30 more valid digits than SUM.
% KnuthLong: As Knuth, but the temporary variables are long doubles, if this
%         is supported by the compiler. If not available, Knuth2 is called.
%         Speed:    350%
%         Accuracy: 21 more valid digits than SUM.
%
% COMMENTS:
%   The accuracy of 'Double', 'Long' and 'KnuthLong' depends on the compiler and
%   the SSE2 and floating point environment flags! In opposite to that, 'Kahan',
%   'Knuth', 'Knuth2' are using standard precision IEEE arithmetics only and
%   reply equivalent values on different platforms.
%
%   'Double', 'Long', 'Kahan' and 'KnuthLong' are included for accademic
%   reasons. For real world problems 'Knuth' is recommended, because it doubles
%   the accuracy without wasting too much time. The higher accuracy of 'Knuth2'
%   is rarely needed (I've used it to check the results of the Knuth method).
%
%   Sorting the input in ascending order increases the accuracy *sometimes*, but
%   it takes a lot of time: e.g. 30 times longer for 1E7 elements. Consider that
%   sorting can decrease the accuracy, when the values have different signs!
%
% COMPILATION:
% This function must be compiled before using. This happens automatically, when
% the M-function is called the first time and the function InstallMex is found.
% See XSum.c for detailed instructions for compiling manually.
%
% REFERENCES, MORE COMMENTS: See XSum.c
%
% Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2010-2014 matlab.THISYEAR(a)nMINUSsimon.de
%
% See also SUM.
% FEX: Alain Barraud: SUMBENCH (#23198), FORWARD_STABLE_SOLVER (#10668),
%      SUMDOTPACK (#8765).
% INTLAB: Prof. Dr. Siegfried M. Rump, http://www.ti3.tu-harburg.de/rump/intlab
% DASUM of XBLAS (see www.netlib.org)

% $JRev: R-u V:020 Sum:iyOaTRus2c+Z Date:01-Mar-2014 12:07:33 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_XSum $
% $File: Tools\GLMath\XSum.m $

% This M-file is a dummy only to carry teh help section and to start the
% automatic compilation.

% This is a dummy code only, which compiles the C-function automatically:
function varargout = XSum(varargin)

Ok = InstallMex('XSum.c', 'uTest_XSum');
if ~Ok
   error('JSimon:XSum:BadCompilation', ...
      'Compilation failed. See help text of XSum.c');
end

[varargout{1:nargout}] = XSum(varargin{:});  % Hopefully not a recursion...
