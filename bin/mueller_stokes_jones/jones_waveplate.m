function JM = jones_waveplate(varargin)
    %function JM = jones_waveplate(varargin)
    %Return the Jones matrix for a linear wave plate with a phase delay given in wavelength units and long axis rotation of 0 degrees.
    %
    %  Syntax:
    %     JM = jones_waveplate()
    %     JM = jones_waveplate(p)
    %     A = jones_waveplate(P)
    %
    %  Description:
    %     JM = jones_waveplate() returns the Jones matrix
    %     for a wave plate with zero phase delay.
    %
    %     JM = jones_waveplate(p) returns the Jones matrix
    %     for a wave plate with phase delay p in wavelength
    %     units, i.e. p is ranging between 0 and 1.
    %
    %     A = jones_waveplate(P), where P is a m-by-n-by-... array of
    %     phase delay values, returns a cell array of
    %     wave plate Jones matrices with size(A)==size(P).
    %
    %  Example:
    %     JM = jones_waveplate(1/6) % Jones matrix for sixth wave plate
    %     JM =
    %
    %        1.00000 + 0.00000i   0.00000 + 0.00000i
    %        0.00000 + 0.00000i   0.50000 - 0.86603i
    %
    %  References:
    %     [1] E. Collett, Field Guide to Polarization,
    %         SPIE Field Guides vol. FG05, SPIE (2005). ISBN 0-8194-5868-6.
    %     [2] R. A. Chipman, "Polarimetry," chapter 22 in Handbook of Optics II,
    %         2nd Ed, M. Bass, editor in chief (McGraw-Hill, New York, 1995)
    %     [3] "Jones calculus", http://en.wikipedia.org/wiki/Jones_calculus,
    %         last retrieved on Jan 13, 2014.
    %
    %  See also:
    %     jones_linretarder  jones_circretarder
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    pilu_defv = 0;
    
    if nargin<1
        phase_in_lambda_units = pilu_defv;
    else
        phase_in_lambda_units = varargin{1};
    end
    
    [phase_in_lambda_units, was_cell] = c2n(phase_in_lambda_units, pilu_defv);
    
    if (numel(phase_in_lambda_units) > 1) || was_cell
        
        JM = cell(size(phase_in_lambda_units));
        JM_subs = cell(1,ndims(JM));
        for jmi=1:numel(JM)
            [JM_subs{:}] = ind2sub(size(JM),jmi);
            JM{JM_subs{:}} = s_waveplate(phase_in_lambda_units(JM_subs{:}));
        end
        
    else
        
        JM = s_waveplate(phase_in_lambda_units);
        
    end
    
end

% helper function
function JM = s_waveplate(phase_in_lambda_units)
    
    JM = zeros(2,2);
    JM(1,1) = 1;
    JM(2,2) = exp(-1i*phase_in_lambda_units*2*pi());
    
end
