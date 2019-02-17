function M = mueller_waveplate(varargin)
    %function M = mueller_waveplate(varargin)
    %Return the Mueller matrix for a linear wave plate with a phase delay given in wavelength units and long axis rotation of 0 degrees.
    %
    %  Syntax:
    %     M = mueller_waveplate()
    %     M = mueller_waveplate(p)
    %     A = mueller_waveplate(P)
    %
    %  Description:
    %     M = mueller_waveplate() returns the Mueller matrix
    %     for a wave plate with zero phase delay.
    %
    %     M = mueller_waveplate(p) returns the Mueller matrix
    %     for a wave plate with phase delay p in wavelength
    %     units, i.e. p is ranging between 0 and 1.
    %
    %     A = mueller_waveplate(P), where P is a m-by-n-by-... array of
    %     phase delay values, returns a cell array of
    %     wave plate Mueller matrices with size(A)==size(P).
    %
    %  Example:
    %     M = mueller_waveplate(1/6) % Mueller matrix for sixth wave plate
    %     M =
    %
    %        1.00000   0.00000   0.00000   0.00000
    %        0.00000   1.00000   0.00000   0.00000
    %        0.00000   0.00000   0.50000  -0.86603
    %        0.00000   0.00000   0.86603   0.50000
    %
    %  References:
    %     [1] E. Collett, Field Guide to Polarization,
    %         SPIE Field Guides vol. FG05, SPIE (2005). ISBN 0-8194-5868-6.
    %     [2] R. A. Chipman, "Polarimetry," chapter 22 in Handbook of Optics II,
    %         2nd Ed, M. Bass, editor in chief (McGraw-Hill, New York, 1995)
    %     [3] "Mueller calculus", http://en.wikipedia.org/wiki/Mueller_calculus,
    %         last retrieved on Dec 17, 2013.
    %
    %  See also:
    %     mueller_linretarder  mueller_circretarder
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
        
        M = cell(size(phase_in_lambda_units));
        M_subs = cell(1,ndims(M));
        for mi=1:numel(M)
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_waveplate(phase_in_lambda_units(M_subs{:}));
        end
        
    else
        
        M = s_waveplate(phase_in_lambda_units);
        
    end
    
end

% helper function
function M = s_waveplate(phase_in_lambda_units)
    
    M = zeros(4,4);
    
    M(1,1) = 1;
    M(2,2) = 1;
    M(3,3) = cos(phase_in_lambda_units*2*pi());
    M(3,4) = -sin(phase_in_lambda_units*2*pi());
    M(4,3) = -M(3,4);
    M(4,4) = M(3,3);
    
end
