function JM = jones_linretarder(varargin)
    %function JM = jones_linretarder(varargin)
    %Return the Jones matrix for a linear retarder with long axis rotation of 0 degrees.
    %
    %  Syntax:
    %     JM = jones_linretarder()
    %     JM = jones_linretarder(p)
    %     A = jones_linretarder(P)
    %     jones_linretarder(..., mode)
    %
    %  Description:
    %     JM = jones_linretarder() returns the Jones matrix
    %     for a linear retarder with zero phase delay.
    %
    %     JM = jones_linretarder(p) returns the Jones matrix
    %     for a linear retarder with phase delay p in radiant
    %     units, i.e. p is ranging between 0 and 2*pi().
    %
    %     A = jones_linretarder(P), where P is a m-by-n-by-... array of
    %     phase delay values, returns a cell array of
    %     linear retarder jones matrices with size(A)==size(P).
    %
    %     jones_linretarder(..., mode) where mode is a string defining
    %     the units for the phase delay: 'radiant' (default),
    %     'degree' (0..360) or 'wavelength' (0..1).
    %
    %  Example:
    %     JM = jones_linretarder(pi()/3)
    %     JM =
    %
    %        1.00000   0.00000   
    %        0.00000   0.50000 - 0.86603i
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
    %     jones_waveplate jones_circretarder
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    phase_defv = 0;
    
    if nargin<1
        phase = phase_defv;
    else
        phase = varargin{1};
    end
    
    [phase, was_cell] = c2n(phase, phase_defv);
    
    if nargin>=2 && ischar(varargin{end})
        if strncmpi(varargin{end},'deg',3)
            phase = phase*pi()/180.0;
        elseif strncmpi(varargin{end},'wav',3)
            phase = phase*2*pi();
        end
    end
    
    if (numel(phase) > 1) || was_cell
        
        JM = cell(size(phase));
        JM_subs = cell(1,ndims(JM));
        for jmi=1:numel(JM)
            [JM_subs{:}] = ind2sub(size(JM),jmi);
            JM{JM_subs{:}} = s_linretarder(phase(JM_subs{:}));
        end
        
    else
        
        JM = s_linretarder(phase);
        
    end
    
end

% helper function
function JM = s_linretarder(phase_in_pi_units)
    
    JM = zeros(2,2);
    JM(1,1) = 1;
    JM(2,2) = exp(-1i*phase_in_pi_units);
    
end
