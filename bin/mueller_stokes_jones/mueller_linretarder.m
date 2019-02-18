function M = mueller_linretarder(varargin)
    %function M = mueller_linretarder(varargin)
    %Return the Mueller matrix for a linear retarder with long axis rotation of 0 degrees.
    %
    %  Syntax:
    %     M = mueller_linretarder()
    %     M = mueller_linretarder(p)
    %     A = mueller_linretarder(P)
    %     mueller_linretarder(..., mode)
    %
    %  Description:
    %     M = mueller_linretarder() returns the Mueller matrix
    %     for a linear retarder with zero phase delay.
    %
    %     M = mueller_linretarder(p) returns the Mueller matrix
    %     for a linear retarder with phase delay p in radiant
    %     units, i.e. p is ranging between 0 and 2*pi().
    %
    %     A = mueller_linretarder(P), where P is a m-by-n-by-... array of
    %     phase delay values, returns a cell array of
    %     linear retarder Mueller matrices with size(A)==size(P).
    %
    %     mueller_linretarder(..., mode) where mode is a string defining
    %     the units for the phase delay: 'radiant' (default),
    %     'degree' (0..360) or 'wavelength' (0..1).
    %
    %  Example:
    %     M = mueller_linretarder(pi()/3)
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
    %     mueller_waveplate mueller_circretarder
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
        
        M = cell(size(phase));
        M_subs = cell(1,ndims(M));
        for mi=1:numel(M)
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_linretarder(phase(M_subs{:}));
        end
        
    else
        
        M = s_linretarder(phase);
        
    end
    
end

% helper function
function M = s_linretarder(phase_in_pi_units)
    
    M = zeros(4,4);
    
    M(1,1) = 1;
    M(2,2) = 1;
    M(3,3) = cos(phase_in_pi_units);
    M(3,4) = -sin(phase_in_pi_units);
    M(4,3) = -M(3,4);
    M(4,4) = M(3,3);
    
end
