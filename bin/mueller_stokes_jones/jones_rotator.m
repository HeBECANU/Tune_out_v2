function JM = jones_rotator(varargin)
    %function JM = jones_rotator(varargin)
    %Return the Jones matrix for a system rotator.
    %
    %     JM = jones_rotator()
    %     JM = jones_rotator(p)
    %     A = jones_rotator(P)
    %     jones_rotator(..., mode)
    %
    %  Description:
    %     JM = jones_rotator() returns the Jones matrix
    %     for a rotator with zero rotation.
    %
    %     JM = jones_rotator(p) returns the Jones matrix
    %     for a rotator at angle p.
    %
    %     A = jones_rotator(P), where P is a m-by-n-by-... array of
    %     angle values, returns a cell array of
    %     rotator Jones matrices with size(A)==size(P).
    %
    %     jones_rotator(..., mode), where mode is a string defining
    %     the units for the angle: 'radiant' (default) or 'degree' (0..360).
    %
    %  Example:
    %     JM = jones_rotator(pi()/3)
    %     JM =
    %
    %        0.50000   0.86603
    %       -0.86603   0.50000
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
    %     jones_rotate
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    angle_defv = 0;
    
    if nargin<1
        angle = angle_defv;
    else
        angle = varargin{1};
    end
    
    [angle, was_cell] = c2n(angle, angle_defv);
    
    if nargin>=2 && ischar(varargin{end})
        if strncmpi(varargin{end},'deg',3)
            angle = angle*pi()/180.0;
        end
    end
    
    if (numel(angle) > 1) || was_cell
        
        JM = cell(size(angle));
        M_subs = cell(1,ndims(M));
        for mi=1:numel(M)
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_rotator(angle(M_subs{:}));
        end
        
    else
        
        JM = s_rotator(angle);
        
    end
    
end

% helper function
function JM = s_rotator(angle_in_radiants)
    
    JM = zeros(2,2);
    
    JM(1,1) = cos(angle_in_radiants);
    JM(1,2) = sin(angle_in_radiants);
    JM(2,1) = -sin(angle_in_radiants);
    JM(2,2) = cos(angle_in_radiants);
    
end
