function M = mueller_rotator(varargin)
    %function M = mueller_rotator(varargin)
    %Return the Mueller matrix for a system rotator.
    %
    %     M = mueller_rotator()
    %     M = mueller_rotator(p)
    %     A = mueller_rotator(P)
    %     mueller_rotator(..., mode)
    %
    %  Description:
    %     M = mueller_rotator() returns the Mueller matrix
    %     for a rotator with zero rotation.
    %
    %     M = mueller_rotator(p) returns the Mueller matrix
    %     for a rotator at angle p.
    %
    %     A = mueller_rotator(P), where P is a m-by-n-by-... array of
    %     angle values, returns a cell array of
    %     rotator Mueller matrices with size(A)==size(P).
    %
    %     mueller_rotator(..., mode), where mode is a string defining
    %     the units for the angle: 'radiant' (default) or 'degree' (0..360).
    %
    %  Example:
    %     M = mueller_rotator(pi()/3)
    %     M =
    %
    %        1.00000   0.00000   0.00000   0.00000
    %        0.00000  -0.50000   0.86603   0.00000
    %        0.00000  -0.86603  -0.50000   0.00000
    %        0.00000   0.00000   0.00000   1.00000
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
    %     mueller_rotate
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
        
        M = cell(size(angle));
        M_subs = cell(1,ndims(M));
        for mi=1:numel(M)
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_rotator(angle(M_subs{:}));
        end
        
    else
        
        M = s_rotator(angle);
        
    end
    
end

% helper function
function M = s_rotator(angle_in_radiants)
    
    M = zeros(4,4);
    
    % short cut to avoid "*2" in each line
    angle_in_radiants = angle_in_radiants*2;
    
    M(1,1) = 1;
    M(2,2) = cos(angle_in_radiants);
    M(2,3) = sin(angle_in_radiants);
    M(3,2) = -sin(angle_in_radiants);
    M(3,3) = cos(angle_in_radiants);
    M(4,4) = 1;
    
end
