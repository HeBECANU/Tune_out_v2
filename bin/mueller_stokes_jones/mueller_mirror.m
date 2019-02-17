function M = mueller_mirror(varargin)
    %function M = mueller_mirror(varargin)
    %Return the Mueller matrix for an ideal mirror.
    %
    %  Syntax:
    %     M = mueller_mirror()
    %     A = mueller_mirror([m, n, ...])
    %     A = mueller_mirror(C)
    %
    %  Description:
    %     M = mueller_mirror() returns the Mueller matrix
    %     for an ideal mirror element.
    %
    %     A = mueller_mirror([m, n, ...]) returns a cell array of
    %     mirror Mueller matrices with size(A)==[m, n, ...].
    %
    %     A = mueller_mirror(C) returns a cell array of
    %     mirror Mueller matrices with size(A)==size(C).
    %
    %  Example:
    %     M = mueller_mirror()
    %     M =
    %
    %        1.00000   0.00000   0.00000   0.00000
    %        0.00000   1.00000   0.00000   0.00000
    %        0.00000   0.00000  -1.00000   0.00000
    %        0.00000   0.00000   0.00000  -1.00000
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
    %     mueller_unity
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    retcell = true;
    if nargin<1
        sc = [1,1];
        retcell = false;
    elseif isnumeric(varargin{1})
        sc = varargin{1};
    else
        sc = size(varargin{1});
    end
    
    if prod(sc) > 1 || retcell
        
        M = cell(sc);
        [M{:}] = deal(s_mirror());
        
    else
        
        M = s_mirror();
        
    end
    
end

% helper function
function M = s_mirror()
    
    M = zeros(4,4);
    M(1,1) = 1;
    M(2,2) = 1;
    M(3,3) = -1;
    M(4,4) = -1;
    
end
