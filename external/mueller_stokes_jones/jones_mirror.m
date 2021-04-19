function JM = jones_mirror(varargin)
    %function JM = jones_mirror(varargin)
    %Return the Jones matrix for an ideal mirror.
    %
    %  Syntax:
    %     JM = jones_mirror()
    %     A = jones_mirror([m, n, ...])
    %     A = jones_mirror(C)
    %
    %  Description:
    %     JM = jones_mirror() returns the Jones matrix
    %     for an ideal mirror element.
    %
    %     A = jones_mirror([m, n, ...]) returns a cell array of
    %     mirror Jones matrices with size(A)==[m, n, ...].
    %
    %     A = jones_mirror(C) returns a cell array of
    %     mirror Jones matrices with size(A)==size(C).
    %
    %  Example:
    %     JM = jones_mirror()
    %     JM =
    %
    %       -1.00000   0.00000
    %        0.00000  -1.00000
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
    %     jones_unity
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
        
        JM = cell(sc);
        [JM{:}] = deal(s_mirror());
        
    else
        
        JM = s_mirror();
        
    end
    
end

% helper function
function JM = s_mirror()
    
    JM = zeros(2,2);
    JM(1,1) = -1;
    JM(2,2) = -1;
    
end

