function JM = jones_unity(varargin)
    %function JM = jones_unity(varargin)
    %Return the unity Jones matrix.
    %
    %  Syntax:
    %     JM = jones_unity()
    %     A = jones_unity([m, n, ...])
    %     A = jones_unity(C)
    %
    %  Description:
    %     JM = jones_unity() returns the unity Jones matrix,
    %     representing a non-polarizing optical element.
    %
    %     A = jones_unity([m, n, ...]) returns a cell array of
    %     unity Jones matrices with size(A)==[m, n, ...].
    %
    %     A = jones_unity(C) returns a cell array of
    %     unity Jones matrices with size(A)==size(C).
    %
    %  Example:
    %     JM = jones_unity()
    %     JM =
    %
    %          1   0
    %          0   1
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
    %     jones_mirror
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
        [JM{:}] = deal(s_unity());
        
    else
        
        JM = s_unity();
        
    end
    
end

% helper function
function JM = s_unity()
    JM = eye(2);
end
