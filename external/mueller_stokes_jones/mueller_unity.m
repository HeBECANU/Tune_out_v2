function M = mueller_unity(varargin)
    %function M = mueller_unity(varargin)
    %Return the unity Mueller matrix.
    %
    %  Syntax:
    %     M = mueller_unity()
    %     A = mueller_unity([m, n, ...])
    %     A = mueller_unity(C)
    %
    %  Description:
    %     M = mueller_unity() returns the unity Mueller matrix,
    %     representing a non-polarizing optical element.
    %
    %     A = mueller_unity([m, n, ...]) returns a cell array of
    %     unity Mueller matrices with size(A)==[m, n, ...].
    %
    %     A = mueller_unity(C) returns a cell array of
    %     unity Mueller matrices with size(A)==size(C).
    %
    %  Example:
    %     M = mueller_unity()
    %     M =
    %
    %        1.00000   0.00000   0.00000   0.00000
    %        0.00000   1.00000   0.00000   0.00000
    %        0.00000   0.00000   1.00000   0.00000
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
    %     mueller_mirror
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
        [M{:}] = deal(s_unity());
        
    else
        
        M = s_unity();
        
    end
    
end

% helper function
function M = s_unity()
    
    M = eye(4);
    
end
