function M = mueller_absorber(varargin)
    %function M = mueller_absorber(varargin)
    %Return the Mueller matrix for a (partial) absorber.
    %
    %  Syntax:
    %     M = mueller_absorber()
    %     M = mueller_absorber(p)
    %     A = mueller_absorber(P)
    %
    %  Description:
    %     M = mueller_absorber() returns the Mueller matrix
    %     for a absorber with zero absorption.
    %
    %     M = mueller_absorber(p) returns the Mueller matrix
    %     for a absorber with absorption p ranging between 0 and 1.
    %
    %     A = mueller_absorber(P), where P is a m x n x ... array of
    %     absorption values, returns a cell array of
    %     absorber Mueller matrices with size(A)==size(P).
    %
    %  Example:
    %     M = mueller_absorber(0.4)
    %     M =
    %
    %        0.60000   0.00000   0.00000   0.00000
    %        0.00000   0.60000   0.00000   0.00000
    %        0.00000   0.00000   0.60000   0.00000
    %        0.00000   0.00000   0.00000   0.60000
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
    %     mueller_lindiattenuator mueller_circdiattenuator
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    absorption_defv = 0;
    
    if nargin<1
        absorption = absorption_defv;
    else
        absorption = varargin{1};
    end
    
    [absorption, was_cell] = c2n(absorption, absorption_defv);
    
    if (numel(absorption) > 1) || was_cell
        
        M = cell(size(absorption));
        M_subs = cell(1,ndims(M));
        for mi=1:numel(M)
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_absorber(absorption(M_subs{:}));
        end
        
    else
        
        M = s_absorber(absorption);
        
    end
    
end

% helper function
function M = s_absorber(absorption)
    M = eye(4)*(1-absorption);
end
