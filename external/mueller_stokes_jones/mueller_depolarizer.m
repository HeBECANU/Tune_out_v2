function M = mueller_depolarizer(varargin)
    %function M = mueller_depolarizer(varargin)
    %Return the Mueller matrix for a (partial) depolarizer.
    %
    %  Syntax:
    %     M = mueller_depolarizer()
    %     M = mueller_depolarizer(p)
    %     A = mueller_depolarizer(P)
    %
    %  Description:
    %     M = mueller_depolarizer() returns the Mueller matrix
    %     for a depolarizer with zero depolarization.
    %
    %     M = mueller_depolarizer(p) returns the Mueller matrix
    %     for a depolarizer with depolarization p ranging between 0 and 1.
    %
    %     A = mueller_depolarizer(P), where P is a m-by-n-by-... array of
    %     depolarization values, returns a cell array of
    %     depolarizer Mueller matrices with size(A)==size(P).
    %
    %  Example:
    %     M = mueller_depolarizer(0.5)
    %     M =
    %
    %        1.00000   0.00000   0.00000   0.00000
    %        0.00000   0.50000   0.00000   0.00000
    %        0.00000   0.00000   0.50000   0.00000
    %        0.00000   0.00000   0.00000   0.50000
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
    %     mueller_linpolarizer
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    depolarization_defv = 0;
    
    if nargin<1
        depolarization = depolarization_defv;
    else
        depolarization = varargin{1};
    end
    
    [depolarization, was_cell] = c2n(depolarization, depolarization_defv);
    
    if (numel(depolarization) > 1) || was_cell
        
        M = cell(size(depolarization));
        M_subs = cell(1,ndims(M));
        for mi=1:numel(M)
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_depolarizer(depolarization(M_subs{:}));
        end
        
    else
        
        M = s_depolarizer(depolarization);
        
    end
    
end

% helper function
function M = s_depolarizer(depolarization)
    
    M = zeros(4,4);
    M(1,1) = 1;
    M(2,2) = depolarization;
    M(3,3) = depolarization;
    M(4,4) = depolarization;
    
end
