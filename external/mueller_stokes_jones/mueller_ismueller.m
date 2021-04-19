function varargout = mueller_ismueller(varargin)
    %function varargout = mueller_ismueller(varargin)
    %Check validity of Mueller matrix or matrices.
    %
    %  Syntax:
    %     t = mueller_ismueller(M)
    %     T = mueller_ismueller(C)
    %     [t,u,...] = mueller_ismueller(M,N,...)
    %     {T,U,...} = mueller_ismueller(C,D,...)
    %
    %  Description:
    %     t = mueller_ismueller(M) returns true if the numeric
    %     matrix M is a valid Mueller matrix.
    %
    %     T = mueller_ismueller(C), where C is a cell array of
    %     potential Mueller matrices, returns a boolean array
    %     of same size, size(T)==size(C).
    %
    %     [t,u,...] = mueller_ismueller(M,N,...) returns multiple boolean
    %     values, depending on whether M,N,... are valid Mueller matrices
    %
    %     [T,U,...] = mueller_ismueller(C,D,...) with cell arrays C,D,... returns
    %     boolean arrays [T,U,...] with corresponding Mueller validity
    %
    %  Example:
    %     t = mueller_ismueller(mueller_linpolarizer())
    %     t = 1
    %
    %     t = mueller_ismueller(zeros(5,5))
    %     t = 0
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
    %     mueller_checkmueller
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    if nargin==0
        varargout = {};
        return;
    end
    
    % loop over parameters
    for vi=1:nargin
        
        M = varargin{vi};
        if iscell(M)
            ismueller = false(size(M));
            M_subs = cell(1,ndims(M));
            for mi=1:numel(M)
                [M_subs{:}] = ind2sub(size(M),mi);
                ismueller(M_subs{:}) = s_ismueller(M{M_subs{:}});
            end
        else
            ismueller = s_ismueller(M);
        end
        
        varargout{vi} = ismueller;
        
    end
    
end

% helper function
function isMueller = s_ismueller(M)
    if ~isnumeric(M)
        isMueller = false;
    elseif ~all(size(M)==[4,4])
        isMueller = false;
    else
        isMueller = true;
    end
    
end


