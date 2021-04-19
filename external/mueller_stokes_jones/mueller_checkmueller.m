function varargout = mueller_checkmueller(varargin)
    %function varargout = mueller_checkmueller(varargin)
    %Check physical validity of Mueller matrix or matrices.
    %
    %  Syntax:
    %     t = mueller_checkmueller(M)
    %     T = mueller_checkmueller(C)
    %     [t,u,...] = mueller_checkmueller(M,N,...)
    %     {T,U,...} = mueller_checkmueller(C,D,...)
    %
    %  Description:
    %     t = mueller_checkmueller(M) returns true if the numeric
    %     matrix M is a physically valid Mueller matrix.
    %
    %     T = mueller_checkmueller(C), where C is a cell array of
    %     potential Mueller matrices, returns a boolean array
    %     of same size, size(T)==size(C).
    %
    %     [t,u,...] = mueller_checkmueller(M,N,...) returns multiple boolean
    %     values, depending on whether M,N,... are physically valid Mueller matrices
    %
    %     [T,U,...] = mueller_checkmueller(C,D,...) with cell arrays C,D,... returns
    %     boolean arrays [T,U,...] with corresponding physical Mueller validity
    %
    %  Example:
    %     t = mueller_checkmueller(mueller_linpolarizer())
    %     t = 1
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
    %     mueller_ismueller
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
            checkmueller = false(size(M));
            M_subs = cell(1,ndims(M));
            for mi=1:numel(M)
                [M_subs{:}] = ind2sub(size(M),mi);
                checkmueller(M_subs{:}) = s_checkmueller(M{M_subs{:}});
            end
        else
            checkmueller = s_checkmueller(M);
        end
        
        varargout{vi} = checkmueller;
        
    end
    
end

% helper function
function checkmueller = s_checkmueller(M)
    if ~isnumeric(M)
        checkmueller = false;
    elseif ~all(size(M)==[4,4])
        checkmueller = false;
        %    elseif det(M) == 0
        %        checkmueller = false;
    else
        checkmueller = true;
    end
end


