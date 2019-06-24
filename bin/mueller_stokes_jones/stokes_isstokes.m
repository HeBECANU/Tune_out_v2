function varargout = stokes_isstokes(varargin)
    %function varargout = stokes_isstokes(varargin)
    %Check validity of Stokes vectors.
    %
    %  Syntax:
    %     t = stokes_isstokes(V)
    %     T = stokes_isstokes(C)
    %     [t,u,...] = stokes_isstokes(V,W,...)
    %     {T,U,...} = stokes_isstokes(C,D,...)
    %
    %  Description:
    %     t = stokes_isstokes(V) returns true if the numeric
    %     vector V is a valid Stokes vector.
    %
    %     T = stokes_isstokes(C), where C is a cell array of
    %     potential Stokes vectors, returns a boolean array
    %     of same size, size(T)==size(C).
    %
    %     [t,u,...] = stokes_isstokes(V,W,...) returns multiple boolean
    %     values, depending on whether V,W,... are valid Stokes vectors
    %
    %     [T,U,...] = stokes_isstokes(C,D,...) with cell arrays C,D,... returns
    %     boolean arrays [T,U,...] with corresponding Stokes validity
    %
    %  Example:
    %     t = stokes_isstokes(stokes_cpleft())
    %     t = 1
    %
    %     t = stokes_isstokes(zeros(5,1))
    %     t = 0
    %
    %  References:
    %     [1] E. Collett, Field Guide to Polarization,
    %         SPIE Field Guides vol. FG05, SPIE (2005). ISBN 0-8194-5868-6.
    %     [2] R. A. Chipman, "Polarimetry," chapter 22 in Handbook of Optics II,
    %         2nd Ed, M. Bass, editor in chief (McGraw-Hill, New York, 1995)
    %     [3] "Stokes parameters", http://en.wikipedia.org/wiki/Stokes_parameters,
    %         last retrieved on Dec 17, 2013.
    %
    %  See also:
    %     stokes_unpolarized
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
    for ni=1:nargin
        
        V = varargin{ni};
        if iscell(V)
            isstokes = false(size(V));
            V_subs = cell(1,ndims(V));
            for vi=1:numel(V)
                [V_subs{:}] = ind2sub(size(V),vi);
                isstokes(V_subs{:}) = s_isstokes(V{V_subs{:}});
            end
        else
            isstokes = s_isstokes(V);
        end
        
        varargout{ni} = isstokes;
        
    end
    
end

% helper function
function isstokes = s_isstokes(V)
    if ~isnumeric(V)
        isstokes = false;
    elseif ~all(size(V)==[4,1])
        isstokes = false;
    else
        isstokes = true;
    end
end


