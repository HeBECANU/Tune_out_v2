function varargout = stokes_intensity(varargin)
    %function varargout = stokes_intensity(varargin)
    %Return intensity of light described by Stokes vectors.
    %
    %  Syntax:
    %     p = stokes_intensity(V)
    %     P = stokes_intensity(C)
    %     [p,q,...] = stokes_intensity(V,W,...)
    %     {P,Q,...} = stokes_intensity(C,D,...)
    %
    %  Description:
    %     p = stokes_intensity(V) returns the light intensity for a Stokes vector V.
    %
    %     P = stokes_intensity(C), where C is a cell array of Stokes vectors,
    %     returns a numeric array of same size, size(P)==size(C).
    %
    %     [p,q,...] = stokes_intensity(V,W,...) returns multiple numeric
    %     values with the light intensity of Stokes vectors V,W,...
    %
    %     [P,Q,...] = stokes_intensity(C,D,...) with cell arrays C,D,... returns
    %     numeric arrays [T,U,...] with intensity values
    %
    %  Example:
    %     t = stokes_intensity(stokes_unpolarized()+stokes_cpleft())
    %     t = 2.0
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
            intensity = zeros(size(V));
            V_subs = cell(1,ndims(V));
            for vi=1:numel(V)
                [V_subs{:}] = ind2sub(size(V),vi);
                intensity(V_subs{:}) = s_intensity(V{V_subs{:}});
            end
        else
            intensity = s_intensity(V);
        end
        
        varargout{ni} = intensity;
        
    end
    
end

% helper function
function intensity = s_intensity(V)
    intensity = V(1);
end
