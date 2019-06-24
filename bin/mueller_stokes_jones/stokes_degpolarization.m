function varargout = stokes_degpolarization(varargin)
    %function varargout = stokes_degpolarization(varargin)
    %Return degree of polarization of light described by Stokes vectors.
    %
    %  Syntax:
    %     p = stokes_degpolarization(V)
    %     P = stokes_degpolarization(C)
    %     [p,q,...] = stokes_degpolarization(V,W,...)
    %     {P,Q,...} = stokes_degpolarization(C,D,...)
    %
    %  Description:
    %     p = stokes_degpolarization(V) returns the degree of
    %     polarization for a Stokes vector V.
    %
    %     P = stokes_degpolarization(C), where C is a cell array of
    %     potential Stokes vectors, returns a numeric array
    %     of same size, size(P)==size(C).
    %
    %     [p,q,...] = stokes_degpolarization(V,W,...) returns multiple numeric
    %     values with the degree of polarization of Stokes vectors V,W,...
    %
    %     [P,Q,...] = stokes_degpolarization(C,D,...) with cell arrays C,D,... returns
    %     numeric arrays [T,U,...] with corresponding degrees of polarization
    %
    %  Example:
    %     t = stokes_degpolarization(stokes_unpolarized()+stokes_cpleft())
    %     t = 0.5
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
            degpolarization = zeros(size(V));
            V_subs = cell(1,ndims(V));
            for vi=1:numel(V)
                [V_subs{:}] = ind2sub(size(V),vi);
                degpolarization(V_subs{:}) = s_degpolarization(V{V_subs{:}});
            end
        else
            degpolarization = s_degpolarization(V);
        end
        
        varargout{ni} = degpolarization;
        
    end
    
end

% helper function
function degpolarization = s_degpolarization(V)
    S13 = V(2:4);
    degpolarization = sqrt(S13'*S13)/V(1);
end
