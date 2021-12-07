function V = stokes_lpplus45(varargin)
    %function V = stokes_lpplus45(varargin)
    %Return the Stokes vector for linearly polarized light at +45 degrees.
    %
    %  Syntax:
    %     V = stokes_lpplus45()
    %     V = stokes_lpplus45(p)
    %     A = stokes_lpplus45(P)
    %
    %  Description:
    %     V = stokes_lpplus45() returns the Stokes vector
    %     for light with linear polarization at +45 degrees and intensity 1.
    %
    %     V = stokes_lpplus45(p) returns the Stokes vector
    %     for light with linear polarization at +45 degrees and intensity p.
    %
    %     A = stokes_lpplus45(P), where P is a m x n x ... array of
    %     intensity values, returns a cell array of
    %     Stokes vectors with size(A)==size(P).
    %
    %  Example:
    %     V = stokes_lpplus45(0.4)
    %     V =
    %
    %        0.40000
    %        0.00000
    %        0.40000
    %        0.00000
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
    %     stokes_lpminus45 stokes_lphorizontal stokes_lpvertical
    %
    %  File information:
    %     Part of package stokes
    %
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    intensity_defv = 1;
    
    if nargin<1
        intensity = intensity_defv;
    else
        intensity = varargin{1};
    end
    
    [intensity, was_cell] = c2n(intensity, intensity_defv);
    
    if (numel(intensity) > 1) || was_cell
        
        V = cell(size(intensity));
        V_subs = cell(1,ndims(V));
        for Vi=1:numel(V)
            [V_subs{:}] = ind2sub(size(V),Vi);
            V{V_subs{:}} = s_lpplus45(intensity(V_subs{:}));
        end
        
    else
        
        V = s_lpplus45(intensity);
        
    end
    
end

% helper function
function V = s_lpplus45(intensity)
    V = [intensity;0;intensity;0];
end

