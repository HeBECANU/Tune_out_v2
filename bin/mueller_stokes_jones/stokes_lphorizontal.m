function V = stokes_lphorizontal(varargin)
    %function V = stokes_lphorizontal(varargin)
    %Return the Stokes vector for horizontal linearly polarized light.
    %
    %  Syntax:
    %     V = stokes_lphorizontal()
    %     V = stokes_lphorizontal(p)
    %     A = stokes_lphorizontal(P)
    %
    %  Description:
    %     V = stokes_lphorizontal() returns the Stokes vector
    %     for light with horizontal linear polarization and intensity 1.
    %
    %     V = stokes_lphorizontal(p) returns the Stokes vector
    %     for light with horizontal linear polarization and intensity p.
    %
    %     A = stokes_lphorizontal(P), where P is a m x n x ... array of
    %     intensity values, returns a cell array of
    %     Stokes vectors with size(A)==size(P).
    %
    %  Example:
    %     V = stokes_lphorizontal(0.4)
    %     V =
    %
    %        0.40000
    %        0.40000
    %        0.00000
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
    %     stokes_lpvertical
    %
    %  File information:
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
            V{V_subs{:}} = s_lphorizontal(intensity(V_subs{:}));
        end
        
    else
        
        V = s_lphorizontal(intensity);
        
    end
    
end

% helper function
function V = s_lphorizontal(intensity)
    V = [intensity;intensity;0;0];
end
