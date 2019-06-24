function V = jones_lpvertical(varargin)
    %function V = jones_lpvertical(varargin)
    %Return the Jones vector for vertical linearly polarized light.
    %
    %  Syntax:
    %     V = jones_lpvertical()
    %     V = jones_lpvertical(p)
    %     A = jones_lpvertical(P)
    %
    %  Description:
    %     V = jones_lpvertical() returns the Jones vector
    %     for light with vertical linear polarization and amplitude 1.
    %
    %     V = jones_lpvertical(p) returns the Jones vector
    %     for light with vertical linear polarization and amplitude p.
    %
    %     A = jones_lpvertical(P), where P is a m x n x ... array of
    %     amplitude values, returns a cell array of
    %     Jones vectors with size(A)==size(P).
    %
    %  Example:
    %     V = jones_lpvertical(0.4)
    %     V =
    %
    %        0.00000
    %        0.40000
    %
    %  References:
    %     [1] E. Collett, Field Guide to Polarization,
    %         SPIE Field Guides vol. FG05, SPIE (2005). ISBN 0-8194-5868-6.
    %     [2] R. A. Chipman, "Polarimetry," chapter 22 in Handbook of Optics II,
    %         2nd Ed, M. Bass, editor in chief (McGraw-Hill, New York, 1995)
    %     [3] "Jones calculus", http://en.wikipedia.org/wiki/Jones_calculus,
    %         last retrieved on Jan 13, 2014.
    %
    %  See also:
    %     jones_lphorizontal
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    amplitude_defv = 1;
    
    if nargin<1
        amplitude = amplitude_defv;
    else
        amplitude = varargin{1};
    end
    
    [amplitude, was_cell] = c2n(amplitude, amplitude_defv);
    
    if (numel(amplitude) > 1) || was_cell
        
        V = cell(size(amplitude));
        V_subs = cell(1,ndims(V));
        for Vi=1:numel(V)
            [V_subs{:}] = ind2sub(size(V),Vi);
            V{V_subs{:}} = s_lpvertical(amplitude(V_subs{:}));
        end
        
    else
        
        V = s_lpvertical(amplitude);
        
    end
    
end

% helper function
function V = s_lpvertical(amplitude)
    V = [0; amplitude];
end
