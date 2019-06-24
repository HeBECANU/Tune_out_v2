function V = jones_lpminus45(varargin)
    %function V = jones_lpminus45(varargin)
    %Return the Jones vector for linearly polarized light at -45 degrees.
    %
    %  Syntax:
    %     V = jones_lpminus45()
    %     V = jones_lpminus45(p)
    %     A = jones_lpminus45(P)
    %
    %  Description:
    %     V = jones_lpminus45() returns the Jones vector
    %     for light with linear polarization at -45 degrees and amplitude 1.
    %
    %     V = jones_lpminus45(p) returns the Jones vector
    %     for light with linear polarization at -45 degrees and amplitude p.
    %
    %     A = jones_lpminus45(P), where P is a m x n x ... array of
    %     amplitude values, returns a cell array of
    %     Jones vectors with size(A)==size(P).
    %
    %  Example:
    %     V = jones_lpminus45(0.4)
    %     V =
    %
    %        0.40000
    %       -0.40000
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
    %     jones_lpplus45 jones_lphorizontal jones_lpvertical
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    amplitude_defv = 1/sqrt(2);
    
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
            V{V_subs{:}} = s_lpminus45(amplitude(V_subs{:}));
        end
        
    else
        
        V = s_lpminus45(amplitude);
        
    end
    
end

% helper function
function V = s_lpminus45(amplitude)
    V = amplitude*[1; -1];
end


