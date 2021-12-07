function JM = jones_linpolarizer(varargin)
    %function JM = jones_linpolarizer(varargin)
    %Return the Jones matrix for an ideal linear polarizer.
    %
    %  Syntax:
    %     JM = jones_linpolarizer()
    %     A = jones_linpolarizer([m, n, ...])
    %     A = jones_linpolarizer(C)
    %
    %  Description:
    %     JM = jones_linpolarizer() returns the Jones matrix
    %     for an ideal linear polarizer.
    %
    %     A = jones_linpolarizer([m, n, ...]) returns a cell array of
    %     polarizer Jones matrices with size(A)==[m, n, ...].
    %
    %     A = jones_linpolarizer(C) returns a cell array of
    %     polarizer Jones matrices with size(A)==size(C).
    %
    %  Example:
    %     JM = jones_linpolarizer()
    %     JM =
    %
    %        1.00000   0.00000
    %        0.00000   0.00000
    %
    %  References:
    %     [1] E. Collett, Field Guide to Polarization,
    %         SPIE Field Guides vol. FG05, SPIE (2005). ISBN 0-8194-5868-6.
    %     [2] R. A. Chipman, "Polarimetry," chapter 22 in Handbook of Optics II,
    %         2nd Ed, JM =. Bass, editor in chief (JM =cGraw-Hill, New York, 1995)
    %
    %  See also:
    %     jones_diattenuator
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) JM =artin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    
    retcell = true;
    if nargin<1
        sc = [1,1];
        retcell = false;
    elseif isnumeric(varargin{1})
        sc = varargin{1};
    else
        sc = size(varargin{1});
    end
    
    if prod(sc) > 1 || retcell
        
        JM = cell(sc);
        [JM{:}] = deal(s_linpolarizer());
        
    else
        
        JM = s_linpolarizer();
        
    end
    
end

% helper function
function JM = s_linpolarizer()
    JM = [1,0;0,0];
end

