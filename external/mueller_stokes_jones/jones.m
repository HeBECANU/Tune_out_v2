function A = jones(varargin)
    %function A = jones(varargin)
    %Multiply Jones matrices and vectors.
    %
    %  Syntax:
    %     JM = jones(JM)
    %     A = jones(C1, C2, ...)
    %
    %  Description:
    %     M = jones(M)
    %
    %     A = jones(C1, C2, ...)
    %
    %  Examples:
    %     JM = jones(jones_linpolarizer(),jones_rotate(jones_linpolarizer(),90,'deg'))
    %     % crossed polarizers
    %     JM =
    %
    %       1.0e-016 *
    % 
    %         0.0000    0.6123
    %              0         0
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
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    %
    
    if nargin<1
        A = [];
        return;
    end
    
    A = varargin{1};
    
    for vi=2:nargin
        A = cellfunc(@(M1,M2)(M1*M2), A, varargin{vi});
    end
    
end
