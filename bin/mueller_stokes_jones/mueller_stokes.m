function A = mueller_stokes(varargin)
    %function A = mueller_stokes(varargin)
    %Multiply Mueller matrices and Stokes vectors.
    %
    %  Syntax:
    %     M = mueller_stokes(M)
    %     A = mueller_stokes(C1, C2, ...)
    %
    %  Description:
    %     M = mueller_stokes(M)
    %
    %     A = mueller_stokes(C1, C2, ...)
    %
    %  Examples:
    %     M = mueller_stokes(mueller_linpolarizer(), mueller_rotate(mueller_linpolarizer(), 90, 'deg'))
    %     % crossed polarizers!
    %     M =
    %
    %         0     0     0     0
    %         0     0     0     0
    %         0     0     0     0
    %         0     0     0     0
    %
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
