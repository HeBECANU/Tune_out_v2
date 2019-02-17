function JM = jones_rotate(varargin)
    %function JM = jones_rotate(varargin)
    %Return the Jones matrix for rotated Jones elements.
    %
    %  Syntax:
    %     JM = jones_rotate()
    %     JM = jones_rotate(JM)
    %     A = jones_rotate(C)
    %     JM = jones_rotate(JM,p)
    %     A = jones_rotate(C,p)
    %     A = jones_rotate(C,P)
    %     jones_rotate(..., mode)
    %
    %  Description:
    %     JM = jones_rotate() returns the matrix for the non-rotated unity Jones element.
    %
    %     JM = jones_rotate(JM) returns the matrix for the non-rotated
    %     element defined by JM.
    %
    %     A = jones_rotate(C), where C is a m-by-n-by-... cell array of
    %     Jones matrices returns the a cell array of Jones
    %     matrices of non-rotated elements with size(A)==size(C).
    %
    %     JM = jones_rotate(JM,p) returns the Jones matrix for the
    %     an element rotated by the angle p.
    %
    %     A = jones_rotate(C,p), where C is a m-by-n-by-... cell array of
    %     Jones matrices, returns the a cell array of Jones matrices with all
    %     elements are rotated by the angle p.
    %
    %     A = jones_rotate(C,P), where C is a m-by-n-by-... cell array of
    %     Jones matrices and P is a numeric array of same size, returns
    %     the a cell array of Jones matrices with all elements are rotated
    %     by the angle p.
    %
    %     jones_rotate(..., mode), where mode is a string defining
    %     the interpretation of the angle value: 'radiants' (default) or 'degree'.
    %
    %  Examples:
    %     JM = jones_rotate(jones_linpolarizer(),pi()/3)
    %     JM =
    %
    %        0.25000   0.43301
    %        0.43301   0.75000
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
    %     jones_rotator
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    angle_defv = 0;
    
    if nargin<1
        JM = jones_unity();
        return;
    elseif nargin<2
        JM = varargin{1};
        return;
    else
        C = varargin{1};
        angle = varargin{2};
    end
    
    [angle, angle_was_cell] = c2n(angle, angle_defv);
    
    if nargin>=3 && ischar(varargin{end})
        if strncmpi(varargin{end},'deg',3)
            angle = angle*pi()/180.0;
        end
    end
    
    if iscell(C) || (numel(angle) > 1) || angle_was_cell
        
        if ~iscell(C)
            C = {C};
        end
        
        % adjust dimensions, i.e. fill missing dimensions with 1
        sizeC = size(C);
        sizeangle = size(angle);
        maxdim = max(length(sizeC),length(sizeangle));
        if length(sizeC) < maxdim
            sizeC = [sizeC, ones(1,maxdim-length(sizeC))];
        end
        if length(sizeangle) < maxdim
            sizeangle = [sizeangle, ones(1,maxdim-length(sizeangle))];
        end
        
        % generate Jones matrices
        maxsize = max([sizeC;sizeangle]);
        JM = cell(maxsize);
        JM_subs = cell(1,ndims(JM));
        numelJM = numel(JM);
        
        % flatten C and angle arrays
        C = C(:);
        angle = angle(:);
        numelC = numel(C);
        numelangle = numel(angle);
        
        for jmi=1:numelJM
            [JM_subs{:}] = ind2sub(size(JM),jmi);
            JM{JM_subs{:}} = s_rotate(C{mod(jmi-1,numelC)+1}, angle(mod(jmi-1,numelangle)+1));
        end
        
    else
        
        JM = s_rotate(C, angle(1));
        
    end
    
end

% helper function
function JM = s_rotate(JM, angle_in_radiants)
    JM = s_rotator(-angle_in_radiants)*JM*s_rotator(angle_in_radiants);
end

% helper function
function JM = s_rotator(angle_in_radiants)
    
    JM = zeros(2,2);
    
    JM(1,1) = cos(angle_in_radiants);
    JM(1,2) = sin(angle_in_radiants);
    JM(2,1) = -sin(angle_in_radiants);
    JM(2,2) = cos(angle_in_radiants);
    
end
