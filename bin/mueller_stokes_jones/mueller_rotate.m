function M = mueller_rotate(varargin)
    %function M = mueller_rotate(varargin)
    %Return the Mueller matrix for rotated Mueller elements.
    %
    %  Syntax:
    %     M = mueller_rotate()
    %     M = mueller_rotate(M)
    %     A = mueller_rotate(C)
    %     M = mueller_rotate(M,p)
    %     A = mueller_rotate(M,P)
    %     A = mueller_rotate(C,p)
    %     A = mueller_rotate(C,P)
    %     mueller_rotate(..., mode)
    %
    %  Description:
    %     M = mueller_rotate() returns the matrix for the non-rotated unity Mueller element.
    %
    %     M = mueller_rotate(M) returns the matrix for the non-rotated
    %     element defined by M.
    %
    %     A = mueller_rotate(C), where C is a m-by-n-by-... cell array of
    %     Mueller matrices returns the a cell array of Mueller
    %     matrices of non-rotated elements with size(A)==size(C).
    %
    %     M = mueller_rotate(M,p) returns the Mueller matrix for the
    %     an element rotated by the angle p.
    %
    %     A = mueller_rotate(M,P) where P is a m-by-n-by-... numeric array of
    %     angle values, returns a cell array of Mueller matrices where M
    %     is rotated by the angle values in P, size(A)==size(P).
    %
    %     A = mueller_rotate(C,p), where C is a m-by-n-by-... cell array of
    %     Mueller matrices, returns a cell array of Mueller matrices with all
    %     elements are rotated by the angle p.
    %
    %     A = mueller_rotate(C,P), where C is a m-by-n-by-... cell array of
    %     Mueller matrices and P is a numeric array of same size, returns
    %     the a cell array of Mueller matrices with all elements are rotated
    %     by the angle p.
    %
    %     mueller_rotate(..., mode), where mode is a string defining
    %     the interpretation of the angle value: 'radiants' (default) or 'degree'.
    %
    %  Examples:
    %     M = mueller_rotate(mueller_linpolarizer(),pi()/3)
    %     M =
    %
    %        0.5000   -0.2500    0.4330         0
    %       -0.2500    0.1250   -0.2165         0
    %        0.4330   -0.2165    0.3750         0
    %             0         0         0         0
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
    %     mueller_rotator
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
        M = mueller_unity();
        return;
    elseif nargin<2
        M = varargin{1};
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
        
        % generate Mueller matrices
        maxsize = max([sizeC;sizeangle]);
        M = cell(maxsize);
        M_subs = cell(1,ndims(M));
        numelM = numel(M);
        
        % flatten C and angle arrays
        C = C(:);
        angle = angle(:);
        numelC = numel(C);
        numelangle = numel(angle);
        
        for mi=1:numelM
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_rotate(C{mod(mi-1,numelC)+1}, angle(mod(mi-1,numelangle)+1));
        end
        
    else
        
        M = s_rotate(C, angle(1));
        
    end
    
end

% helper function
function M = s_rotate(M, angle_in_radiants)
    M = s_rotator(-angle_in_radiants)*M*s_rotator(angle_in_radiants);
end

% helper function
function M = s_rotator(angle_in_radiants)
    
    M = zeros(4,4);
    
    % short cut to avoid "*2" in each line
    angle_in_radiants = angle_in_radiants*2;
    
    M(1,1) = 1;
    M(2,2) = cos(angle_in_radiants);
    M(2,3) = sin(angle_in_radiants);
    M(3,2) = -sin(angle_in_radiants);
    M(3,3) = cos(angle_in_radiants);
    M(4,4) = 1;
end

