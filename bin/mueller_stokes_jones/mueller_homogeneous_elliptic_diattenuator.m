function M = mueller_homogeneous_elliptic_diattenuator(varargin)
    %function M = mueller_homogeneous_elliptic_diattenuator(varargin)
    %Return the Mueller matrix for a homogeneous elliptic diattenuator (see refs 3 and 4).
    %
    %  Syntax:
    %     M = mueller_homogeneous_elliptic_diattenuator(t0, d, azimuth, ellipticity)
    %     A = mueller_homogeneous_elliptic_diattenuator(T0, D, Azimuth, Ellipticity)
    %     mueller_homogeneous_elliptic_diattenuator(..., mode)
    %
    %  Description:
    %     M = mueller_homogeneous_elliptic_diattenuator(t0, d, azimuth, ellipticity)
    %     returns the Mueller matrix for a homogeneous elliptic diattenuator with
    %     total transmission T0 (default: 1), diattenuation d (default: 0), and azimuth
    %     (default: 0) and ellipticity (default: 0) describe the two orthogonal polarization eigenstates.
    %
    %     A = mueller_homogeneous_elliptic_diattenuator(T0, D, Azimuth, Ellipticity)
    %     returns a cell array of Mueller matrices for homogeneous elliptic diattenuators if
    %     any of the four numeric parameters is a numeric array. Then, scalar parameters are
    %     expanded to uniform arrays and the size of A is given by
    %     min(size(T0),size(D),size(Azimuth),size(Ellipticity)).
    %
    %     mueller_homogeneous_elliptic_diattenuator(..., mode), where mode is a string defining
    %     the interpretation of the azimuth value: 'radiants' (default) or 'degree'.
    %
    %  Example:
    %     M = mueller_homogeneous_elliptic_diattenuator(1.0,0.5,45,0.5,'deg')
    %     M =
    %
    %        1.0000e+00   1.6542e-17   2.7015e-01   4.2074e-01
    %        1.6542e-17   8.6603e-01   2.3948e-18   3.7297e-18
    %        2.7015e-01   2.3948e-18   9.0514e-01   6.0911e-02
    %        4.2074e-01   3.7297e-18   6.0911e-02   9.6089e-01
    %
    %  References:
    %     [1] E. Collett, Field Guide to Polarization,
    %         SPIE Field Guides vol. FG05, SPIE (2005). ISBN 0-8194-5868-6.
    %     [2] R. A. Chipman, "Polarimetry," chapter 22 in Handbook of Optics II,
    %         2nd Ed, M. Bass, editor in chief (McGraw-Hill, New York, 1995)
    %     [3] "Mueller calculus", http://en.wikipedia.org/wiki/Mueller_calculus,
    %         last retrieved on Dec 17, 2013.
    %     [4] F. Boulvert et al., "Decomposition algorithm of an experimental
    %         Mueller matrix", Opt.Comm. 282(2009):692-704
    %     [5] TODO: check ref[2] in [3] for decomposition in homogeneous elliptic
    %         retarder und diattenuator
    %
    %  See also:
    %     mueller_homogeneous_elliptic_retarder
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    switch(nargin)
        case 0
            t0 = 1;
            d = 0;
            azimuth = 0;
            ellipticity = 0;
            was_cell = false;
            
        case 1
            [t0, was_cell] = c2n(varargin{1}, 1);
            d = 0;
            azimuth = 0;
            ellipticity = 0;
            
        case 2
            [t0, was_cell_t0] = c2n(varargin{1}, 1);
            [d, was_cell_d] = c2n(varargin{2}, 0);
            azimuth = 0;
            ellipticity = 0;
            was_cell = was_cell_t0 || was_cell_d;
            
        case 3
            [t0, was_cell_t0] = c2n(varargin{1}, 1);
            [d, was_cell_d] = c2n(varargin{2}, 0);
            [azimuth, was_cell_azimuth] = c2n(varargin{3}, 0);
            ellipticity = 0;
            was_cell = was_cell_t0 || was_cell_d || was_cell_azimuth;
            
        otherwise
            [t0, was_cell_t0] = c2n(varargin{1}, 1);
            [d, was_cell_d] = c2n(varargin{2}, 0);
            [azimuth, was_cell_azimuth] = c2n(varargin{3}, 0);
            % check special case of homogeneous_elliptic_diattenuator(t0,d,a,'degree')
            if ischar(varargin{4})
                ellipticity = 0;
                was_cell_ellipticity = false;
            else
                [ellipticity, was_cell_ellipticity] = c2n(varargin{4}, 0);
            end
            was_cell = was_cell_t0 || was_cell_d || was_cell_azimuth || was_cell_ellipticity;
            
    end
    
    % convert azimuth angle if necessary
    if nargin>=4 && ischar(varargin{end})
        if strncmpi(varargin{end},'deg',3)
            azimuth = azimuth*pi()/180.0;
        end
    end
    
    % any matrix in parameters?
    if (any([numel(t0),numel(d),numel(azimuth),numel(ellipticity)] > 1)) ...
            || was_cell
        
        % adjust dimensions, i.e. fill missing dimensions with 1
        st0 = size(t0);
        sd = size(d);
        sazimuth = size(azimuth);
        sellipticity = size(ellipticity);
        
        maxdim = max([length(st0),length(sd),length(sazimuth),length(sellipticity)]);
        if length(st0) < maxdim
            st0 = [st0, ones(1,maxdim-length(st0))];
        end
        if length(sd) < maxdim
            sd = [sd, ones(1,maxdim-length(sd))];
        end
        if length(sazimuth) < maxdim
            sazimuth = [sazimuth, ones(1,maxdim-length(sazimuth))];
        end
        if length(sellipticity) < maxdim
            sellipticity = [sellipticity, ones(1,maxdim-length(sellipticity))];
        end
        
        % generate Mueller matrices
        maxsize = max([st0;sd;sazimuth;sellipticity]);
        M = cell(maxsize);
        M_subs = cell(1,ndims(M));
        numelM = numel(M);
        
        % flatten parameter arrays
        t0 = t0(:);
        d = d(:);
        azimuth = azimuth(:);
        ellipticity = ellipticity(:);
        numelt0 = numel(t0);
        numeld = numel(d);
        numelazimuth = numel(azimuth);
        numelellipticity = numel(ellipticity);
        
        for mi=1:numelM
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_homogenenous_elliptic_diattenuator(t0(mod(mi-1,numelt0)+1),...
                d(mod(mi-1,numeld)+1),...
                azimuth(mod(mi-1,numelazimuth)+1), ...
                ellipticity(mod(mi-1,numelellipticity)+1));
        end
        
    else
        
        % generate Mueller matrix
        M = s_homogenenous_elliptic_diattenuator(t0,d,azimuth,ellipticity);
        
    end
    
end

% helper function
function M = s_homogenenous_elliptic_diattenuator(t0,d,azimuth,ellipticity)
    M = zeros(4,4);
    
    a = cos(2*ellipticity)*cos(2*azimuth);
    b = cos(2*ellipticity)*sin(2*azimuth);
    c = sin(2*ellipticity);
    
    M(1,1) = 1;
    M(1,2) = a*d;
    M(1,3) = b*d;
    M(1,4) = c*d;
    M(2,1) = M(1,2);
    M(2,2) = sqrt(1-d^2)+a^2*(1-sqrt(1-d^2));
    M(2,3) = a*b*(1-sqrt(1-d^2));
    M(2,4) = a*c*(1-sqrt(1-d^2));
    M(3,1) = M(1,3);
    M(3,2) = M(2,3);
    M(3,3) = sqrt(1-d^2)+b^2*(1-sqrt(1-d^2));
    M(3,4) = b*c*(1-sqrt(1-d^2));
    M(4,1) = M(1,4);
    M(4,2) = M(2,4);
    M(4,3) = M(3,4);
    M(4,4) = sqrt(1-d^2)+c^2*(1-sqrt(1-d^2));
    
    M = M*t0;
    
end
