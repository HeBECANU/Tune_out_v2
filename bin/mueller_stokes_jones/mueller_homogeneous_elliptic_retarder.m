function M = mueller_homogeneous_elliptic_retarder(varargin)
    %function M = mueller_homogeneous_elliptic_retarder(varargin)
    %Return the Mueller matrix for a homogeneous elliptic retarder (see refs 3 and 4).
    %
    %  Syntax:
    %     M = mueller_homogeneous_elliptic_retarder(t0, delay, azimuth, ellipticity)
    %     A = mueller_homogeneous_elliptic_retarder(T0, Delay, Azimuth, Ellipticity)
    %     mueller_homogeneous_elliptic_retarder(..., delaymode)
    %     mueller_homogeneous_elliptic_retarder(..., delaymode, azimuthmode)
    %
    %  Description:
    %     M = mueller_homogeneous_elliptic_retarder(t0, delay, azimuth, ellipticity)
    %     returns the Mueller matrix for a homogeneous elliptic retarder with
    %     total transmission T0 (default: 1), retardation delay (default: 0), and azimuth
    %     (default: 0) and ellipticity (default: 0) describe the two orthogonal polarization eigenstates.
    %
    %     A = mueller_homogeneous_elliptic_retarder(T0, Delay, Azimuth, Ellipticity)
    %     returns a cell array of Mueller matrices for homogeneous elliptic retarders if
    %     any of the four numeric parameters is a numeric array. Then, the size of A is
    %     given by max(size(T0),size(Delay),size(Azimuth),size(Ellipticity))
    %     and elements of smaller matrices of T0, Delay, Azimuth or Ellipticity are
    %     used in a loop-over manner.
    %
    %     mueller_homogeneous_elliptic_retarder(..., delaymode)
    %     mueller_homogeneous_elliptic_retarder(..., delaymode, azimuthmode)
    %     where delaymode and azimuthmode are strings defining the interpretation of the
    %     retardation delay and azimuth values: 'radiants' (default) or 'degree' or 'wavelength' (for delaymode)
    %
    %  Example:
    %     M = mueller_homogeneous_elliptic_retarder(1.0,0.5,45,0.5,'wav','deg')
    %     M =
    %
    %        1.00000   0.00000   0.00000   0.00000
    %        0.00000  -1.00000   0.00000  -0.00000
    %        0.00000   0.00000  -0.41615   0.90930
    %        0.00000  -0.00000   0.90930   0.41615
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
    %         retarder und retarder
    %
    %  See also:
    %     mueller_homogeneous_elliptic_diattenuator
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
            delay = 0;
            azimuth = 0;
            ellipticity = 0;
            was_cell = false;
            
        case 1
            [t0, was_cell] = c2n(varargin{1}, 1);
            delay = 0;
            azimuth = 0;
            ellipticity = 0;
            
        case 2
            [t0, was_cell_t0] = c2n(varargin{1}, 1);
            [delay, was_cell_delay] = c2n(varargin{2}, 0);
            azimuth = 0;
            ellipticity = 0;
            was_cell = was_cell_t0 || was_cell_delay;
            
        case 3
            [t0, was_cell_t0] = c2n(varargin{1}, 1);
            [delay, was_cell_delay] = c2n(varargin{2}, 0);
            % check special case of homogeneous_elliptic_retarder(t0,d,'degree')
            if ischar(varargin{3})
                azimuth = 0;
                was_cell_azimuth = false;
            else
                [azimuth, was_cell_azimuth] = c2n(varargin{3}, 0);
            end
            ellipticity = 0;
            was_cell = was_cell_t0 || was_cell_delay || was_cell_azimuth;
            
        otherwise
            [t0, was_cell_t0] = c2n(varargin{1}, 1);
            [delay, was_cell_delay] = c2n(varargin{2}, 0);
            [azimuth, was_cell_azimuth] = c2n(varargin{3}, 0);
            % check special case of homogeneous_elliptic_retarder(t0,d,a,'degree')
            if ischar(varargin{4})
                ellipticity = 0;
                was_cell_ellipticity = false;
            else
                [ellipticity, was_cell_ellipticity] = c2n(varargin{4}, 0);
            end
            was_cell = was_cell_t0 || was_cell_delay || ...
                was_cell_azimuth || was_cell_ellipticity;
            
    end
    
    % convert delay and azimuth angle if necessary
    delaycvt = [];
    azimuthcvt = [];
    if nargin>=3
        if ischar(varargin{end})
            if ischar(varargin{end-1})
                delaycvt = nargin-1;
                azimuthcvt = nargin;
            else
                delaycvt = nargin;
            end
        end
    end
    if ~isempty(delaycvt)
        if strncmpi(varargin{delaycvt},'deg',3)
            delay = delay*pi()/180.0;
        elseif strncmpi(varargin{delaycvt},'wav',3)
            delay = delay*2*pi();
        end
    end
    if ~isempty(azimuthcvt)
        if strncmpi(varargin{azimuthcvt},'deg',3)
            azimuth = azimuth*pi()/180.0;
        end
    end
    
    % any matrix in parameters?
    if (any([numel(t0),numel(delay),numel(azimuth),numel(ellipticity)] > 1)) || ...
            was_cell
        
        % adjust dimensions, i.e. fill missing dimensions with 1
        st0 = size(t0);
        sdelay = size(delay);
        sazimuth = size(azimuth);
        sellipticity = size(ellipticity);
        
        maxdim = max([length(st0),length(sdelay),length(sazimuth),length(sellipticity)]);
        if length(st0) < maxdim
            st0 = [st0, ones(1,maxdim-length(st0))];
        end
        if length(sdelay) < maxdim
            sdelay = [sdelay, ones(1,maxdim-length(sdelay))];
        end
        if length(sazimuth) < maxdim
            sazimuth = [sazimuth, ones(1,maxdim-length(sazimuth))];
        end
        if length(sellipticity) < maxdim
            sellipticity = [sellipticity, ones(1,maxdim-length(sellipticity))];
        end
        
        % generate Mueller matrices
        maxsize = max([st0;sdelay;sazimuth;sellipticity]);
        M = cell(maxsize);
        M_subs = cell(1,ndims(M));
        
        % flatten parameter arrays
        t0 = t0(:);
        delay = delay(:);
        azimuth = azimuth(:);
        ellipticity = ellipticity(:);
        numelt0 = numel(t0);
        numeldelay = numel(delay);
        numelazimuth = numel(azimuth);
        numelellipticity = numel(ellipticity);
        
        for mi=1:numel(M)
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_homogenenous_elliptic_retarder(t0(mod(mi-1,numelt0)+1),...
                delay(mod(mi-1,numeldelay)+1),...
                azimuth(mod(mi-1,numelazimuth)+1), ...
                ellipticity(mod(mi-1,numelellipticity)+1));
        end
        
    else
        
        % generate Mueller matrix
        M = s_homogenenous_elliptic_retarder(t0,delay,azimuth,ellipticity);
        
    end
    
end

% helper function
function M = s_homogenenous_elliptic_retarder(t0,delay,azimuth,ellipticity)
    
    M = zeros(4,4);
    
    d = cos(2*ellipticity)*cos(2*azimuth)*sin(delay/2);
    e = cos(2*ellipticity)*sin(2*azimuth)*sin(delay/2);
    f = sin(2*ellipticity)*sin(delay/2);
    g = cos(delay/2);
    
    M(1,1) = 1;
    M(2,2) = d^2-e^2-f^2+g^2;
    M(2,3) = 2*(d*e+f*g);
    M(2,4) = 2*(d*f-e*g);
    M(3,2) = M(2,3);
    M(3,3) = -d^2+e^2-f^2+g^2;
    M(3,4) = 2*(e*f-d*g);
    M(4,2) = M(2,4);
    M(4,3) = M(3,4);
    M(4,4) = -d^2-e^2+f^2+g^2;
    
    M = M*t0;
    
end
