function M = mueller_circdiattenuator(varargin)
    %function M = mueller_circdiattenuator(varargin)
    %Return the Mueller matrix for a circular diattenuator at zero rotation.
    %
    %  Syntax:
    %     M = mueller_circdiattenuator()
    %     M = mueller_circdiattenuator(d)
    %     A = mueller_circdiattenuator(D)
    %     M = mueller_circdiattenuator(pr,pl)
    %     A = mueller_circdiattenuator(PR,pl)
    %     A = mueller_circdiattenuator(pr,PL)
    %     A = mueller_circdiattenuator(PR,PL)
    %     mueller_circdiattenuator(..., mode)
    %
    %  Description:
    %     M = mueller_circdiattenuator() returns the Mueller matrix for a circular
    %     diattenuator with zero (di)attenuation.
    %
    %     M = mueller_circdiattenuator(d) returns the Mueller matrix for a circular
    %     diattenuator with pr=d, i.e. transmission in 'right' circular direction
    %     is (1-d)((1+d), whereas transmission in 'left' circular direction is 1.
    %
    %     A = mueller_circdiattenuator(D), where D is a m-by-n-by-... array of
    %     diattentuation values, returns a cell array of
    %     diattenuator Mueller matrices with size(A)==size(P).
    %
    %     M = mueller_circdiattenuator(pr,pl) returns the Mueller matrix for a circular
    %     diattenuator with transmittances pr and pl in 'right' and 'left' circular direction, respectively.
    %
    %     A = mueller_circdiattenuator(PR,pl), where PR is a m-by-n-by-... array of
    %     transmittance values in 'right' circular direction, and pl is a scalar value for the
    %     transmittance in 'left' circular direction, returns a cell array of
    %     diattenuator Mueller matrices with size(A)==size(pr).
    %
    %     A = mueller_circdiattenuator(pr,PL), where PL is a m-by-n-by-... array of
    %     transmittance values in 'left' circular direction, and PR is a scalar value for the
    %     transmittance in 'right' circular direction, returns a cell array of
    %     diattenuator Mueller matrices with size(A)==size(pl).
    %
    %     A = mueller_circdiattenuator(PR,PL), where PR and PL are a m-by-n-by-... arrays of
    %     transmittance values, returns a cell array of
    %     diattenuator Mueller matrices with size(A)==size(pr)==size(pl).
    %
    %     mueller_lindiattentuator(..., mode) where mode is a string defining
    %     the interpretation of pr and transmittance values:
    %     'intensity' (default) or 'amplitude'.
    %
    %  Examples:
    %     M = mueller_circdiattenuator(1,0.5)
    %     M =
    %
    %         0.75000   0.00000   0.00000   0.25000
    %         0.00000   0.70711   0.00000   0.00000
    %         0.00000   0.00000   0.70711   0.00000
    %         0.25000   0.00000   0.00000   0.75000
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
    %     mueller_lindiattenuator mueller_absorber
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    if nargin<1
        
        pr = 1;
        pl = 1;
        was_cell = false;
        
    elseif nargin==1
        
        pr = varargin{1};
        [pr, was_cell] = c2n(pr, 0);
        
        pl = (1-pr)./(1+pr);
        
    else
        
        pr = varargin{1};
        [pr, was_cellr] = c2n(pr, 1);
        
        pl = varargin{2};
        [pl, was_celll] = c2n(pl, 1);
        
        was_cell = was_cellr ||was_celll;
        
    end
    
    % check mode
    s_function = @s_circdiattenuator_int;
    if nargin>=2 && ischar(varargin{end})
        if strncmpi(varargin{end},'amp',3)
            s_function = @s_circdiattenuator_amp;
        end
    end
    
    % any matrix in parameters?
    if (any([numel(pr),numel(pl)] > 1)) || was_cell
        
        % adjust dimensions, i.e. fill missing dimensions with 1
        spr = size(pr);
        spl = size(pl);
        
        maxdim = max([length(spr),length(spl)]);
        if length(spr) < maxdim
            spr = [spr, ones(1,maxdim-length(spr))];
        end
        if length(spl) < maxdim
            spl = [spl, ones(1,maxdim-length(spl))];
        end
        
        % generate Mueller matrices
        maxsize = max([spr;spl]);
        M = cell(maxsize);
        M_subs = cell(1,ndims(M));
        numelM = numel(M);
        
        % flatten parameter arrays
        pr = pr(:);
        pl = pl(:);
        numelpr = numel(pr);
        numelpl = numel(pl);
        
        for mi=1:numelM
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_function(pr(mod(mi-1,numelpr)+1), pl(mod(mi-1,numelpl)+1));
        end
        
    else
        
        M = s_function(pr, pl);
        
    end
    
end

% helper function
function M = s_circdiattenuator_amp(pr_in_amplitude_domain,pl_in_amplitude_domain)
    
    M = zeros(4,4);
    
    M(1,1) = (pr_in_amplitude_domain^2+pl_in_amplitude_domain^2)/2;
    M(1,4) = (pr_in_amplitude_domain^2-pl_in_amplitude_domain^2)/2;
    M(2,2) = pr_in_amplitude_domain*pl_in_amplitude_domain;
    M(3,3) = M(2,2);
    M(4,1) = M(1,4);
    M(4,4) = M(1,1);
    
end

% helper function 2
function M = s_circdiattenuator_int(kr_in_intensity_domain,kl_in_intensity_domain)
    
    M = zeros(4,4);
    
    M(1,1) = (kr_in_intensity_domain+kl_in_intensity_domain)/2;
    M(1,4) = (kr_in_intensity_domain-kl_in_intensity_domain)/2;
    M(2,2) = sqrt(kr_in_intensity_domain*kl_in_intensity_domain);
    M(3,3) = M(2,2);
    M(4,1) = M(1,4);
    M(4,4) = M(1,1);
    
end


