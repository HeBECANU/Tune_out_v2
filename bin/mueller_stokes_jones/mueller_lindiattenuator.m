function M = mueller_lindiattenuator(varargin)
    %function M = mueller_lindiattenuator(varargin)
    %Return the Mueller matrix for a linear diattenuator at zero rotation.
    %
    %  Syntax:
    %     M = mueller_lindiattenuator()
    %     M = mueller_lindiattenuator(d)
    %     A = mueller_lindiattenuator(D)
    %     M = mueller_lindiattenuator(px,py)
    %     A = mueller_lindiattenuator(PX,py)
    %     A = mueller_lindiattenuator(px,PY)
    %     A = mueller_lindiattenuator(PX,PY)
    %     mueller_lindiattenuator(..., mode)
    %
    %  Description:
    %     M = mueller_lindiattenuator() returns the Mueller matrix for a linear
    %     diattenuator with zero (di)attenuation.
    %
    %     M = mueller_lindiattenuator(d) returns the Mueller matrix for a linear
    %     diattenuator with px=d, i.e. transmission in y direction
    %     is (1-d)((1+d), whereas transmission in x direction is 1.
    %
    %     A = mueller_lindiattenuator(D), where D is a m-by-n-by-... array of
    %     diattentuation values, returns a cell array of
    %     diattenuator Mueller matrices with size(A)==size(P).
    %
    %     M = mueller_lindiattenuator(px,py) returns the Mueller matrix for a linear
    %     diattenuator with transmittances px and py in x and y direction, respectively.
    %
    %     A = mueller_lindiattenuator(PX,py), where PX is a m-by-n-by-... array of
    %     transmittance values in x direction, and py is a scalar value for the
    %     transmittance in y direction, returns a cell array of
    %     diattenuator Mueller matrices with size(A)==size(PX).
    %
    %     A = mueller_lindiattenuator(px,PY), where PY is a m-by-n-by-... array of
    %     transmittance values in y direction, and px is a scalar value for the
    %     transmittance in x direction, returns a cell array of
    %     diattenuator Mueller matrices with size(A)==size(PY).
    %
    %     A = mueller_lindiattenuator(PX,PY), where PX and PY are a m-by-n-by-... arrays of
    %     transmittance values, returns a cell array of
    %     diattenuator Mueller matrices with size(A)==size(PX)==size(PY).
    %
    %     mueller_lindiattenuator(..., mode), where mode is a string defining
    %     the interpretation of px and transmittance values:
    %     'intensity' (default) or 'amplitude'.
    %
    %  Examples:
    %     M = mueller_lindiattenuator(1,0.9)
    %     M =
    %
    %         0.9500    0.0500         0         0
    %         0.0500    0.9500         0         0
    %              0         0    0.9487         0
    %              0         0         0    0.9487
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
    %     mueller_circdiattenuator mueller_absorber
    %
    %  File information:
    %     version 1.0 (jan 2014)
    %     (c) Martin Vogel
    %     email: matlab@martin-vogel.info
    %
    %  Revision history:
    %     1.0 (jan 2014) initial release version
    
    if nargin<1
        
        px = 1;
        py = 1;
        was_cell = false;
        
    elseif nargin==1
        
        px = varargin{1};
        [px, was_cell] = c2n(px, 0);
        py = (1-px)./(1+px);
        px(:) = 1;
        
    else
        
        px = varargin{1};
        [px, was_cellx] = c2n(px, 1);
        
        py = varargin{2};
        [py, was_celly] = c2n(py, 1);
        
        was_cell = was_cellx || was_celly;
        
    end
    
    % check mode
    s_function = @s_lindiattenuator_int;
    if nargin>=2 && ischar(varargin{end})
        if strncmpi(varargin{end},'amp',3)
            s_function = @s_lindiattenuator_amp;
        end
    end
    
    % any matrix in parameters?
    if (any([numel(px),numel(py)] > 1)) || was_cell
        
        % adjust dimensions, i.e. fill missing dimensions with 1
        spx = size(px);
        spy = size(py);
        
        maxdim = max([length(spx),length(spy)]);
        if length(spx) < maxdim
            spx = [spx, ones(1,maxdim-length(spx))];
        end
        if length(spy) < maxdim
            spy = [spy, ones(1,maxdim-length(spy))];
        end
        
        % generate Mueller matrices
        maxsize = max([spx;spy]);
        M = cell(maxsize);
        M_subs = cell(1,ndims(M));
        numelM = numel(M);
        
        % flatten parameter arrays
        px = px(:);
        py = py(:);
        numelpx = numel(px);
        numelpy = numel(py);
        
        for mi=1:numelM
            [M_subs{:}] = ind2sub(size(M),mi);
            M{M_subs{:}} = s_function(px(mod(mi-1,numelpx)+1), py(mod(mi-1,numelpy)+1));
        end
        
    else
        
        M = s_function(px, py);
        
    end
    
end

% helper function
function M = s_lindiattenuator_amp(px_in_amplitude_domain,py_in_amplitude_domain)
    
    M = zeros(4,4);
    
    M(1,1) = (px_in_amplitude_domain^2+py_in_amplitude_domain^2)/2;
    M(1,2) = (px_in_amplitude_domain^2-py_in_amplitude_domain^2)/2;
    M(2,1) = M(1,2);
    M(2,2) = M(1,1);
    M(3,3) = px_in_amplitude_domain*py_in_amplitude_domain;
    M(4,4) = M(3,3);
    
end

% helper function 2
function M = s_lindiattenuator_int(kx_in_intensity_domain,ky_in_intensity_domain)
    
    M = zeros(4,4);
    
    M(1,1) = (kx_in_intensity_domain+ky_in_intensity_domain)/2;
    M(1,2) = (kx_in_intensity_domain-ky_in_intensity_domain)/2;
    M(2,1) = M(1,2);
    M(2,2) = M(1,1);
    M(3,3) = sqrt(kx_in_intensity_domain*ky_in_intensity_domain);
    M(4,4) = M(3,3);
    
end


