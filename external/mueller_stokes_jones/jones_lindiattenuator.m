function JM = jones_lindiattenuator(varargin)
    %function JM = jones_lindiattenuator(varargin)
    %Return the Jones matrix for a linear diattenuator at zero rotation.
    %
    %  Syntax:
    %     JM =jones_lindiattenuator()
    %     JM =jones_lindiattenuator(d)
    %     A = jones_lindiattenuator(D)
    %     JM =jones_lindiattenuator(px,py)
    %     A = jones_lindiattenuator(PX,py)
    %     A = jones_lindiattenuator(px,PY)
    %     A = jones_lindiattenuator(PX,PY)
    %     jones_lindiattenuator(..., mode)
    %
    %  Description:
    %     JM =jones_lindiattenuator() returns the Jones matrix for a linear
    %     diattenuator with zero (di)attenuation.
    %
    %     JM =jones_lindiattenuator(d) returns the Jones matrix for a linear
    %     diattenuator with px=d, i.e. transmission in y direction
    %     is (1-d)((1+d), whereas transmission in x direction is 1.
    %
    %     A = jones_lindiattenuator(D), where D is a m-by-n-by-... array of
    %     diattentuation values, returns a cell array of
    %     diattenuator Jones matrices with size(A)==size(P).
    %
    %     JM =jones_lindiattenuator(px,py) returns the Jones matrix for a linear
    %     diattenuator with transmittances px and py in x and y direction, respectively.
    %
    %     A = jones_lindiattenuator(PX,py), where PX is a m-by-n-by-... array of
    %     transmittance values in x direction, and py is a scalar value for the
    %     transmittance in y direction, returns a cell array of
    %     diattenuator Jones matrices with size(A)==size(PX).
    %
    %     A = jones_lindiattenuator(px,PY), where PY is a m-by-n-by-... array of
    %     transmittance values in y direction, and px is a scalar value for the
    %     transmittance in x direction, returns a cell array of
    %     diattenuator Jones matrices with size(A)==size(PY).
    %
    %     A = jones_lindiattenuator(PX,PY), where PX and PY are a m-by-n-by-... arrays of
    %     transmittance values, returns a cell array of
    %     diattenuator Jones matrices with size(A)==size(PX)==size(PY).
    %
    %     jones_lindiattenuator(..., mode), where mode is a string defining
    %     the interpretation of px and transmittance values:
    %     'intensity' (default) or 'amplitude'.
    %
    %  Examples:
    %     JM =jones_lindiattenuator(1,0.9)
    %     JM =
    %
    %         1.00000   0.00000
    %         0.00000   0.94868
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
    %     jones_linpolarizer
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
        
        % generate Jones matrices
        maxsize = max([spx;spy]);
        JM = cell(maxsize);
        JM_subs = cell(1,ndims(JM));
        numelJM = numel(JM);
        
        % flatten parameter arrays
        px = px(:);
        py = py(:);
        numelpx = numel(px);
        numelpy = numel(py);
        
        for jmi=1:numelJM
            [JM_subs{:}] = ind2sub(size(JM),jmi);
            JM{JM_subs{:}} = s_function(px(mod(jmi-1,numelpx)+1), py(mod(jmi-1,numelpy)+1));
        end
        
    else
        
        JM = s_function(px, py);
        
    end
    
end

% helper function
function JM = s_lindiattenuator_amp(px_in_amplitude_domain,py_in_amplitude_domain)
    JM = zeros(2,2);
    JM(1,1) = px_in_amplitude_domain;
    JM(2,2) = py_in_amplitude_domain;
end

% helper function 2
function JM = s_lindiattenuator_int(kx_in_intensity_domain,ky_in_intensity_domain)
    JM = zeros(2,2);
    JM(1,1) = sqrt(kx_in_intensity_domain);
    JM(2,2) = sqrt(ky_in_intensity_domain);
end



