function A = cellfunc(varargin)
    
    % TODO: add some documentation for this internal function
    
    if nargin == 0
        A = [];
        return;
    end
    
    fh = varargin{1};
    if ~isa(fh, 'function_handle')
        A = fh;
        return;
    end
    
    % are there additional parameters?
    if nargin == 1
        A = fh();
        return;
    end
    
    % collect all parameters in cell array M
    fhparameter_offset = 1;
    margin = nargin-fhparameter_offset;
    returncellarray = false;
    Msizes = zeros(1,0);
    M = cell(1,margin);
    numelM = zeros(1,margin);
    for mi=1:margin
        mai = varargin{mi+fhparameter_offset};
        if iscell(mai)
            returncellarray = true;
        else
            % convert numeric array to cell(1,1) arrays to ease handling below
            mai = {mai};
        end
        % adjust dimensions, i.e. fill missing dimensions with 1
        sizemai = size(mai);
        lengthmai = length(sizemai);
        if lengthmai > size(Msizes,2)
            Msizes(:,end+1:lengthmai) = 1;
        end
        Msizes(end+1,:) = sizemai; %#ok<AGROW>
        % flatten arrays
        M{mi} = mai(:);
        numelM(mi) = numel(M{mi});
    end
    
    % generate output array A
    Asize = max(Msizes);
    A = cell(Asize);
    A_subs = cell(1,ndims(A));
    numelA = numel(A);
    for ai=1:numelA
        
        % collect specific parameters for current output cell
        b = cell(1, margin);
        for mi=1:margin
            m = M{mi};
            % we are using parameters in a looped manner
            b{mi} = m{mod(ai-1,numelM(mi))+1};
        end
        
        [A_subs{:}] = ind2sub(Asize,ai);
        A{A_subs{:}} = fh(b{:});
    end
    
    if ~returncellarray
        A = A{1};
    end
    
end
