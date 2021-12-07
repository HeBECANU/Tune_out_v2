function [array, was_cell] = c2n(array, defaultvalue)
    
    % if possible, convert cell array to numeric array, else return an
    % array of default values of same size
    
    was_cell = false;
    if ~isnumeric(array)
        if iscell(array)
            was_cell = true;
            array2 = cell2mat(array);
            if (ndims(array2)==ndims(array)) && ...
                    all(size(array2)==size(array)) && ...
                    isnumeric(array2)
                array = array2;
            else
                array = ones(size(array))*defaultvalue;
            end
        else
            array(:) = defaultvalue;
        end
    end
    
end

