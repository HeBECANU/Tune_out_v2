function out = nucellf(fn,cell)
    % A shortcut to call a cellfun with non-uniform output
    % Replaces A = cellfun(@(x) foo(x),cell, 'UniformOutput',false)
    %   with   A = nucellf(@(x) foo(x),cell)
    out = cellfun(fn,cell, 'UniformOutput',false);
end