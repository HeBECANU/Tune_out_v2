function out=row_col_fun_mat(funhandle,mat_in,dirn)
%apply a function to the rows or column of a matrix
%TODO
%   - speedup improvements
%   - comprehensive test script
%   - make a version to handle arb dimensionality

warning('this wraper will be removed in a future release')
if dirn==2
    dirn=1;
elseif dirn==1
    dirn=2;
end
out=col_row_fun_mat(funhandle,mat_in,dirn);
end