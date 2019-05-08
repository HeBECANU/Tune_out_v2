function out=row_col_fun_mat(funhandle,mat_in,dirn)
%apply a function to the rows or column of a matrix
%TODO
%   - speedup improvements
%   - comprehensive test script
%   - make a version to handle arb dimensionality

if dirn==1
    row_slice=mat_in(1,:); %TODO use rotation
    first_out=funhandle(row_slice(:));
    first_out=first_out(:);
    first_out_size=size(first_out);
    %fun_return_row_or_col=first_out_size(2)==1;
    %fun_return_scalar=isscalar(first_out_size);
    iimax=size(mat_in,1); 
    out=zeros(iimax,max(first_out_size)); %TODO use fun_return_row_or_col to chose first_out_size
    out(1,:)=first_out;
    for ii=2:iimax
        row_slice=mat_in(ii,:);
        out_slice=funhandle(row_slice(:));
        out(ii,:)=out_slice(:);
    end
elseif dirn==2
    col_slice=mat_in(:,1);
    first_out=funhandle(col_slice(:));
    first_out=first_out(:);
    first_out_size=size(first_out);
    %fun_return_row_or_col=first_out_size(2)==1;
    %fun_return_scalar=isscalar(first_out_size);
    iimax=size(mat_in,2); 
    out=zeros(max(first_out_size),iimax); %TODO use fun_return_row_or_col to chose first_out_size
    out(:,1)=first_out;
    for ii=2:iimax
        col_slice=mat_in(:,ii);
        out_slice=funhandle(col_slice(:));
        out(:,ii)=out_slice(:);  %TODO use optional rotation here
    end
else 
    error('dirn not 1,2')
end
end