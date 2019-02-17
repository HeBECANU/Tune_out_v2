function out=vec_cellfun(fun,x)
%given some function fun which takes a cell and returns a [1 x A] vector
%apply this function to each element of x [B x 1] and make a matrix out of the result [B x A]
%is about 10x faster than
%a=arrayfun(fun,x,'UniformOutput',0);
%out=cat(1,a{:});

first_out=fun(x(1));
in_vecl_len=size(x,1);
out=repmat(nan*first_out,[in_vecl_len,1]);
out(1,:)=first_out;
for ii=2:in_vecl_len
    out(ii,:)=fun(x(ii));
end

end