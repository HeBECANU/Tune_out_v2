function check=is_alternating_logical_vec(in)
%checks that the input is an alternating logical vector eg [0,1,0,1,0,1]
%example usage:
%   is_alternating_logical_vec(repmat(logical([1,0]'),[1000,1]))
%   is_alternating_logical_vec([0,1,0,1,1])
%   is_alternating_logical_vec([0,0,1,0,1])

in=logical(in);
in=in(:)';
check=islogical(in) && isempty(strfind(in, logical([1,1]))) && isempty(strfind(in, logical([0,0])));

end