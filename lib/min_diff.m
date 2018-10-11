function [min_val,idxy]=min_diff(in)
%fing the minimum absolute difference between every element and all other elements exluding itself
%bugs/improvements
%could get a factor of 2 speedup by only looking for the min in the upper or lower diag;
diff_mat=abs(bsxfun(@minus, in ,in'));
sd=size(diff_mat);
n=sd(1);
diff_mat=diff_mat+diag(NaN(n,1));
[min_val,idx]=nanmin(diff_mat(:));
[idx1,idx2] = ind2sub(sd,idx);
idxy=sort([idx1,idx2] );
end