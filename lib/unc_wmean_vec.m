function mean_ste_std_vec=unc_wmean_vec(x,unc)
% return the mean and ste as a vector
[mean_val,ste_mean]=unc_wmean(x,unc);
mean_ste_std_vec=cat(1,mean_val,ste_mean,std(x));
end