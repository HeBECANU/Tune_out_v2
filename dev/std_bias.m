%test standard deviation bias



samp_size=col_vec(unique(round(logspace(0,4,1e2))));
resample=1e4;
std_mean=@(m,n) mean(arrayfun(@(x) std(randn(round(m),1)),1:n));
stdc4_mean=@(m,n) mean(arrayfun(@(x) std_c4(randn(round(m),1)),1:n));
bias_with_n=arrayfun(@(x) std_mean(x,resample),samp_size);
ubias_with_n=arrayfun(@(x) stdc4_mean(x,resample),samp_size);

%%
stfig('bias comparison');
clf
loglog(samp_size,bias_with_n)
hold on
loglog(samp_size,ubias_with_n)


%% test with a unifor dist
samp_size=col_vec(unique(round(logspace(0,4,1e2))));
resample=1e4;
std_mean=@(m,n) mean(arrayfun(@(x) std(rand(round(m),1)),1:n));
stdc4_mean=@(m,n) mean(arrayfun(@(x) std_c4(rand(round(m),1)),1:n));
bias_with_n=arrayfun(@(x) std_mean(x,resample),samp_size);
ubias_with_n=arrayfun(@(x) stdc4_mean(x,resample),samp_size);

%%
stfig('bias comparison');
clf
loglog(samp_size,bias_with_n)
hold on
loglog(samp_size,ubias_with_n)
