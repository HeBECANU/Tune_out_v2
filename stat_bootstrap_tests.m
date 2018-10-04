%statistical bootstraping tests
%developing an approach to find the undertiantly of analysis based on performing the analysis on
%many smaller subsets
% will also investigate the scalling of the uncertianty with the size of this subset

%there are two ways to bootstrap
% 1. with replacements
%   * estimate the population sd by std(anal(data_sample_with_rep))*sqrt(n_sample)
% 2. without replacements
%  * estimate the population sd by
%  std(anal(data_sample_no_rep))*sqrt(n_sample)/finite_pop_correction
%calculates finite population correction
%https://www.jstor.org/stable/2340569?origin=crossref&seq=2#metadata_info_tab_contents

%I am still getting problems with small data population size

rng(round(pi*exp(1)*1e3));
data= normrnd(0,1,[1e2,1]);
anal=@(x) mean(x);
%for simple things tell us what the actual stadnard error in the mean is
real_ste=1/sqrt(numel(data));

sample_frac_vec=linspace(1e-2,0.9,10);
repeat_samp_prefactor=1e4;
unc_frac=NaN(numel(sample_frac_vec),2);
unc_pop=NaN(size(sample_frac_vec));
n_total=numel(data);

fprintf('%02u',0)
for ii=1:numel(sample_frac_vec)
    sample_frac=sample_frac_vec(ii);
    n_sample=floor(sample_frac*n_total);
    %the scaling of std(std(x)) is 1/n so lets take more data at smaller sample_frac
    repeat_samp=round(repeat_samp_prefactor*1/sample_frac);
    anal_sample_no_rep=NaN(repeat_samp,1);
    anal_sample_with_rep=anal_sample_no_rep;
    for jj=1:repeat_samp
        anal_sample_with_rep(jj)=anal(randsample(data,n_sample,true));
        anal_sample_no_rep(jj)=anal(randsample(data,n_sample));
    end
    %estimape the pop std for with replacements
    unc_frac(ii,1)=std(anal_sample_with_rep,1)*sqrt(n_sample-1);
    finite_pop_corr=sqrt((n_total-n_sample)/(n_total-1));
    unc_frac(ii,2)=std(anal_sample_no_rep,1)*sqrt(n_sample-1)/finite_pop_corr;
    fprintf('\b\b%02u',ii);
end
fprintf('..Done\n')

%%
figure(1);
clf
plot(sample_frac_vec,unc_frac(:,1)/sqrt(numel(data)));
hold on
plot(sample_frac_vec,unc_frac(:,2)/sqrt(numel(data)));
xl=xlim(gca);
line(xl,[real_ste,real_ste],'Color','k','LineWidth',2)
hold off
legend('with replacement','without replacement')
xlabel('frac data')
ylabel('std subset * sqrt(n)')