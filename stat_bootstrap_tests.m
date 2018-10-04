%stat_bootstrap_tests
%develop statistical bootstraping to deterine the standard error of complicated functions
%and check it converges to the correct values for known distributions

%there are two ways to bootstrap
% 1. with replacements
%   * estimate the population sd by std(anal(data_sample_with_rep))*sqrt(n_sample)
% 2. without replacements
%  * estimate the population sd by
%  std(anal(data_sample_no_rep))*sqrt(n_sample)/finite_pop_correction
%calculates finite population correction
%https://www.jstor.org/stable/2340569?origin=crossref&seq=2#metadata_info_tab_contents

% Known BUGS/ Possible Improvements
%   -make a function that integerates all this into something that can be easily wrapped arround a
%   dataset
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-10-04

%repeatability
%rng(round(pi*exp(1)*1e3));

%normal data operation=mean
% data= normrnd(0,1,[1e2,1]);
% anal_opp=@(x) mean(x);
% %for simple things tell us what the actual stadnard error in the mean is
% real_ste=1/sqrt(numel(data));

%unit interval uniform
data=rand([1e2,1]);
anal_opp=@(x) mean(x);
real_ste=sqrt((1/12))/sqrt(numel(data));


sample_ste=std(data,1)/sqrt(numel(data)-1);
sample_frac_vec=linspace(1e-2,0.95,10);
repeat_samp_prefactor=1e5;
unc_frac=NaN(numel(sample_frac_vec),2);
unc_pop=NaN(size(sample_frac_vec));
n_total=numel(data);

fprintf('%02u',0)
for ii=1:numel(sample_frac_vec)
    sample_frac=sample_frac_vec(ii);
    n_sample=floor(sample_frac*n_total);
    if n_sample>3
        %the scaling of std(std(x)) is 1/n so lets take more data at smaller sample_frac
        repeat_samp=round(repeat_samp_prefactor*1/sample_frac);
        anal_sample_no_rep=NaN(repeat_samp,1);
        anal_sample_with_rep=anal_sample_no_rep;
        for jj=1:repeat_samp
            anal_sample_with_rep(jj)=anal_opp(randsample(data,n_sample,true));
            anal_sample_no_rep(jj)=anal_opp(randsample(data,n_sample));
        end
        %estimape the pop std for with replacements
        unc_frac(ii,1)=std(anal_sample_with_rep)*sqrt(n_sample);
        finite_pop_corr=sqrt((n_total-n_sample)/(n_total-1));
        unc_frac(ii,2)=std(anal_sample_no_rep)*sqrt(n_sample)/finite_pop_corr;
    end
    fprintf('\b\b%02u',ii);
end
fprintf('..Done\n')

%%



figure(1);
clf
est_anal_unc(:,1)=unc_frac(:,1)/sqrt(numel(data));
mean_est_anal_unc(1)=nanmean(est_anal_unc(:,1));
est_anal_unc(:,2)=unc_frac(:,2)/sqrt(numel(data));
mean_est_anal_unc(2)=nanmean(est_anal_unc(:,2));
plot(sample_frac_vec,est_anal_unc(:,1));
hold on
plot(sample_frac_vec,est_anal_unc(:,2));
xl=xlim(gca);
line(xl,[real_ste,real_ste],'Color','k','LineWidth',2)
line(xl,[sample_ste,sample_ste],'Color','m','LineWidth',2)

hold off
legend('with replacement','without replacement','distribution SE','sample SE')
xlabel('frac data')
ylabel('std subset * sqrt(n)')

%% Matlabs inbuilt

