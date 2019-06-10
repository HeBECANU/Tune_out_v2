%corr_cancel_play

%given two vectors u and v find w=u-a*v such that w,v has no correlation

vec_len=1e3;
time=1:vec_len;
vec1=rand(vec_len,1);
vec1=vec1-mean(vec1);
vec1=gaussfilt(time,vec1,10);
vec2=rand(vec_len,1);
vec2=gaussfilt(time,vec2,10);
vec2=vec2-mean(vec2);
%%
leakage_fac=randn(1)
corr_time=50;
vec1_leak=vec1+leakage_fac.*vec2;

[corrvec12,corr_lags]=xcorr(vec1_leak,vec2);
corrvec12_normv2=corrvec12./sum(vec2.*vec2);
subplot(3,1,1)
plot(corr_lags,corrvec12)
subplot(3,1,2)
plot(corr_lags,corrvec12_normv2)
[~,idx]=min(abs(corr_lags))
mean_away_from_peak=mean(corrvec12_normv2(abs(corr_lags)>corr_time));
%[~,idx]=max(abs(corrvec12));
est_xtalk=corrvec12_normv2(idx);
est_xtalk=est_xtalk-mean_away_from_peak

vec1_anti_xtalk=vec1_leak-est_xtalk*vec2;

[corrvec12_antileak,corr_lags]=xcorr(vec1_anti_xtalk,vec2);
subplot(3,1,3)
plot(corr_lags,corrvec12_antileak)