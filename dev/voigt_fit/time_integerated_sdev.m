function [sdev,unc]=time_integerated_sdev(t,x,tau)
% this function finds how the standard deviation of x(t>t0 & t<t0+tau) averaged over all t
% not sure what the formal name is
% Bryce Henson 2020-04-24
[~,sort_idx]=sort(t);
t=t(sort_idx);
x=x(sort_idx);

iimax=numel(tau);
sdev=nan(iimax,1);
unc=nan(iimax,1);
trange=range(t);
for ii=1:iimax 
    this_tau=tau(ii);
    jjmax=floor(trange/this_tau);
    if jjmax>3
        tmp_sdevs=nan(jjmax,1);
        for jj=1:jjmax
            tmin=this_tau*(jj-1); 
            tmax=this_tau*jj;
            mask_idx= fast_sorted_mask(t,tmin,tmax);
            if diff(mask_idx)>2
                tmp_sdevs(jj)=std(x(mask_idx(1):mask_idx(2)));
            end
        end
        sdev(ii)=nanmean(tmp_sdevs);
        unc(ii)=nanstd(tmp_sdevs)/sqrt(sum(~isnan(tmp_sdevs)));
    end
end


end

