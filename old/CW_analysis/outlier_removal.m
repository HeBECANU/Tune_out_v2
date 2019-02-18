function trim_data = outlier_removal(signal,cutoff)
%INPUTS
%   List of unique set frequencies
%   [f,phasor] pairs
%   # of standard deviations from mean at which to truncate
%OUTPUT
% [setpt,mean,std] list of [float,complex,float]
    set_freqs = unique(signal(:,1));
    num_wls = length(set_freqs);
    resp_means = zeros(num_wls,1);
    resp_std = zeros(num_wls,1);
    keep_means = zeros(num_wls,1);
    keep_std = zeros(num_wls,1);
    for i=1:num_wls
        mask = find(signal(:,1)== set_freqs(i));
        select_responses = signal(mask,2);
        response_mags = abs(select_responses);
        mu = mean(abs(select_responses));
        sigma = std(abs(select_responses));
        resp_means(i) = mu;
        resp_std(i) = sigma;
        data_range = [mu-cutoff*sigma,mu+cutoff*sigma];
        data_mask = (response_mags > data_range(1)) & (response_mags < data_range(2));
        data_keep = select_responses(data_mask);
        keep_means(i) = mean(data_keep);
        keep_std(i) = std(data_keep);
    end
    trim_data = [set_freqs,keep_means,keep_std];
end