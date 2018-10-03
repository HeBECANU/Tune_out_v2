%test_mask_sorted_txy

dataset=rand(round(10^6),3);
[~,I]=sort(dataset(:,1));
dataset=dataset(I,:);
window=rand(3,2);
window(1,:)=(sort(window(1,:))-0.5)*1e-4+0.5;
window(2,:)=sort(window(2,:));
window(3,:)=sort(window(3,:));

%%
iimax=10;
time_search=nan(iimax,1);
time_brute=nan(iimax,1);
for ii=1:iimax

tic;
counts_masked_search=mask_sorted_txy(dataset,window);
time_search(ii)=toc;

tic;
counts_masked_brute=masktxy(dataset,window);
time_brute(ii)=toc;
end


mean(time_search)/mean(time_brute)


if ~isequal(counts_masked_search,counts_masked_brute)
    size(counts_masked_search)
    size(counts_masked_brute)
    error('not equal')
end