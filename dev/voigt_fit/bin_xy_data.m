function bined_data=bin_xy_data(opts)

x_dat=opts.x;
y_dat=opts.y;

if ~isfield(opts,'weights')
    opts.weights = ones(size(x_dat));
end
weights=opts.weights;
weights = weights/sum(weights);

x_dat=col_vec(x_dat);
y_dat=col_vec(y_dat);
weights=col_vec(weights);


if ~isfield(opts,'xbin_edges')
    opts.xbin_edges=linspace(min(x_dat),max(x_dat),ceil(numel(x_dat)/10));
else
xbin_edges=opts.xbin_edges;


num_bin=numel(xbin_edges)-1;
bined_data = [];
bined_data.x.mean=nan(num_bin,1);
bined_data.x.std=nan(num_bin,1);
bined_data.x.lims=nan(num_bin,2);
bined_data.y.mean=nan(num_bin,1);
bined_data.y.sd=nan(num_bin,1);
bined_data.y.se=nan(num_bin,1);

for jj=1:num_bin
    if jj==num_bin
        mask=x_dat>=xbin_edges(jj) & x_dat<=xbin_edges(jj+1);
    else
        mask=x_dat>=xbin_edges(jj) & x_dat<xbin_edges(jj+1);
    end
    if sum(mask)>0
        bined_data.x.mean(jj) = nanmean(x_dat(mask));
        bined_data.y.mean(jj)=  wmean(y_dat(mask),weights(mask));
    end
    if sum(mask)>3
        bined_data.x.std(jj) = std(x_dat(mask));
        bined_data.y.sd(jj)=sqrt(nanvar(y_dat(mask),weights(mask)));
        bined_data.y.se(jj)=sewm(y_dat(mask),weights(mask));
    end
    try
        bined_data.x.lims(jj,:) = [min(x_dat(mask)); max(x_dat(mask))];
    catch
        bined_data.x.lims(jj,:) = [0;0];
    end
end

mask=~isnan(bined_data.x.mean);

bined_data.x.mean=bined_data.x.mean(mask);
bined_data.x.std=bined_data.x.std(mask);
bined_data.x.lims=bined_data.x.lims(mask,:);
bined_data.y.mean=bined_data.y.mean(mask);
bined_data.y.sd=bined_data.y.sd(mask);
bined_data.y.se=bined_data.y.se(mask);


end