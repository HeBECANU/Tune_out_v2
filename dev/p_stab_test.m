fwtext('Running')
profile on
clear all
dir = 'V:\Documents\K\Tune_out_v2_trap_freq\diagnostics';

f{1}.fname = 'DATA04.CSV';
f{1}.label = '320\circ';

f{2}.fname = 'DATA05.CSV';
f{2}.label = '140\circ';

nbins = 25;

num_items = numel(f);
labels = cellfun(@(x) x.label, f, 'UniformOutput',false);

Pdata = cell(num_items,1);
for pt = 1:num_items
    Pdata{pt} = importpowermeter(fullfile(dir,f{pt}.fname));
end

sfigure(66254);
clf;
subplot(2,1,1)
for pt = 1:numel(f)
plot(Pdata{pt}.time,Pdata{pt}.data)
hold on
end
legend(labels)

subplot(2,1,2)
for pt=1:numel(Pdata)
    histogram(Pdata{pt}.data,nbins,'FaceAlpha',0.4)
    hold on
end
legend(labels)
fwtext('Done')
profile off
profile viewer