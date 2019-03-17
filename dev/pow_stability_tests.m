fwtext('Running')
clear all
dir = 'V:\Documents\K\Tune_out_v2_trap_freq\diagnostics';

f{1}.fname = 'DATA04.CSV';
f{1}.label = '320\circ';

f{2}.fname = 'DATA05.CSV';
f{2}.label = '140\circ';

num_items = numel(f)
labels = cellfun(@(x) x.label, f, 'UniformOutput',false);

nbins = 25;

Pdata{1} = importpowermeter(fullfile(dir,f4));
% labels{1}='320\circ';

Pdata{2} = importpowermeter(fullfile(dir,f5));
% labels{2}='140\circ';


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