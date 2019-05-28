function blue_probe=calc_probe_blue(wm_log_proc,aom_shift)
%output in HZ

blue_probe=[];

%set pt
blue_probe.set=wm_log_proc.probe.freq.set;
blue_probe.set(~isnan(blue_probe.set))=(blue_probe.set(~isnan(blue_probe.set))*1e6*2) +aom_shift*1e6;


blue_probe.act.mean=wm_log_proc.probe.freq.act.mean;
blue_probe.act.mean(~isnan(blue_probe.act.mean))=(blue_probe.act.mean(~isnan(blue_probe.act.mean))*1e6*2) +aom_shift*1e6;


blue_probe.act.std=wm_log_proc.probe.freq.act.std;
blue_probe.act.std(~isnan(blue_probe.act.std))=blue_probe.act.std*1e6*2;

end