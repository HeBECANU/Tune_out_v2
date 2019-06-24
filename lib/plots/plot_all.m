folder  = 'C:\Users\BEC Machine\Documents\K\Tune_out_v2_trap_freq\lib\plots';
list    = dir(fullfile(folder, '*.m'));
nFile   = length(list);
success = false(1, nFile);
for k = 1:nFile
  file = list(k).name;
  if ~strcmp(file,'plot_all.m')
      try
        run(fullfile(folder, file));
        success(k) = true;
      catch
        fprintf('failed: %s\n', file);
      end
  end
end