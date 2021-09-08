folder  = 'C:\Users\jaker\GoogleDrive\HEBEC\Projects\Tune_out_v2\lib\plots';
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