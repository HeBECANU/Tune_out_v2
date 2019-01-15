%fix the mislabeling of shots
%because probe and cal shots were swapped
%how it was taken
% shot    type    setpt index	label
% 1       cal       1           probe
% 2       probe     1           cal
% 3       cal       2           probe
% 4       probe     2           cal
% 5       cal       3           probe

% it should be 
% shot    type    setpt index
% 1       probe     1
% 2       cal       1
% 3       probe     2
% 4       cal       2
% 5       probe     3

lv_log_dir_in='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181122_alignment_dep_34_5\log_LabviewMatlab.txt';
lv_log_dir_out='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\20181122_alignment_dep_34_5\log_LabviewMatlab_fixed.txt';
fid = fopen(lv_log_dir_in );
lv_log_cell_in=textscan(fid,'%s','Delimiter','\n');
fclose(fid);
lv_log_cell_in=cellfun(@(x) strsplit(x,','),lv_log_cell_in{1},'UniformOutput',0);

%% 
lv_log_cell_out={};
for ii=1:size( lv_log_cell_in,1)
    line_cells=lv_log_cell_in{ii};
    line_cells_fixed={};
    if isequal(line_cells{4},'measure_probe')
        line_cells_fixed(1:3)=line_cells(1:3);
        line_cells_fixed{4}='calibrate';
        line_cells_fixed{5}=line_cells{6};
        line_cells_fixed{6}=line_cells{7};
    elseif isequal(line_cells{4},'calibrate')
        line_cells_fixed(1:3)=line_cells(1:3);
        line_cells_fixed{4}='measure_probe';
        line_cells_fixed{5}=lv_log_cell_in{ii-1}{5};
        line_cells_fixed{6}=line_cells{5};
        line_cells_fixed{7}=line_cells{6};
    end
    lv_log_cell_out{ii}=line_cells_fixed;
end
%%
lv_log_cell_out=cellfun(@(x) strjoin(x,',') ,lv_log_cell_out,'UniformOutput',0);

fid = fopen(lv_log_dir_out,'w');
for ii=1:size( lv_log_cell_out,2)
    fprintf(fid,'%s\n',lv_log_cell_out{ii});
end
fclose(fid);

