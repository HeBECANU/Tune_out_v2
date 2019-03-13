json_dir = 'C:\Users\BEC Machine\Documents\K\Tune_out_v2_trap_freq\dev\json_files\';
json_dir = 'V:\Documents\K\Tune_out_v2_trap_freq\dev\json_files\';
%% HWP
hwp_dir = {
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_51_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_29_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190201_to_hwp_131_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_171_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_191_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_208_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\\20190211_to_hwp_194_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_187_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_181_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_177_nuller_reconfig_part_b\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_177_nuller_reconfig_part_a\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_171_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_165_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_160_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_155_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_145_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_140_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190208_to_hwp_120_nuller_reconfig_okish\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190208_to_hwp_99_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190207_to_hwp_80_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_121_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_24_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_46_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_61_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_92_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_230_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_217_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_199_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190206_to_hwp_100_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_141_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_160_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_180_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_70_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_240_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_250_nuller_reconfig_tenma_setpoint\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190201_to_hwp_111_nuller_reconfig\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_151_nuller_reconfig\'
};
%% QWP
qwp_dir = {
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190304_qwp_280_pure',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190302_qwp_283_pure_long',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190301_qwp_283_pure_run',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190301_qwp_246_2',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_270',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_286',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_310',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_286_no_analog',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_254',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_246',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_234',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_226',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_220',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_202',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_187',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_177',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_162',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_154',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_130',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_134',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_138',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_142',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_146',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190223_qwp_150',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_260'
};
%%
pol_data_post = zeros(numel(qwp_dir)+numel(hwp_dir),7);
%% hwp
polz_data = {
    0.34,111.0,79.0,20.0,NaN,NaN,304,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190201_to_hwp_111_nuller_reconfig\';
    0.37,130.5,120.0,40.0,nan,nan,314,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190201_to_hwp_131_nuller_reconfig\';
    0.14,151.0,98.0,63.0,nan,nan,325.0,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_151_nuller_reconfig\';
    0.13,170,113,84,nan,nan,335,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_171_nuller_reconfig\';
    0.68,190,133,104,nan,nan,344,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_191_nuller_reconfig\';
    0.94,29,89,120,85,119,354,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_29_nuller_reconfig\';
    0.96,50.5,75,146,nan,nan,6,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_51_nuller_reconfig\';
    0.59,70,120,165,nan,nan,15,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_70_nuller_reconfig\';
    0.06,92,77,181,nan,nan,26,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_92_nuller_reconfig\';
    0.52,61,84,154,nan,nan,10,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_61_nuller_reconfig\';
    0.86,46,74,137,nan,nan,3,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_46_nuller_reconfig\';
    0.54,24.5,72,112,nan,nan,352,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_24_nuller_reconfig\';
    0.19,180,93,89,nan,nan,339,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_180_nuller_reconfig\';
    0.12,160,120,70.5,nan,nan,319,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_160_nuller_reconfig\';
    0.33,140.5,112,50,118,50,320,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_140_nuller_reconfig\';
    0.36,121,86,32,nan,nan,310,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_121_nuller_reconfig\';
    0.17,100,72,191,nan,nan,300,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190206_to_hwp_100_nuller_reconfig\';
    0.07,80,107.5,165,nan,nan,290,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190207_to_hwp_80_nuller_reconfig\';
    0.25,99.5,72,188,72,187.5,299,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190208_to_hwp_99_nuller_reconfig\';
    0.37,120,78,221,nan,nan,309,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190208_to_hwp_120_nuller_reconfig_okish\';
    0.24,140,62,230,nan,nan,320,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_141_nuller_reconfig\';
    0.23,147,87.6,238,nan,nan,324,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_145_nuller_reconfig\';
    0.11,155,60.4,248,nan,nan,327,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190209_to_hwp_155_nuller_reconfig\';
    0.14,160,76,249,nan,nan,329,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_160_nuller_reconfig\';
    0.20,165,80,255,77,254,332,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_165_nuller_reconfig\';
    0.16,171,79,264,nan,nan,335,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_171_nuller_reconfig\';
    0.23,176.5,65,266,nan,nan,338.5,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_177_nuller_reconfig_part_b\';
    0.23,176.5,65,266,nan,nan,338.5,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_177_nuller_reconfig_part_a\';
    0.33,181,71,271,nan,nan,340,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_181_nuller_reconfig\';
    0.52,187,92,278,nan,nan,342,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_187_nuller_reconfig\';
    0.66,194,78,282,nan,nan,347,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\\20190211_to_hwp_194_nuller_reconfig\';
    0.72,199,82,290,80,291,350,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_199_nuller_reconfig\';
    0.65,208,57,298,nan,nan,353,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_208_nuller_reconfig\';
    0.90,217,78,308,nan,nan,358,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_217_nuller_reconfig\';
    0.70,231,70.1,135,nan,nan,4,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_230_nuller_reconfig\';
    0.64,239.5,53,150,60,151,9.5,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_240_nuller_reconfig\';
    0.58,249,75.1,160,nan,nan,14.5,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_250_nuller_reconfig_tenma_setpoint\'};

for ii = 1:length(polz_data)
    clear pol_data
    pol_data.power.min = polz_data{ii,1};
    pol_data.power.max = polz_data{ii,3};
    pol_data.phi.min = polz_data{ii,2};
    pol_data.phi.max = polz_data{ii,4};
    pol_data.qwp_ang = nan;
    pol_data.hwp_ang = polz_data{ii,7};
    pol_data.handedness = 2*polz_data{ii,8}-1;
    indx = find(strcmp(hwp_dir,polz_data{ii,9}));
    pol_data_post(indx,:) = [pol_data.qwp_ang,pol_data.hwp_ang,pol_data.power.max,pol_data.phi.max,pol_data.power.min,pol_data.phi.min,pol_data.handedness];
    jsonStr = jsonencode(pol_data);
    fid = fopen([json_dir,'pol_data_',num2str(ii+100),'.json'], 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, jsonStr, 'char');
    fclose(fid);
end
%% qwp
polz_data = {
   150,90.6,245,63.9,152,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190223_qwp_150';
   146,118,260,74,344,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_146';
   142,77.5,278,53,9,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_142';
   138,134,284,58,3,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_138';
   134,156,292,59,18,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_134';
   130,119,284,31,12,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190224_qwp_130';
   154,113,236,58,328,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_154';
   162,132,242,35,329,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_162';
   177,150,248,8.7,338,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_177';
   187,165,253,0.46,343,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_187';
   202,178,262,11.5,355,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_202';
   220,110,273,44.1,192,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_220';
   226,87,290,73,199,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190225_qwp_226';
   234,95,226,82,319,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_234';
   246,101,50,42,316,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_246';
   254,115,240,22.1,143,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_254';
   260,124,51,18.1,329,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190226_qwp_260';
   286,1.2,201,0.05,352,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_286_no_analog'; %again don't use in final value
   310,140.7,117,40.8,205,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_310';
   286.5,165,88.5,2.18,179,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_286';
   270,167.3,248,5.37,159,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190227_qwp_270';
   246,127,29,45,136,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190301_qwp_246_2';
   283,195,83,0.6,352,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190301_qwp_283_pure_run';
   283,195,83,0.6,352,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190302_qwp_283_pure_long';
   280,168,257,0.13,171,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190304_qwp_280_pure'
   };%
temp = cell2mat(polz_data(:,1));
sign_flip = ones(size(polz_data,1),1).*(1-2.*(45>abs(temp-235)));
for ii = 1:length(polz_data)
    clear pol_data
    pol_data.power.min = polz_data{ii,4};
    pol_data.power.max = polz_data{ii,2};
    pol_data.phi.min = polz_data{ii,5};
    pol_data.phi.max = polz_data{ii,3};
    pol_data.qwp_ang = polz_data{ii,1};
    pol_data.handedness = sign_flip(ii);
    pol_data.hwp_ang = 333;
    indx = find(strcmp(qwp_dir,polz_data{ii,6}));
    pol_data_post(indx+numel(hwp_dir),:) = [pol_data.qwp_ang,pol_data.hwp_ang,pol_data.power.max,pol_data.phi.max,pol_data.power.min,pol_data.phi.min,pol_data.handedness];
    jsonStr = jsonencode(pol_data);
    fid = fopen([json_dir,'pol_data_',num2str(ii),'.json'], 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, jsonStr, 'char');
    fclose(fid);
end