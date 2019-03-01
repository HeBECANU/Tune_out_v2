json_dir = 'C:\Users\BEC Machine\Documents\K\Tune_out_v2_trap_freq\dev\json_files\';
%% qwp
polz_data = [
   150,90.6,245,63.9,152;
   146,118,260,74,344;
   142,77.5,278,53,9;
   138,134,284,58,3;
   134,156,292,59,18;
   130,119,284,31,12;
   154,113,236,58,328;
   162,132,242,35,329;
   177,150,248,8.7,338;
   187,165,253,0.46,343;
   202,178,262,11.5,355;
   220,110,273,44.1,192;
   226,87,290,73,199;
   234,95,226,82,319;
   246,101,50,42,316;
   254,115,240,22.1,143
   ];%
for ii = 1:length(polz_data)
    clear pol_data
    pol_data.power.min = polz_data(ii,4);
    pol_data.power.max = polz_data(ii,2);
    pol_data.phi.min = polz_data(ii,5);
    pol_data.phi.max = polz_data(ii,3);
    pol_data.qwp_ang = polz_data(ii,1);
    pol_data.hwp_ang = 333;
    
    jsonStr = jsonencode(pol_data);
    fid = fopen([json_dir,'pol_data_',num2str(ii),'.json'], 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, jsonStr, 'char');
    fclose(fid);
end
%% hwp
polz_data = [
    0.34,111.0,79.0,20.0,NaN,NaN,304,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190201_to_hwp_111_nuller_reconfig\';
    0.37,130.5,120.0,40.0,nan,nan,314,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190201_to_hwp_131_nuller_reconfig\';
    0.14,151.0,98.0,63.0,nan,nan,325.0,1;
    0.13,170,113,84,nan,nan,335,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_171_nuller_reconfig\';
    0.68,190,133,104,nan,nan,344,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_191_nuller_reconfig\';
    0.94,29,89,120,85,119,354,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_29_nuller_reconfig\';
    0.96,50.5,75,146,nan,nan,6,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190204_to_hwp_51_nuller_reconfig\';
    0.59,70,120,165,nan,nan,15,0;
    0.06,92,77,181,nan,nan,26,1;
    0.52,61,84,154,nan,nan,10,0;
    0.86,46,74,137,nan,nan,3,0;
    0.54,24.5,72,112,nan,nan,352,0;
    0.19,180,93,89,nan,nan,339,0;
    0.12,160,120,70.5,nan,nan,319,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_160_nuller_reconfig\';
    0.33,140.5,112,50,118,50,320,1;
    0.36,121,86,32,nan,nan,310,1;
    0.17,100,72,191,nan,nan,300,1;
    0.07,80,107.5,165,nan,nan,290,0;
    0.25,99.5,72,188,72,187.5,299,1;
    0.37,120,78,221,nan,nan,309,1;
    0.24,140,62,230,nan,nan,320,1;
    0.23,147,87.6,238,nan,nan,324,1;
    0.11,155,60.4,248,nan,nan,327,1;
    0.14,160,76,249,nan,nan,329,1;
    0.20,165,80,255,77,254,332,1,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_165_nuller_reconfig\';
    0.16,171,79,264,nan,nan,335,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_171_nuller_reconfig\';
    0.23,176.5,65,266,nan,nan,338.5,0;
    0.23,176.5,65,266,nan,nan,338.5,0;
    0.33,181,71,271,nan,nan,340,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_181_nuller_reconfig\';
    0.52,187,92,278,nan,nan,342,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_187_nuller_reconfig\';
    0.66,194,78,282,nan,nan,347,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\\20190211_to_hwp_194_nuller_reconfig\';
    0.72,199,82,290,80,291,350,0;
    0.65,208,57,298,nan,nan,353,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_208_nuller_reconfig\';
    0.90,217,78,308,nan,nan,358,0;
    0.70,231,70.1,135,nan,nan,4,0;
    0.64,239.5,53,150,60,151,9.5,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_240_nuller_reconfig\';
    0.58,249,75.1,160,nan,nan,14.5,0,'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190211_to_hwp_250_nuller_reconfig_tenma_setpoint\'];

loop_config.dir = {
    ,
    ,
    ,
    ,
    ,
    ,
    ,
    ,
    ,
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_177_nuller_reconfig_part_b\',
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190210_to_hwp_177_nuller_reconfig_part_a\',
    ,
    ,
    ,
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
    ,
    ,
    ,
    'Z:\EXPERIMENT-DATA\2018_Tune_Out_V2\20190203_to_hwp_151_nuller_reconfig\'
    };

for ii = 1:length(polz_data)
    clear pol_data
    pol_data.power.min = polz_data(ii,1);
    pol_data.power.max = polz_data(ii,3);
    pol_data.phi.min = polz_data(ii,2);
    pol_data.phi.max = polz_data(ii,4);
    pol_data.qwp_ang = nan;
    pol_data.hwp_ang = polz_data(ii,7);
    
    jsonStr = jsonencode(pol_data);
    fid = fopen([json_dir,'pol_data_',num2str(ii+100),'.json'], 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, jsonStr, 'char');
    fclose(fid);
end