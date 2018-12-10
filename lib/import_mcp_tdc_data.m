function [mcp_tdc_data,import_opts]=import_mcp_tdc_data(import_opts)
%import_mcp_tdc_data -a wrapper that allows the import_mcp_tdc_data_core to be cached
cache_opts=[];
cache_opts.verbose=0;
cache_opts.force_cache_load=import_opts.force_load_save;
import_opts=rmfield(import_opts,'force_load_save');
cache_opts.force_recalc=import_opts.force_reimport;
import_opts=rmfield(import_opts,'force_reimport');
outputs=function_cache(cache_opts,@import_mcp_tdc_data_core,{import_opts});
[mcp_tdc_data,import_opts]=outputs{1:2};
end

function  [mcp_tdc_data,import_opts]=import_mcp_tdc_data_core(import_opts)
%import_mcp_tdc_data_core - imports mcp-dld-tdc data into a convineint strucure
%designed to decrease time-to-results this function deals with the tedious importing of data
%the output is a well aranged structure with everything you could want
%
% Syntax:  [data,import_opts]=import_data(import_opts)
%
% Inputs:
%    import_opts.dir            - string,directory of data files
%    import_opts.file_name      -string, file name before number eg d543.txt would be 'd'
%    import_opts.force_forc;    -logical,remake the forc preconverted files
%    import_opts.dld_xy_rot     -float, rotation angle in radians
%    import_opts.txylim         -[3x2]float array,t,x,y limits in seconds,meters
%    import_opts.shot_num       -[1,shots] numbers of data files to import
%
% Outputs:
%    import_opts - struct, only changes made should be to add a \ to the end of import_opts.dir if required
%    mcp_tdc_data - struct
%    mcp_tdc_data.time_create_write-[shots,2]float array, posix time of when the data file was writen
%    mcp_tdc_data.num_counts -[1,shots]float array, number of counts in reconstucted data
%    mcp_tdc_data.counts_txy -{1,shots}cell array containing [counts,3] float array of time,x position, y position in seconds
%                               ,meters
%    mcp_tdc_data.shot_num   -[1,shots] numbers of data files that were imported, note not ness. the same as import_opts.shot_num 
%                               if for example there were missing files
%
% Example: 
%     import_opts.dir='C:\User\data\'
%     import_opts.file_name='d';
%     import_opts.force_reimport=0;
%     import_opts.force_forc=0;
%     import_opts.dld_xy_rot=0.61;
%     xlim=[-20e-3, 20e-3];     %tight XY lims to eliminate hot spot from destroying pulse widths
%     ylim=[-20e-3, 20e-3];
%     tlim=[.4,1.4];
%     import_opts.txylim=[tlim;xlim;ylim];
%     import_opts.shot_num=find_data_files(import_opts);
%     [data,import_opts]=import_data(import_opts);
%     vertcat(data.txy{data.total_num>1e3}) %combine all the data for shots with more than 1e3 counts

% Other m-files required: dld_raw_to_txy,masktxy,data_tcreate,dld_read_5channels_reconst_multi_imp,txy_importer
% Also See:find_data_files
% Subfunctions: none
% MAT-files required: none
%
% Known BUGS/ Possible Improvements
%    -more commenting
%
% Author: Bryce Henson
% email: Bryce.Henson@live.com
% Last revision:2018-08-19

%------------- BEGIN CODE --------------

%do some basic checks on input
%mandatory
if ~isfield(import_opts, 'dir') ,error('Bad Input:dir'), end
if ~isa(import_opts.dir,'char') ,error('Bad Input:dir'), end

%optional
if ~isfield(import_opts, 'mat_save') ,import_opts.mat_save=true; end
if ~isfield(import_opts, 'mod_wait') ,import_opts.mod_wait=1; end
if ~isfield(import_opts, 'force_load_save') ,import_opts.force_load_save=false; end
if ~isfield(import_opts, 'force_reimport') ,import_opts.force_reimport=false; end
if ~isfield(import_opts, 'force_forc') ,import_opts.force_forc=false; end
% these defaults should work for most cases
if ~isfield(import_opts, 'file_name') ,import_opts.file_name='d'; end
if ~isfield(import_opts, 'dld_xy_rot') ,import_opts.dld_xy_rot=0.61; end
if ~isfield(import_opts, 'txylim') ,import_opts.txylim=[[0,10];[-30e-3, 30e-3];[-30e-3, 30e-3]]; end
%if the shot numbers are not specified import everythin in the directory
if ~isfield(import_opts,'shot_num') 
    import_opts.shot_num=find_data_files(import_opts);
end

%check optionals
if ~isa(import_opts.file_name,'char') ,error('Bad Input:file_name'), end
if ~isa(import_opts.force_load_save,'logical') ,error('Bad Input:force_load_save'), end
if ~isa(import_opts.force_reimport,'logical') ,error('Bad Input:force_reimport'), end
if ~isa(import_opts.force_forc,'logical') ,error('Bad Input:force_forc'), end
if ~isa(import_opts.dld_xy_rot,'double') ,error('Bad Input:dld_xy_rot'), end
if ~isa(import_opts.file_name, 'char') ,error('Bad Input:file_name'), end
if size(import_opts.txylim)~=[3,2], error('Bad Input:txylim size'), end
if ~(isa(import_opts.shot_num,'int') || isa(import_opts.shot_num,'double')) ,error('Bad Input:shot_num'), end
if numel(import_opts.shot_num)==0
    error('Bad Input:no shots to be imported')
end

%fix if there is not a trailing file seperator on the directory, use filesep for linux compatablility
if import_opts.dir(end)~=filesep
    import_opts.dir=[import_opts.dir,filesep];
end

mcp_tdc_data=[];

fprintf('importing mcp-tdc files %04i:%04i',size(import_opts.shot_num,2),0)
for ii=1:size(import_opts.shot_num,2)
    convert_dld_to_txy=import_opts.force_forc;
    mcp_tdc_data.shot_num(ii)=nan;
    mcp_tdc_data.num_counts(ii)=nan;
    mcp_tdc_data.counts_txy{ii}={};
    if ~(exist([import_opts.dir,import_opts.file_name,num2str(import_opts.shot_num(ii)),'.txt'],'file')==2)
        fprintf('\n no_file %04i \n %04i\n',ii,ii)
    elseif ~is_dld_done_writing(import_opts.dir,[import_opts.file_name,num2str(import_opts.shot_num(ii)),'.txt'],import_opts.mod_wait)
        fprintf(2,'\n data file not done writing will not process %04i \n %04i\n',import_opts.shot_num(ii),ii)
    else
         mcp_tdc_data.time_create_write(ii,:)=data_tcreate([import_opts.dir,import_opts.file_name],num2str(import_opts.shot_num(ii)));
         %if the txy_forc does not exist, if import_opts.force_forc, or the forc file was earlier than the dld file (re) make it
         if ~convert_dld_to_txy && ...
                 ~(exist([import_opts.dir,import_opts.file_name,'_txy_forc',num2str(import_opts.shot_num(ii)),'.txt'],'file')==2)
             convert_dld_to_txy=true;
         elseif ~convert_dld_to_txy
            %check that the _txy_forc file was created after the raw dld file
            time_forc=data_tcreate([import_opts.dir,import_opts.file_name,'_txy_forc'],num2str(import_opts.shot_num(ii)));
            if time_forc(2) < mcp_tdc_data.time_create_write(ii,2)
            	convert_dld_to_txy=true;
                
            end
         else
         end

         if convert_dld_to_txy
             dld_raw_to_txy([import_opts.dir,import_opts.file_name],import_opts.shot_num(ii),import_opts.shot_num(ii));
         end
         %ineffecient to read back what whas just written
         txydata=txy_importer([import_opts.dir,import_opts.file_name],num2str(import_opts.shot_num(ii)));
         txydata=masktxy(txydata,import_opts.txylim); %mask for counts in the window txylim     
         %rotate the counts into the trap axis
         alpha=-import_opts.dld_xy_rot;
         mcp_tdc_data.counts_txy{ii}=txydata*[1 0 0;0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
         mcp_tdc_data.num_counts(ii)=size(txydata,1);
         mcp_tdc_data.shot_num(ii)=import_opts.shot_num(ii);
    end %file exists condition
    fprintf('\b\b\b\b%04i',ii)
end
import_opts_old=import_opts;
fprintf('\b\b\b\b...Done\n')
    


end


% if isfile('import_mcp_tdc_save.mat') && (import_opts.force_load_save || (~import_opts.force_forc && ~import_opts.force_reimport) )
%     load('import_mcp_tdc_save.mat','import_opts_old')
%     import_opts_old.force_reimport=import_opts.force_reimport; %prevents reimport after import_opts.force_reimport changed 1->0
%     import_opts_old.force_forc=import_opts.force_forc;
%     if isequal(import_opts_old,import_opts) || import_opts.force_load_save
%         import_data=false; 
%         fprintf('import_opts the same loading old data...')
%         load('import_mcp_tdc_save.mat','mcp_tdc_data','import_opts_old')
%         import_opts=import_opts_old;
%         fprintf('Done\n')
%     else
%         clear('data')
%     end
% end


%if ~isa(import_opts.force_reimport,'logical') ,error('Bad Input:force_reimport'), end
%if ~isfield(import_opts, 'force_reimport') ,error('Bad Input:force_reimport'), end
%     %saving the data takes a while, some comparisons:
%     if import_opts.mat_save
%         fprintf('Saving mat file...')
%         save('import_mcp_tdc_save.mat','mcp_tdc_data','import_opts_old','-v7.3');
%         fprintf('Done\n')
%     end
