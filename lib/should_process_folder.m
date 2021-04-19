function process_folder=should_process_folder(import_dir,reprocess_folder_if_older_than,active_process_mod_time)


process_folder=false;
currently_processing=false;
has_valid_done=false;
% see if there is a done file in this directory
%TODO mask out folders with out/something/done
%check if this data folder has an output dir
out_dir_folder=fullfile(import_dir,'out');
if exist(out_dir_folder,'dir')
    % check that the process time is newer than reprocess_folder_if_older_than
    out_sub_folders=dir(out_dir_folder);
    out_sub_folders=out_sub_folders(3:end);
    out_sub_folders=out_sub_folders(cat(1,out_sub_folders.isdir));
    times_posix=cellfun(@(x) posixtime(datetime(datenum(x,'yyyymmddTHHMMSS'),'TimeZone','local','ConvertFrom','datenum')),{out_sub_folders.name});
    newer_than_reprocess_time=times_posix>reprocess_folder_if_older_than;
    if sum(newer_than_reprocess_time)==0
        process_folder=true;
    else
        %check that that out dir does not have a mod time that is close enough to the current date
        up_to_date_subfolders=out_sub_folders(newer_than_reprocess_time);
        for ii=1:numel(up_to_date_subfolders)
            out_sub_folder_contents=dir(fullfile(out_dir_folder,up_to_date_subfolders(ii).name));
            out_sub_folder_contents=out_sub_folder_contents(3:end);
            file_mod_posix=posixtime(datetime(max([out_sub_folder_contents.datenum]),'TimeZone','local','ConvertFrom','datenum'));
            if abs(file_mod_posix-posixtime(datetime('now','TimeZone','local')))<active_process_mod_time
                currently_processing=true;
            end
            %check each dir thats new enough if it has a done file
            if sum(strcmp({out_sub_folder_contents.name},'Done.txt'))>0 || sum(strcmp({out_sub_folder_contents.name},'done.txt'))>0
                has_valid_done=true;
            end   
        end
        if has_valid_done || currently_processing
            process_folder=false;
        else
            process_folder=true;
        end
    end

else
    process_folder=true;
end

end