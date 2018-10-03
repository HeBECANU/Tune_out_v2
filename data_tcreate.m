function time_posix=data_tcreate(data_path,filenum_str)
%just return the write time of a data file
%i think the creation time would be closer to what i want however there is no way to have matlab get
%this with inbuilts so ill compromise and use the write time
data_file_path=[data_path,filenum_str,'.txt'];
if isequal(computer,'PCWIN64')
    reply=GetFileTime(data_file_path);
    time_posix=[posixtime(datetime(reply.Creation)),posixtime(datetime(reply.Write))];
else
    file_dir=dir(data_file_path);
    time_posix=[nan,posixtime(datetime(file_dir.datenum,'ConvertFrom','datenum'))];
end
end