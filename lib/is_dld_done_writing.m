function[logic_file_done_writing]= is_dld_done_writing(dir,new_file,wait_for_mod)
%this function checks that the modification date of a file is
%wait_for_mod seconds in the past and then checks that the last charater of the file is a newline

file_pointer=fullfile(dir,new_file);
%get the last mod time
if isequal(computer,'PCWIN64')   
    %only works for win64
    reply=GetFileTime(file_pointer);
    time_posix_write=posixtime(datetime(reply.Write));
else
    file_dir=dir(file_pointer);
    time_posix_write=posixtime(datetime(file_dir.datenum,'ConvertFrom','datenum'));
end
time_posix_now=posixtime(datetime('now'));
logic_file_done_writing=(time_posix_write+wait_for_mod)<time_posix_now;

if logic_file_done_writing
    %checks that the last charater of the file is a newline
    %modified from https://stackoverflow.com/questions/2659375/matlab-command-to-access-the-last-line-of-each-file
    fid = fopen(file_pointer,'r');    %# Open the file as a binary
    offset = 30;                      %# Offset from the end of file (bytes)
    fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
    new_char = fread(fid,1,'*char');  %# Read one character
    last_char=''; %initalize to deal with empty file case
    while  numel(new_char)>0
        last_char = new_char;   %# Add the character to a string
        new_char = fread(fid,1,'*char');  %# Read one character
    end
    fclose(fid);  %# Close the file
    logic_file_done_writing=isequal(last_char,newline);
end

end