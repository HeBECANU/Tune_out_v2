function fwtext(args)
    % A function that makes nice headers in the console
    % Works OK with numeric substitution, but needs a length checking
    % function.
    if ischar(args)
        msg = args;
        vals = [];
    elseif iscell(args)
        msg = args{1};
        vals = args{2:end};
    end
    msg_out = sprintf(msg,vals);
    msg_len=length(msg_out);
    blank = '\n============================================================\n';
    width = length(blank);
        
    if msg_len ==0 
    elseif mod(msg_len,2) ==0 && msg_len ~=0
        msg0 = width/2 - msg_len/2;
        msg1 = width/2 + msg_len/2;
        blank(msg0:msg1+1) = [' ',msg_out,' '];
    else
        msg0 = width/2 - (msg_len-1)/2;
        msg1 = width/2 + (msg_len+1)/2;
        blank(msg0:msg1+1) = [' ',msg_out,' '];
    end
    fprintf(blank,vals);
end