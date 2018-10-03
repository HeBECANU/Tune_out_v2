function struct=clean_log_structure(struct,path)
if isstruct(struct)
    fnames=fieldnames(struct);
    for fn=fnames'
        struct.(fn{1})=clean_log_structure(struct.(fn{1}),[path,'.',fn{1}]);
    end
else
    if numel(struct)~=1 && iscell(struct) && isnumeric(struct{1})
        empty_mask=cellfun(@(x) isequal(x,[]),struct);
        struct(empty_mask)={nan};
        struct=cell2mat(struct);
    else
        try
        values=unique(struct);
        catch
            warning('debug here')
        end
        if numel(values)==2 %turn 'on' 'off' to logical
            values=unique(struct);
            if sum(strcmp(values,'on')+strcmp(values,'off'))==2
                struct=strcmp(struct,'on'); %convert to logical
            end
        elseif numel(values)==1 %if there is only one unique value set it to that
            struct=values;
        end
        
end
end