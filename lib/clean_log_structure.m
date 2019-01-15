function struct=clean_log_structure(struct,struct_path)
if isstruct(struct)
    fnames=fieldnames(struct);
    for fn=fnames'
        struct.(fn{1})=clean_log_structure(struct.(fn{1}),[struct_path,'.',fn{1}]);
    end
else
    if numel(struct)~=1 && iscell(struct) && isnumeric(struct{1}) && sum(~cellfun(@(x) isnumeric(x),struct))==0
        empty_mask=cellfun(@(x) isequal(x,[]),struct);
        struct(empty_mask)={nan};
        struct=cell2mat(struct);
    else
        try
            if iscell(struct)
                %make empty cells empty strings
                mask=cellfun(@(x) ~ischar(x),struct);
                struct(mask)={''};
            end
            values=unique(struct);
            if iscell(values)
                %make empty cells empty strigns
                mask=cellfun(@(x) ~isequal(x,''),values);
                values_not_empty=values(mask);
            end
        catch
            warning('debug here')
            values=[];
        end
        %relies on shortcut evaluation
        if iscell(values) && sum(ismember(values_not_empty,{'on','off'}))==numel(values_not_empty) %turn 'on' 'off' to logical
            if sum(strcmp(values,'on')+strcmp(values,'off'))==2
                %empty becomes off here, could use int as logical with empty as nan
                struct=strcmp(struct,'on'); %convert to logical
            end
        elseif numel(values)==1 %if there is only one unique value set it to that
            struct=values;
        end
        
    end
end