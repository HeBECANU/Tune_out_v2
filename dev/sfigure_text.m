function handle_out = sfigure_text(in,add_stack_to_text,show_handle)
% SFIGURE  Create figure window with text (minus annoying focus-theft).
%
% based on sfig by Daniel Eaton, 2005
%
% Bryce Henson 2019-05-15
% See also figure
%
max_text_width=20;

fun_stack_msb_left=1; %should the stack text have the most signifigant function to the left(true)
add_str_left=0; %should the input text be added to the left(true) or right (false) Sensible to be ~fun_stack_msb_order
clip_str_left=1; %clip the end(right)(true) or start(left)(false) of the title text. Sensible to be ~add_str_left
sep_str=':';

if nargin>=1 
    if ischar(in)||isstring(in)
    	if nargin<3 || isempty(show_handle)
            show_handle=false;
        end
        if nargin<2 || isempty(add_stack_to_text)
            add_stack_to_text=false;
        end
        if add_stack_to_text
            fun_stack_str=make_fun_stack_str(fun_stack_msb_left);
            if add_str_left
                fig_str=[in,sep_str,fun_stack_str];
            else
                fig_str=[fun_stack_str,sep_str,in];
            end
        else
            fig_str=in;
        end
        fig_str=clip_str(fig_str,max_text_width,clip_str_left);
        %find all the existing figures
        exisiting_figs=findall(groot, 'Type', 'figure');  
        %was 4.9+1.4% of 0.029s
        %strings=arrayfun(@(x) x.Name,exisiting_figs,'UniformOutput',false);
        %match_idx=find(strcmp(strings,fig_str));
        % slightly faster than making a cell vec of strings
        str_match_vec=arrayfun(@(x) strcmp(x.Name,fig_str),exisiting_figs);
        match_idx=find(str_match_vec);
        if numel(match_idx)>1
            warning('multiple figs with that name exist, will use the first instance')
            match_idx=match_idx(1);
        end    
        if numel(match_idx)==1
            set(0, 'CurrentFigure', exisiting_figs(match_idx).Number)
            exisiting_figs(match_idx).NumberTitle=show_handle;
            handle_out=exisiting_figs(match_idx);
        else
            if show_handle
                handle_out=figure('Name',fig_str);
            else
                handle_out=figure('Name',fig_str,'NumberTitle','off');
            end
            set(handle_out,'color','w')
        end
    elseif isnumeric(in)
        if round(in)~=in
            error('input must be integer')
        end
        if ishandle(in)
            set(0, 'CurrentFigure', in);
            handle_out=gcf;
        else
            handle_out = figure(in);
            set(handle_out,'color','w');
        end
    else
        warning('invalid input as first argument must be string or int')
    end
else
	%if no arguments call the function but set the fig name to be the stack trace
    fun_stack_str=clip_str(make_fun_stack_str(fun_stack_msb_left),max_text_width,clip_str_left);
    sfigure_text(fun_stack_str)
end
end


function str_out=make_fun_stack_str(fun_stack_msb_left)
% fun_stack_msb_left %should the stack text have the most signifigant function to the left(true) or right(false)

sep_str=':';

if nargin==0
    fun_stack_msb_left=1; 
end

fun_stack_obj=dbstack;
if fun_stack_msb_left
    fun_stack_obj=flipud(fun_stack_obj);
    fun_stack_txt=sprintf(['%s',sep_str],fun_stack_obj(1:end-2).name);
else
    fun_stack_txt=sprintf(['%s',sep_str],fun_stack_obj(3:end).name);
end
str_out=fun_stack_txt(1:end-numel(sep_str));   
end


function str_out=clip_str(str_in,max_text_width,clip_left)
% clip_right clip the right side of the text(true) or left(false)

cont_str='...';
if numel(str_in)>max_text_width
    if clip_left
    	str_out=str_in(end-max_text_width-numel(cont_str):end);
    	str_out=['...',str_out];
    else
        str_out=str_in(1:max_text_width-numel(cont_str));
        str_out=[str_out,'...'];
    end
else
    str_out=str_in;
end

end