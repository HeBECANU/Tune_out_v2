
function sliders_layout(fig, position)
    global src data cmap vidobj format hImage hImageAxes settings;
    global running background setBackgroundData popup_value;
    
    font_size = str2double(settings.font_size);
    margin = str2double(settings.margin);
    
    % Camera properties
    info = imaqhwinfo(settings.adaptor, str2double(settings.device_id));
    range = info.SupportedFormats;
    
    % For cases popup_value is not set (sometimes the file
    % SparrowCam_format.mat is not available)
    if isempty(popup_value)
        popup_value = 1;
    end
    
    
    % Gain
    prop_info_gain = propinfo(src,'Gain');
    gain_info = prop_info_gain.ConstraintValue;

    gain_min = gain_info(1);
    gain_max = gain_info(2);
    
    if gain_min < src.Gain && src.Gain < gain_max
        init_gain_val = src.Gain;
    else
        init_gain_val = gain_min;
    end

    % Exposure
    
    prop_info_exp = propinfo(src,'Exposure');
    exp_info = prop_info_exp.ConstraintValue;

    exp_min = exp_info(1);
    exp_max = exp_info(2);

    if exp_min < src.Exposure && src.Exposure < exp_max
        init_exp_val = src.Exposure;
    else
        init_exp_val = exp_min;
    end
    
    % ------------
    win_pos = get(fig, 'Position');
    
    % Layout variables

    
    label_width = 45;
    label_height = 20;

    panel_cam_size_x = position(3) - 2*margin;
    panel_cam_size_y = position(4) - 2*margin;

    panel_cam_left = position(1) + margin;
    panel_cam_bottom = position(2) + margin;
    
    slider_length = win_pos(3)*panel_cam_size_x - 2*label_width;
    slider_height = label_height;

    % GAIN SLIDER
    
    base_pos = 5;
    panel_cam = uipanel(fig,'Title','SparrowCam','FontSize',font_size,...
                        'Position',...
                        [panel_cam_left panel_cam_bottom panel_cam_size_x panel_cam_size_y]);

    label_gain_val = uicontrol(panel_cam, 'Style','Text', 'String', ['Gain:' num2str(init_gain_val)],...
                                'Position', [label_width base_pos + label_height slider_length label_height],...
                                'FontSize',font_size);

    label_gain_min = uicontrol(panel_cam, 'Style','Text', 'String', gain_min,...
                                'Position', [0 base_pos label_width label_height],...
                                'FontSize',font_size);
    label_gain_max = uicontrol(panel_cam, 'Style','Text', 'String', gain_max,...
                                'Position', [label_width + slider_length base_pos label_width label_height],...
                                'FontSize',font_size);
    
    slider_gain = uicontrol(panel_cam, 'Style','slider',...
                            'Min',gain_min,'Max',gain_max,'Value',init_gain_val,...
                            'Position',[label_width base_pos slider_length slider_height],'SliderStep',[0.005 0.05],...
                            'Callback', @changeGain);
    function changeGain(source,~)
        % here logic of changing camera gain
        src.Gain = source.Value;
        set(label_gain_val, 'String', ['Gain:' num2str(source.Value)])
    end

    % EXPOSURE SLIDER

    base_pos = base_pos + 2*label_height;
    
    label_exp_val = uicontrol(panel_cam, 'Style','Text', 'String', ['Exposure:' num2str(init_exp_val)],...
                                'Position', [label_width base_pos + label_height slider_length label_height],...
                                'FontSize',font_size);
    
    label_exp_min = uicontrol(panel_cam, 'Style','Text', 'String', num2str(exp_min),...
                                'Position', [0 base_pos label_width label_height],...
                                'FontSize',font_size);
    label_exp_max = uicontrol(panel_cam, 'Style','Text', 'String', num2str(exp_max),...
                                'Position', [label_width + slider_length base_pos label_width label_height],...
                                'FontSize',font_size);
    
    slider_exp = uicontrol(panel_cam, 'Style','slider',...
                            'Min',exp_min,'Max',exp_max,'Value', init_exp_val,...
                            'Position',[label_width base_pos slider_length slider_height],'SliderStep',[0.01 0.1],...
                            'Callback', @changeExp);
                        
    function changeExp(source,~)
        % here logic of changing camera exposure
        10^(source.Value)
        src.Exposure = source.Value;
        set(label_exp_val, 'String', ['Exposure[s]:' num2str(10^(source.Value))])
    end
    
    

   

    % BACKGROUND TOGGLE

    base_pos = base_pos + 2*label_height;
    
    button_toggle_background = uicontrol(panel_cam, 'Style','pushbutton', 'String', 'Subtract Background',...
                                'Position', [label_width base_pos slider_length slider_height],...
                                'FontSize',font_size,...
                                'Callback', @toggle_background);
                            
    function toggle_background(source, ~)
        if background
            set(source,'String', 'Subtract Background');
            background = 0;
        else
            set(source,'String', 'Return Background');
            background = 1;
            setBackgroundData = 1;
        end
    end


    base_pos = base_pos + 2*label_height;
    
    % BUTTON TOGGLE PREVIEW
    button_toggle_preview = uicontrol(panel_cam, 'Style','pushbutton', 'String', 'Pause',...
                                'Position', [label_width base_pos slider_length*0.5 slider_height],...
                                'FontSize',font_size,...
                                'Callback', @toggle_preview);
                            
    function toggle_preview(source, ~)
        if running
            stoppreview(vidobj);
            set(source,'String', 'Go');
            running = 0;
        else
            preview(vidobj, hImage);
            set(source,'String', 'Pause');
            running = 1;
        end
    end
    
    % BUTTON SAVE

    button_save = uicontrol(panel_cam, 'Style','pushbutton', 'String', 'Save',...
                                'Position', [label_width+slider_length*0.5 base_pos slider_length*0.5 slider_height],...
                                'FontSize',font_size,...
                                'Callback', @toggle_save);
                            
     function toggle_save(~, ~)
     dd = datestr(now, 'yyyymmdd_HHMMSS');
     filename = [dd,'.tif'];
     imwrite(data, cmap, filename,'tif')
     end
 
   
    
    % POPUP FORMAT


    base_pos = base_pos + 2*label_height;

    popup_format = uicontrol(panel_cam, 'Style','popup', 'String', range,...
                                'Position', [label_width base_pos slider_length slider_height],...
                                'FontSize',font_size,...
                                'Value', popup_value,...
                                'Callback', @change_format);
    function change_format(~, ~)
        val = popup_format.Value;
        maps = popup_format.String;
        format = char(maps(val));
        popup_value = val;
        %save('SparrowCam_format.mat','format', 'popup_value');
        close; 
        SparrowCam;
  
    end

    base_pos = base_pos + 2*label_height;
    
    label_format = uicontrol(panel_cam, 'Style','Text', 'String', ['Current format: ' format],...
                                'Position', [label_width base_pos slider_length slider_height],...
                                'FontSize',font_size);

 
end 


