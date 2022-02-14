function labels_layout(fig, position)
    global settings
    global edit_BRan_x edit_BRan_y edit_FWHM_x edit_FWHM_y edit_maxpos_x edit_maxpos_y edit_maxval edit_fps;
    % ------------
    % To get the width and height for position calculations
    win_pos = get(fig, 'Position');
    
    % Layout variables
    font_size = str2double(settings.font_size);
    margin = str2double(settings.margin);
    
    padding = win_pos(3)*0.01;

    panel_cam_size_x = position(3) - 2*margin;
    panel_cam_size_y = position(4) - 2*margin;

    panel_cam_left = position(1) + margin;
    panel_cam_bottom = position(2) + margin;
    

    label_width = win_pos(3)*panel_cam_size_x/2;
    label_height = 20;
    
    
    panel_beam = uipanel(fig,'Title','Beam','FontSize',font_size,...
                        'Position',...
                        [panel_cam_left panel_cam_bottom panel_cam_size_x panel_cam_size_y]);

    % MOST BOTTOM Field
    
    base_pos = padding;
    
    label_BRan_y = uicontrol(panel_beam, 'Style', 'Text', 'String', '2Wu@13.5% (Y):',...
                                'Position', [padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);

    edit_BRan_y = uicontrol(panel_beam, 'Style','Edit', 'String', 0,...
                                'Position', [label_width + padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);
                            
    
    
    % NEXT Field UPWARDS
    
    base_pos = base_pos + 1.5*label_height;
    
    label_BRan_x = uicontrol(panel_beam, 'Style', 'Text', 'String', '2Wu@13.5% (X):',...
                                'Position', [padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);

    edit_BRan_x = uicontrol(panel_beam, 'Style','Edit', 'String', 0,...
                                'Position', [label_width + padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);
    
    % NEXT Field UPWARDS
    
    base_pos = base_pos + 1.5*label_height;
    
    label_FWHM_y = uicontrol(panel_beam, 'Style', 'Text', 'String', '2Wu@50.0% (Y):',...
                                'Position', [padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);

    edit_FWHM_y = uicontrol(panel_beam, 'Style','Edit', 'String', 0,...
                                'Position', [label_width + padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);
                            
    % NEXT Field UPWARDS

    base_pos = base_pos + 1.5*label_height;

    label_FWHM_x = uicontrol(panel_beam, 'Style', 'Text', 'String', '2Wu@50.0% (X):',...
                                'Position', [padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);

    edit_FWHM_x = uicontrol(panel_beam, 'Style','Edit', 'String', 0,...
                                'Position', [label_width + padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);
                            
    % NEXT Field UPWARDS

    base_pos = base_pos + 1.5*label_height;

    label_maxpos_y = uicontrol(panel_beam, 'Style', 'Text', 'String', 'Max Position (Y):',...
                                'Position', [padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);

    edit_maxpos_y = uicontrol(panel_beam, 'Style','Edit', 'String', 0,...
                                'Position', [label_width + padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);                            
    % NEXT Field UPWARDS

    base_pos = base_pos + 1.5*label_height;

    label_maxpos_x = uicontrol(panel_beam, 'Style', 'Text', 'String', 'Max Position (X):',...
                                'Position', [padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);

    edit_maxpos_x = uicontrol(panel_beam, 'Style','Edit', 'String', 0,...
                                'Position', [label_width + padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);
                            
    % NEXT Field UPWARDS

    base_pos = base_pos + 1.5*label_height;

    label_maxval = uicontrol(panel_beam, 'Style', 'Text', 'String', 'Max Intensity:',...
                                'Position', [padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);

    edit_maxval = uicontrol(panel_beam, 'Style','Edit', 'String', 0,...
                                'Position', [label_width + padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size);
                                      
    % FPS Label

    base_pos = base_pos + 2.5*label_height;

    label_fps = uicontrol(panel_beam, 'Style', 'Text', 'String', 'FPS:',...
                                'Position', [padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size-2);

    edit_fps = uicontrol(panel_beam, 'Style','Edit', 'String', 0,...
                                'Position', [label_width + padding base_pos label_width - 2*padding label_height],...
                                'FontSize',font_size-2);
                            
end
% slider_exposure = 