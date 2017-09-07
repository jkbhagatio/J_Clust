%keyboard_plot_control
%
%Description: This script allows for controlling the plot displays of both axes in the main GUI via the keyboard. The arrow pad adjusts the
%'Channel_Scatter' axis (up and down change the feature displayed, and right and left change the wires displayed). The number pad adjussts the
%'Time_Scatter' axis (8 (up) and 2 (down) change the feature displayed, and 4 (left) and 6 (right) change the wire displayed. 
%

switch eventdata.Key
    case 'rightarrow' %channel_wires plus
        if handles.chan_combo_val == 12
            handles.chan_combo_val = 3;
        else
            handles.chan_combo_val = handles.chan_combo_val + 1;
        end
        handles.chan_disp1 = str2num(handles.chan_combo_contents{handles.chan_combo_val});
        plot_on_channel_scatter;
        
        handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection
        
    case 'leftarrow' %channel_wires minus
        if handles.chan_combo_val == 3
            handles.chan_combo_val = 12;
        else
            handles.chan_combo_val = handles.chan_combo_val - 1;
        end
        handles.chan_disp1 = str2num(handles.chan_combo_contents{handles.chan_combo_val});
        plot_on_channel_scatter;
        
        handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection
        
    case 'downarrow' %feature_wires plus
        if handles.feat_wires_val == 7 %change to '== 9' once wavelets are added
            handles.feat_wires_val = 3;
        else
            handles.feat_wires_val = handles.feat_wires_val + 1;
        end
        handles.feat_disp1 = handles.feat_wires_contents{handles.feat_wires_val};
        plot_on_channel_scatter;
        
        handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection
        
    case 'uparrow' %feature_wires minus
        if handles.feat_wires_val == 3 
            handles.feat_wires_val = 7; %change 7 to 9 once wavelets are added
        else
            handles.feat_wires_val = handles.feat_wires_val - 1; 
        end
        handles.feat_disp1 = handles.feat_wires_contents{handles.feat_wires_val};
        plot_on_channel_scatter;
        
        handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection
        
    case 'numpad6' %time_wires plus
        if handles.chan_time_val == 6
            handles.chan_time_val = 3;
        else
            handles.chan_time_val = handles.chan_time_val + 1;
        end
        handles.chan_disp2 = str2num(handles.chan_time_contents{handles.chan_time_val});
        plot_on_time_scatter;
        
        handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection
    case 'numpad4' %time_wires minus
        if handles.chan_time_val == 3
            handles.chan_time_val = 6;
        else
            handles.chan_time_val = handles.chan_time_val - 1;
        end
        handles.chan_disp2 = str2num(handles.chan_time_contents{handles.chan_time_val});
        plot_on_time_scatter;
        
        handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection
    case 'numpad2'
        if handles.feat_time_val == 7 %change 7 to 9 once wavelets are added
            handles.feat_time_val = 3; 
        else
            handles.feat_time_val = handles.feat_time_val + 1;
        end
        handles.feat_disp2 = handles.feat_time_contents{handles.feat_time_val};
        plot_on_time_scatter;
        
        handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection

    case 'numpad8'
        if handles.feat_time_val == 3
            handles.feat_time_val = 7; %change 7 to 9 once wavelets are added
        else
            handles.feat_time_val = handles.feat_time_val - 1;
        end
        handles.feat_disp2 = handles.feat_time_contents{handles.feat_time_val};
        plot_on_time_scatter;
        
        handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection
end