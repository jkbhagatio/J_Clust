%keyboard_plot_control
%
%Description: This script allows for controlling the plot displays of both axes in the main GUI via the keyboard. The arrow pad adjusts the
%'Channel_Scatter' axis (up and down change the feature displayed, and right and left change the wires displayed). The number pad adjussts the
%'Time_Scatter' axis (8 (up) and 2 (down) change the feature displayed, and 4 (left) and 6 (right) change the wire displayed. 
%

switch eventdata.Key
    case 'd'
        selected_pts = get_selected_pts();
        if ~isempty(selected_pts)
            k = questdlg('Are you sure you want to completely delete these spikes from the dataset?', 'Delete Spikes');
            if ~isempty(k) && strcmp(k,'Yes')
                [~, indxs_d] = intersect(handles.overlaps, selected_pts);
                handles.overlaps(indxs_d) = [];
                if ~handles.preload
                    handles.ts(selected_pts) = [];
                else
                    handles.cur_ts(selected_pts) = [];
                end
                handles.waveforms(:,:,selected_pts) = [];
                
                if ~isempty(handles.unit_pts)
                    for i = 1:length(handles.unit_pts)
                        cur_unit = handles.unit_pts{i};
                        [~, indxs_d] = intersect(cur_unit, selected_pts);
                        cur_unit(indxs_d) = [];
                        handles.unit_pts{i} = cur_unit;
                    end
                    
                    %can't do it within the same loop because we're removing based on indices
                    for i = 1:length(handles.unit_pts)
                        for j = 1:length(selected_pts)
                            cur_unit = handles.unit_pts{i};
                            cur_unit(cur_unit > selected_pts(j)) = cur_unit(cur_unit > selected_pts(j)) - 1;
                            handles.unit_pts{i} = cur_unit;
                        end
                    end
                            
                end
                for i = 1:8
                    if i==6 || i==7 %wavelets and wavelets_c (currently not being calculated)
                        continue;
                    end
                    if i == 4 %PC scores
                        cur_feat = handles.features{i,1}; 
                        cur_feat(:, selected_pts, :) = [];
                        handles.features{i,1} = cur_feat;
                    elseif i == 5 %Concatenated PC scores
                        cur_feat = handles.features{i,1};
                        cur_feat(selected_pts, :) = [];
                        handles.features{i,1} = cur_feat;
                    else %direct waveform features
                        cur_feat = handles.features{i,1};
                        cur_feat(:, selected_pts) = [];
                        handles.features{i,1} = cur_feat;
                    end
                end
                
                plot_on_channel_scatter;
                plot_on_time_scatter;
                
                if handles.preload
                    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts(handles.first_spk:handles.last_spk), handles.overlaps, handles.cluster_info);
                else
                    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts, handles.overlaps, handles.cluster_info);
                end
                
            end   
        end
        


            
    
    
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