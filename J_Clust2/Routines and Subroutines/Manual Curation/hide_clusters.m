function [hide_cl_scat_wires, hide_cl_scat_time] = hide_clusters(unit_pts, edit_cl, val, chan_disp1, feature_wires, chan_disp2, feature_time, hide_cl_scat_wires, hide_cl_scat_time, ts, Channel_Scatter, Time_Scatter)

cl_to_hide = find(edit_cl);

if val == 1
    for i = 1:length(cl_to_hide)
        if cl_to_hide(i) > 13
            continue;
        elseif cl_to_hide(i) > 12
            cur_cl_to_hide = get_selected_pts();
        else
            cur_cl_to_hide = unit_pts{cl_to_hide(i)};
        end
        
        axes(Channel_Scatter)
        hold on
        if length(chan_disp1) > 2
            hide_cl_scat_wires{i} = scatter3(feature_wires(chan_disp1(1), cur_cl_to_hide), feature_wires(chan_disp1(2), cur_cl_to_hide), feature_wires(chan_disp1(3), cur_cl_to_hide), 36, [0 .447 .741], '.');
        else %length(chan_disp) == 2
            hide_cl_scat_wires{i} = scatter(feature_wires(chan_disp1(1), cur_cl_to_hide), feature_wires(chan_disp1(2), cur_cl_to_hide), 36, [0 .447 .741], '.');
        end
        hold off
        
        axes(Time_Scatter)
        hold on
        hide_cl_scat_time{i} = scatter(ts(cur_cl_to_hide), feature_time(chan_disp2, cur_cl_to_hide), 36, [0 .447 .741], '.');
        hold off
    end
else
    for i = 1:length(cl_to_hide)
        delete(hide_cl_scat_wires{i});
        delete(hide_cl_scat_time{i});
    end
end



