function [hide_cl_scat_time] = hide_cl_time(unit_pts, edit_cl, val, chan_disp, feature, hide_cl_scat_time, ts, Time_Scatter)

%Description: This .m file is called via the 'Hide Cluster' toggle button in the main GUI, and depending on the toggle state, either 1) hides spikes 
%belonging to a paritcular cluster(s) on the 'Time_Scatter' axis by substituting their cluster-associated color(s) with the 'All Spikes' color, or 
%2) reveals cluster(s) that have been hidden
%
%Input: 'unit_pts' = points belonging to all current units,  'edit_cl' = vector containing info on which cluster(s) to hide, 'val' = toggle state
%of 'Hide Cluster' button, 'chan_disp' = current channel displayed on axis, 'feature' = current feature displayed on axis, 'hide_cl_scat_time' = 
%scatter plot that hides spikes on axis, 'ts' = spike timestamps, 'Time_Scatter' = axis
%
%Output: 'unit_pts' = updated 'unit_pts'
%

cl_to_hide = find(edit_cl);

axes(Time_Scatter)

if val == 1
    hold on
    for i = 1:length(cl_to_hide)
        if cl_to_hide(i) > 13
            continue;
        elseif cl_to_hide(i) > 12
            cur_cl_to_hide = get_selected_pts();
        else
            cur_cl_to_hide = unit_pts{cl_to_hide(i)};
        end
            hide_cl_scat_time{i} = scatter(ts(cur_cl_to_hide), feature(chan_disp, cur_cl_to_hide), 36, [0 .447 .741], '.');
    end
    hold off
else
    for i = 1:length(cl_to_hide)
        delete(hide_cl_scat_time{i});
    end
end

end