function [hide_spks_scat_time] = hide_spks_time(unit_pts, edit_cl, val, chan_disp, feature, hide_spks_scat_time, ts, Time_Scatter)

%Description: This .m file is called via the 'Hide Spikes' toggle button in the main GUI, and depending on the toggle state, either 1) hides spikes 
%belonging to a paritcular cluster(s) on the 'Time_Scatter' axis by removing their cluster-associated color(s), or 2) reveals cluster(s) that have 
%been hidden
%
%Input: 'unit_pts' = points belonging to all current units,  'edit_cl' = vector containing info on which cluster(s) to hide, 'val' = toggle state
%of 'Hide Cluster' button, 'chan_disp' = current channel displayed on axis, 'feature' = current feature displayed on axis, 'hide_spks_scat_time' = 
%scatter plot that hides spikes on axis, 'ts' = spike timestamps, 'Time_Scatter' = axis
%
%Output: 'unit_pts' = updated 'unit_pts'
%

cl_to_hide = find(edit_cl);

axes(Time_Scatter)

if val
    hold on
    for i = 1:length(cl_to_hide)
        if cl_to_hide(i) > 13
            all_pts = [1:length(feature)];
            unit_pts_all = [];
            for j = 1:length(unit_pts)
                unit_pts_all = [unit_pts_all, unit_pts{j}];
            end
            cur_cl_to_hide = all_pts;
            cur_cl_to_hide(unit_pts_all) = [];
        elseif cl_to_hide(i) > 12
            cur_cl_to_hide = get_selected_pts();
        else
            cur_cl_to_hide = unit_pts{cl_to_hide(i)};
        end
            hide_spks_scat_time{i} = scatter(ts(cur_cl_to_hide), feature(chan_disp, cur_cl_to_hide), 36, [0 .447 .741], 'w.');
    end
    hold off
else
    for i = 1:length(cl_to_hide)
        delete(hide_spks_scat_time{i});
    end
end

end