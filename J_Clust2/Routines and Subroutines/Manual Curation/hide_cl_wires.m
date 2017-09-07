function [hide_cl_scat_wires] = hide_cl_wires(unit_pts, edit_cl, val, chan_disp, feature, hide_cl_scat_wires, Channel_Scatter)

%Description: This .m file is called via the 'Hide Cluster' toggle button in the main GUI, and depending on the toggle state, either 1) hides spikes 
%belonging to a paritcular cluster(s) on the 'Time_Scatter' axis by substituting their cluster-associated color(s) with the 'All Spikes' color, or 
%2) reveals cluster(s) that have been hidden
%
%Input: 'unit_pts' = points belonging to all current units,  'edit_cl' = vector containing info on which cluster(s) to hide, 'val' = toggle state
%of 'Hide Cluster' button, 'chan_disp' = current channels displayed on axis, 'feature' = current feature displayed on axis, 'hide_cl_scat_wires' = 
%scatter plot that hides spikes on axis, 'Channel_Scatter' = axis
%
%Output: 'unit_pts' = updated 'unit_pts'
%

cl_to_hide = find(edit_cl);

axes(Channel_Scatter)

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
        if length(chan_disp) > 2
            hide_cl_scat_wires{i} = scatter3(feature(chan_disp(1), cur_cl_to_hide), feature(chan_disp(2), cur_cl_to_hide), feature(chan_disp(3), cur_cl_to_hide), 36, [0 .447 .741], '.');
        else %length(chan_disp) == 2
            hide_cl_scat_wires{i} = scatter(feature(chan_disp(1), cur_cl_to_hide), feature(chan_disp(2), cur_cl_to_hide), 36, [0 .447 .741], '.');
        end
    end
    hold off
else
    for i = 1:length(cl_to_hide)
        delete(hide_cl_scat_wires{i});
    end
end

end