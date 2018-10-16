function selected_pts = get_selected_pts()

%Description: This .m file gets selected points (on either of the main GUI axes) via 'brush' and '2-D Poly Tool'
%
%Ouput: 'selected_pts' = spikes selected from main GUI axes

J_Clust_tag = findobj('Tag','J_Clust_fig');
data_from_main = guidata(J_Clust_tag);
polygon_indxs = data_from_main.polygon_indxs;
Channel_Scatter = data_from_main.Channel_Scatter;
Time_Scatter = data_from_main.Time_Scatter;
unit_pts = data_from_main.unit_pts;
view_1st_cl = find(data_from_main.view_1st_cl);
edit_cl = find(data_from_main.edit_cl);

comp_metric_indxs = [];
single_metric_indxs = [];

comp_metrics_tag = findobj('Tag', 'comparative_metrics_fig');
if ~isempty(comp_metrics_tag)
    data_from_comp = guidata(comp_metrics_tag);
    if isfield(data_from_comp, 'remove_selected_indxs') && ~isempty(data_from_comp.remove_selected_indxs)
        unit_pts_cur = unit_pts{view_1st_cl};
        comp_metric_indxs = unit_pts_cur(data_from_comp.remove_selected_indxs);
    end
end

single_metrics_tag = findobj('Tag', 'single_metrics_fig');
if ~isempty(single_metrics_tag)
    data_from_sing = guidata(single_metrics_tag);
    if isfield(data_from_sing, 'remove_selected_indxs')
        unit_pts_cur = unit_pts{view_1st_cl};
        single_metric_indxs = unit_pts_cur(data_from_sing.remove_selected_indxs);
        comp_metric_indxs = []; %in case 'single_metric_indxs' overlap with 'comp_metric_indxs' and there's a subset of the latter we do NOT want to remove
    end
end

child_scatters_wires = get(Channel_Scatter, 'Children');
brushed_indxs_wires = find(get(child_scatters_wires(length(child_scatters_wires)), 'BrushData'));

child_scatters_time = get(Time_Scatter, 'Children');
brushed_indxs_time = find(get(child_scatters_time(length(child_scatters_time)), 'BrushData'));

selected_pts = unique(horzcat(brushed_indxs_wires, brushed_indxs_time, polygon_indxs, comp_metric_indxs, single_metric_indxs));
%selected_pts = horzcat(brushed_indxs_wires, brushed_indxs_time, polygon_indxs];

end