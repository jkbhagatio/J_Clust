function polygon_indxs = use_poly_tool(chan_disp1, chan_disp2, feature_wires, feature_time, Channel_Scatter, Time_Scatter, ts)

%Description: This .m file opens an interactive polygon tool the user can use to select spikes on either of the main GUI axes
%
%Input: 'chan_disp1' = channels displayed on 'Channel_Scatter' axis, 'chan_disp2' = channel displayed on 'Time_Scatter' axis, 
%'feature_wires' = feature displayed on 'Channel_Scatter' axis, 'feature_time' = feature displayed on 'Time_Scatter' axis,
%'Channel_Scatter' = axis, 'Time_Scatter' = axis, 'ts' = spike timestamps
%
%Output: 'polygon_indxs' = points located within the user-created polygon
%

select_plot = questdlg('On which plot do you want to use the polygon tool?', 'Select plot for creating polygon', 'Feature Plot', 'Time Plot', 'Cancel', 'Feature Plot');

if isempty(select_plot)
    polygon_indxs = [];
    return;
end

switch select_plot
    case 'Cancel'
        polygon_indxs = [];
        return;
    case 'Feature Plot'
        h_poly = impoly(Channel_Scatter);
        k = waitforbuttonpress;
        while k == 0
            k = waitforbuttonpress;
        end
        if k == 1
            pos = h_poly.getPosition;
            polygon_indxs = find(inpolygon(feature_wires(chan_disp1(1),:), feature_wires(chan_disp1(2), :), pos(:,1), pos(:,2)));
            hold on, scatter(feature_wires(chan_disp1(1), polygon_indxs), feature_wires(chan_disp1(2), polygon_indxs), 'k.'), hold off
            delete(h_poly)
        end
    case 'Time Plot'
        h_poly = impoly(Time_Scatter);
        k = waitforbuttonpress;
        while k == 0
            k = waitforbuttonpress;
        end
        if k == 1
            pos = h_poly.getPosition;
            polygon_indxs = find(inpolygon(ts, feature_time(chan_disp2, :), pos(:,1), pos(:,2)));
            hold on, scatter(ts(polygon_indxs), feature_time(chan_disp2, polygon_indxs), 'k.'), hold off
            delete(h_poly)
        end
end


end