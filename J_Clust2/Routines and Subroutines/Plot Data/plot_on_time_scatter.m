%Plot_on_Time_Scatter
%
%Description: This script is called to save a few lines of code when updating the 'Time_Scatter' axis in the main GUI

axes(handles.Time_Scatter)
if handles.preload
    if isempty(handles.unit_pts)
        handles.feature_time = time_scatter_plot(handles.chan_disp2, handles.feat_disp2, handles.features, handles.cur_ts);
    else %some units have already been defined
        handles.feature_time = time_scatter_plot(handles.chan_disp2, handles.feat_disp2, handles.features, handles.cur_ts, handles.unit_pts, handles.colors);
    end
else
    if isempty(handles.unit_pts)
        handles.feature_time = time_scatter_plot(handles.chan_disp2, handles.feat_disp2, handles.features, handles.ts);
    else %some units have already been defined
        handles.feature_time = time_scatter_plot(handles.chan_disp2, handles.feat_disp2, handles.features, handles.ts, handles.unit_pts, handles.colors);
    end
end