%Plot_on_Channel_Scatter
%
%Description: This script is called to save a few lines of code when updating the 'Channel_Scatter' axis in the main GUI

axes(handles.Channel_Scatter)
if isempty(handles.unit_pts)
    handles.feature_wires = channel_scatter_plot(handles.chan_disp1, handles.feat_disp1, handles.features);
else %some units have already been defined
    handles.feature_wires = channel_scatter_plot(handles.chan_disp1, handles.feat_disp1, handles.features, handles.unit_pts, handles.colors);
end