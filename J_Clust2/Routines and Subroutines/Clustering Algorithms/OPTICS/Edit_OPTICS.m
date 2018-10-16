function varargout = Edit_OPTICS(varargin)
% Edit_OPTICS MATLAB code file for Edit_OPTICS.fig
%
%Description: This .m file opens a GUI with two axes each containing the reachability-distance (RD) plot returned by the OPTICS algorithm. The axes 
%on the right contains the default extracted clusters, while the GUI on the left allows the user to manually select clusters directly from the RD plot.
%
% Begin initialization code - DO NOT EDIT
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Edit_OPTICS_OpeningFcn, ...
                   'gui_OutputFcn',  @Edit_OPTICS_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Edit_OPTICS is made visible.
function Edit_OPTICS_OpeningFcn(edit_optics_obj, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% edit_optics_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Edit_OPTICS (see VARARGIN)

% Choose default command line output for Edit_OPTICS
handles.output = edit_optics_obj;
set(edit_optics_obj,'menubar','figure')

% Import data from main gui and update handles structure
J_Clust_tag = findobj('Tag','J_Clust_fig');
data_from_main = guidata(J_Clust_tag);

handles.unit_pts = data_from_main.unit_pts;
handles.RD = data_from_main.RD;
handles.order = data_from_main.optics_order;
handles.g_max = max(findpeaks(handles.RD(handles.order)));
handles.colors = data_from_main.colors;
handles.new_unit_pts = [];

axes(handles.view_cl_RD_ax)
bar(handles.RD(handles.order)) %bar graph of RD
title('RD Plot of Current Clusters')
ylabel('Mahalanobis Distance^2')

hold on
for i = 1:length(handles.unit_pts)
    [~, order_indxs, ~] = intersect(handles.order,handles.unit_pts{i});
    scatter(order_indxs, (ones(length(handles.unit_pts{i}),1) * handles.g_max), 'Marker', 'x', 'MarkerEdgeColor', handles.colors(i,:));
end
hold off

warning 'off'
legend({'RD bar graph', 'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5', 'Cluster6', 'Cluster7', 'Cluster8', 'Cluster9', 'Cluster10', 'Cluster11', 'Cluster12'}, 'Location', 'NorthWest', 'FontSize', 8);
warning 'on'

axes(handles.edit_cl_RD_ax)
bar(handles.RD(handles.order))
title('RD Plot to Edit Clusters')


guidata(edit_optics_obj, handles);

% UIWAIT makes Edit_OPTICS wait for user response (see UIRESUME)
% uiwait(handles.edit_optics_fig);


% --- Outputs from this function are returned to the command line.
function varargout = Edit_OPTICS_OutputFcn(edit_optics_obj, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% edit_optics_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in add_cl.
function add_cl_Callback(edit_optics_obj, eventdata, handles)
% edit_optics_obj    handle to add_cl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.edit_cl_RD_ax)
hold on

[x,~] = ginput(2);
x = round(x);
if x(1) < 0
    x(1) = 0;
end
if x(2) > length(handles.RD)
    x(2) = length(handles.RD);
end

m = length(handles.new_unit_pts) + 1;
handles.new_unit_pts{m} = handles.order(x(1):x(2));
n = length(handles.new_unit_pts{m});
[~, order_indxs, ~] = intersect(handles.order,handles.new_unit_pts{m});

field = ['unit', num2str(m)];

handles.new_scat.(field) = scatter(order_indxs, (ones(n,1)*handles.g_max), 'Marker', 'x', 'MarkerEdgeColor', handles.colors(m,:)); 
hold off

leg_name = [];
for i = 1:m
    leg_name = [leg_name, ['Cluster', num2str(i)]];
end
warning 'off'
legend({'RD bar graph', 'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5', 'Cluster6', 'Cluster7', 'Cluster8', 'Cluster9', 'Cluster10', 'Cluster11', 'Cluster12'}, 'Location', 'NorthWest', 'FontSize', 8);
warning 'on'

guidata(edit_optics_obj, handles);

% --- Executes on button press in delete_cl.
function delete_cl_Callback(edit_optics_obj, eventdata, handles)
% edit_optics_obj    handle to delete_cl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

del_cl = inputdlg('Enter cluster number you wish to delete in order of creation','Delete a cluster from RD Plot',[1 45]);
if ~isempty(del_cl)
    del_cl = str2num(cell2mat(del_cl));
else
    return;
end
handles.new_unit_pts{del_cl} = [];
handles.new_unit_pts = handles.new_unit_pts(~cellfun('isempty', handles.new_unit_pts));
field = ['unit', num2str(del_cl)];
delete(handles.new_scat.(field));

warning 'off'
legend({'RD bar graph', 'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5', 'Cluster6', 'Cluster7', 'Cluster8', 'Cluster9', 'Cluster10', 'Cluster11', 'Cluster12'}, 'Location', 'NorthWest', 'FontSize', 8);
warning 'on'

guidata(edit_optics_obj, handles);

% --- Executes on button press in view_all_cl.
function view_all_cl_Callback(edit_optics_obj, eventdata, handles)
% edit_optics_obj    handle to view_all_cl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

J_Clust_tag = findobj('Tag','J_Clust_fig');
data_from_main = guidata(J_Clust_tag);

axes(data_from_main.Channel_Scatter)
if length(data_from_main.chan_disp1) > 2
    scatter3(data_from_main.feature_wires(data_from_main.chan_disp1(1),:), data_from_main.feature_wires(data_from_main.chan_disp1(2),:), data_from_main.feature_wires(data_from_main.chan_disp1(3),:), '.');
    hold on
    for i = 1:length(handles.new_unit_pts)
        scatter3(data_from_main.feature_wires(data_from_main.chan_disp1(1), handles.new_unit_pts{i}), data_from_main.feature_wires(data_from_main.chan_disp1(2), handles.new_unit_pts{i}), data_from_main.feature_wires(data_from_main.chan_disp1(3), handles.new_unit_pts{i}), 'Marker', '.', 'MarkerEdgeColor', handles.colors(i,:));
    end
    hold off
    zlabel(num2str(data_from_main.chan_disp1(3)))
else %length(chan_disp) == 2
    scatter(data_from_main.feature_wires(data_from_main.chan_disp1(1),:), data_from_main.feature_wires(data_from_main.chan_disp1(2),:), '.');
    hold on
    for i = 1:length(handles.new_unit_pts)
        scatter(data_from_main.feature_wires(data_from_main.chan_disp1(1), handles.new_unit_pts{i}), data_from_main.feature_wires(data_from_main.chan_disp1(2), handles.new_unit_pts{i}), 'Marker', '.', 'MarkerEdgeColor', handles.colors(i,:));
    end
    hold off
end
xlabel(num2str(data_from_main.chan_disp1(1)))
ylabel(num2str(data_from_main.chan_disp1(2)))
title(data_from_main.feat_disp1)

%Commented out below to currently only display OPTICS selections on channel_scatter

% no_time_disp = 0;
% if strcmp(data_from_main.feat_disp2, 'X vs. Y') || strcmp(data_from_main.feat_disp2, 'View All Waveforms')
%     no_time_disp = 1;
% end
% if ~no_time_disp
%     axes(data_from_main.Time_Scatter)
%     scatter(data_from_main.ts, data_from_main.feature_time(data_from_main.chan_disp2, :), '.')
%     hold on
%     for i = 1:length(handles.new_unit_pts)
%         scatter(data_from_main.ts(handles.new_unit_pts{i}), data_from_main.feature_time(data_from_main.chan_disp2, handles.new_unit_pts{i}), 'Marker', '.', 'MarkerEdgeColor', handles.colors(i,:));
%     end
%     hold off
% end
% 
% xlabel(num2str(data_from_main.chan_disp1(1)))
% title(data_from_main.feat_disp2)

guidata(edit_optics_obj, handles);


% --- Executes on button press in save_all_cl.
function save_all_cl_Callback(edit_optics_obj, eventdata, handles)
% edit_optics_obj    handle to save_all_cl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

save_cl = questdlg('Are you sure you want to save your created clusters? This will overwrite any previous clusters.', ...
    'Confirm Save Clusters');

if strcmp(save_cl, 'Yes')
    handles.new_unit_pts = handles.new_unit_pts(~cellfun('isempty', handles.new_unit_pts));
    global new_unit_pts
    new_unit_pts = handles.new_unit_pts;
end

guidata(edit_optics_obj, handles);

close(Edit_OPTICS);
