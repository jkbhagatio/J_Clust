function varargout = Manual_Cluster_3D(varargin)
% Manual_Cluster_3D MATLAB code file for Manual_Cluster_3D.fig
%
%Description: This .m file opens a GUI with four axes containing the currently selected spike feature on the 'Channel_Scatter' axis of the main GUI. 
%Each axis contains one of the four possible triples of combinations of tetrode wires. The user can add a cluster by selecting vertices of the 
%putative cluster in the 3-D space of any of the axes via the 'brush' tool.
%
% Begin initialization code - DO NOT EDIT
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Manual_Cluster_3D_OpeningFcn, ...
                   'gui_OutputFcn',  @Manual_Cluster_3D_OutputFcn, ...
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


% --- Executes just before Manual_Cluster_3D is made visible.
function Manual_Cluster_3D_OpeningFcn(man_cluster_obj, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% man_cluster_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Manual_Cluster_3D (see VARARGIN)

% Choose default command line output for Manual_Cluster_3D
handles.output = man_cluster_obj;
set(man_cluster_obj,'menubar','figure')

% Import data from main gui and update handles structure
J_Clust_tag = findobj('Tag','J_Clust_fig');
data_from_main = guidata(J_Clust_tag);

handles.feature = data_from_main.feature_wires;
handles.feat_disp = data_from_main.feat_disp1;
handles.colors = data_from_main.colors;
handles.unit_pts = [];
handles.cur_unit = 1;

axes(handles.axes1)
scatter3(handles.feature(1,:), handles.feature(2,:), handles.feature(3,:), '.');
title(handles.feat_disp)
xlabel('Wire 1')
ylabel('Wire 2')
zlabel('Wire 3')

axes(handles.axes2)
scatter3(handles.feature(1,:), handles.feature(2,:), handles.feature(4,:), '.');
title(handles.feat_disp)
xlabel('Wire 1')
ylabel('Wire 3')
zlabel('Wire 4')

axes(handles.axes3)
scatter3(handles.feature(1,:), handles.feature(3,:), handles.feature(4,:), '.');
title(handles.feat_disp)
xlabel('Wire 1')
ylabel('Wire 3')
zlabel('Wire 4')

axes(handles.axes4)
scatter3(handles.feature(2,:), handles.feature(3,:), handles.feature(4,:), '.');
title(handles.feat_disp)
xlabel('Wire 2')
ylabel('Wire 3')
zlabel('Wire 4')

guidata(man_cluster_obj, handles);

% UIWAIT makes Manual_Cluster_3D wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Manual_Cluster_3D_OutputFcn(man_cluster_obj, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% man_cluster_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in add_cluster.
function add_cluster_Callback(man_cluster_obj, eventdata, handles)
% man_cluster_obj    handle to add_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ax = gca;
fig = gcf;

scat = get(ax, 'Children');
scat = scat(length(scat)); %original scatter with all points
pts = [scat.XData; scat.YData; scat.ZData];

%% Get Vertices

% fin = false;
h_brush = brush(fig);
set(h_brush,'Color',handles.colors(handles.cur_unit,:),'Enable','on');
vertices = [];
k = 0;
hold (ax, 'on')
while k ~= 2
    
    k = waitforbuttonpress;
    brushed_indxs = find(get(scat, 'BrushData'));
    
    delete(findobj(fig, 'type', 'legend'))
    scatter3(pts(1,brushed_indxs), pts(2,brushed_indxs), pts(3,brushed_indxs),'Marker', 'o', 'MarkerEdgeColor', handles.colors(handles.cur_unit,:), 'SizeData', 80);
    
    if strcmp(get(fig,'SelectionType'),'alt')
        scat_objs = get(ax, 'Children');
        delete(scat_objs(1:length(scat_objs)-1));
        brushed_indxs = [];
        set(scat, 'BrushData', 0);
        
        if handles.cur_unit > 1 %re-scatter old unit points if clustering multiple units
            for i = 1:handles.cur_unit
                field = ['unit', num2str(i)];
                hold(ax, 'on')
                cur_unit_pts = handles.unit_pts{i};
                scatter3(pts(1, cur_unit_pts), pts(2, cur_unit_pts), pts(3, cur_unit_pts), 'Marker', '.', 'MarkerFaceColor', handles.colors(i,:), 'MarkerEdgeColor', handles.colors(i,:));
            end
        end
        
    end
    
    if strcmp(get(fig,'SelectionType'),'open')
        vertices = pts(:, brushed_indxs);
        if length(vertices) > 4
            set(scat, 'BrushData', 0);
            k = 2;
        else
            fprintf('\nYou must select at least 5 vertices');
        end     
    end
    
end

hold (ax, 'off')

% while ~fin
%     %Choose vertex points using brush tool
%     
%     if strcmp(get(fig,'SelectionType'),'open')
%         brushed_indxs = find(get(scat, 'BrushData'));
%         hold(ax, 'on')
%         brush_scat = scatter3(pts(1,brushed_indxs), pts(2,brushed_indxs), pts(3,brushed_indxs),'Marker','o', 'MarkerEdgeColor', handles.colors(handles.cur_unit,:), 'SizeData', 80);
%         hold(ax, 'off')
%         
%         set(scat, 'BrushData', 0);
%         vertices = pts(:,brushed_indxs);
%         
%         if length(vertices) < 4
%             fprintf('\nMust select at least 4 vertices for clustering\n')
%             delete(brush_scat);
%         end
%     end
%     
%     if strcmp(get(fig,'SelectionType'),'extend')
%         delete(brush_scat);
%         brushed_indxs = [];
%         vertices = [];
%     end
%     
%     if strcmp(get(fig,'SelectionType'),'alt')
%         if length(vertices) > 4
%             fin = true;
%         end
%     end
%     
% end

%% Create cluster from vertices

tetra = delaunayn(vertices'); % Generate delaunay triangulization
faces = freeBoundary(triangulation(tetra,vertices')); % use free boundary as triangulation
handles.unit_pts{handles.cur_unit} = intriangulation(vertices',faces,pts');
handles.unit_pts{handles.cur_unit} = [find(handles.unit_pts{handles.cur_unit}); brushed_indxs'];

%re-plot all clusters on all axes
for x = 1:4 %where x = each axis
    switch x
        case 1
            chans = [1 2 3];
        case 2 
            chans = [1 2 4];
        case 3
            chans = [1 3 4];
        case 4
            chans = [2 3 4];
    end
       
    field_ax = ['axes',num2str(x)]; %the current axis
    axes(handles.(field_ax))
    scatter3(handles.feature(chans(1),:), handles.feature(chans(2),:), handles.feature(chans(3),:), '.');
    scat = get(handles.(field_ax), 'Children');
    pts = [scat.XData; scat.YData; scat.ZData];
    
    hold on

    for i = 1:handles.cur_unit
        field_unit = ['unit', num2str(i)];
        handles.unit_scatter.(field_ax).(field_unit) = scatter3(pts(1, handles.unit_pts{i}),pts(2, handles.unit_pts{i}), pts(3, handles.unit_pts{i}),'Marker', '.', 'MarkerFaceColor', handles.colors(handles.cur_unit,:), 'MarkerEdgeColor', handles.colors(i,:));
    end
    hold off
    warning 'off'
    handles.leg.(field_ax) = legend({'All Spikes', 'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4', 'Cluster5', 'Cluster6', 'Cluster7', 'Cluster8', 'Cluster9', 'Cluster10', 'Cluster11', 'Cluster12'}, 'Location', 'SouthWest', 'FontSize', 8);
    warning 'on'
end

handles.cur_unit = handles.cur_unit + 1;


guidata(man_cluster_obj, handles);           


% --- Executes on button press in delete_cluster.
function delete_cluster_Callback(man_cluster_obj, eventdata, handles)
% man_cluster_obj    handle to delete_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

del_cl = str2num(cell2mat(inputdlg('Enter number of cluster you wish to delete: ', 'Delete Cluster', [1 50])));

handles.unit_pts{del_cl} = [];

for x = 1:4 %where x = each axis
    field_ax = ['axes', num2str(x)];
    field_unit = ['unit', num2str(del_cl)];
    
    delete(handles.unit_scatter.(field_ax).(field_unit));
end


guidata(man_cluster_obj, handles);


% --- Executes on button press in save_clusters.
function save_clusters_Callback(man_cluster_obj, eventdata, handles)
% man_cluster_obj    handle to save_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

save_cl = questdlg('Are you satisfied with these clusters?', 'Confirm Save Clusters');

if strcmp(save_cl, 'Yes')
    global new_unit_pts
    new_unit_pts = handles.unit_pts;
end

guidata(man_cluster_obj, handles);

close(Manual_Cluster_3D)
