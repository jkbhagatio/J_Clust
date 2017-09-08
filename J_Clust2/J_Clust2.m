function varargout = J_Clust2(varargin)
%J_Clust2 MATLAB code file for J_Clust2.fig
%
%Description: This .m file opens the main GUI for J_Clust2. All nested functions in this file are callback functions of the main GUI, or
%subfunctions within the software package J_Clust2.
%
%Note: The command, 'J_Clust2', entered in the MATLAB Command Window, launches the program.
%
% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @J_Clust2_OpeningFcn, ...
                   'gui_OutputFcn',  @J_Clust2_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


function J_Clust2_OpeningFcn(J_Clust_obj, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% J_Clust_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
% --- Executes just before J_Clust2 is made visible.

% Choose default command line output for J_Clust2
handles.output = J_Clust_obj;


% Update handles structure
guidata(J_Clust_obj, handles);
set(J_Clust_obj,'menubar','figure')


function varargout = J_Clust2_OutputFcn(J_Clust_obj, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% J_Clust_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Outputs from this function are returned to the command line.

% Get default command line output from handles structure
varargout{1} = handles.output;

function load_data_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on selection change in load_data.
% Hints: contents = cellstr(get(J_Clust_obj,'String')) returns load_data contents as cell array
%        contents{get(J_Clust_obj,'Value')} returns selected item from load_data

%get file info
load_data_contents = get(J_Clust_obj, 'String');
load_data_val = get(J_Clust_obj, 'Value');

[handles.Fs, handles.uV_conversion, var_outs]  = load_sig(load_data_contents, load_data_val); %type of data loaded: (raw signal, filtered signal, or spike waveforms & timestamps)
if length(var_outs) < 2 %raw signal loaded into GUI
    handles.filt_sig = var_outs{1};
    handles.preload = 0;
else %waveforms preloaded into GUI
    handles.waveforms = var_outs{1}; handles.ts = var_outs{2}; handles.num_spks = var_outs{3}; handles.new_spk_mrkr = var_outs{4}; handles.num_samples = var_outs{5}; handles.overlaps = var_outs{6}; 
    handles.filt_sig = 0;
    handles.preload = 1;
end

%-Variable Initialization

%initialize clusters ('unit_pts') and cluster colors
handles.unit_pts = [];
handles.threshold = [];
handles.colors = [1 0 0; 0 1 0; 1 0 1; 1 1 0; .6 .6 .4; .7 .3 .3; .7 .7 .3; .55 .45 .55; .85 .85 .15; .5 0 .5; 0 .5 .5; .25 .75 .25; 0 0 0; 0 .447 .741]; 

%set colors for clusters in cluster menu
set(handles.all_spikes_cluster_menu, 'backgroundcolor', handles.colors(14,:))
set(handles.brush_spikes_cluster_menu, 'backgroundcolor', 'k')
set(handles.brush_spikes_cluster_menu, 'foregroundcolor', [1 1 1]);

for i = 1:12
    var_name = ['cl_', num2str(i), '_cluster_menu'];
    set(handles.(var_name), 'backgroundcolor', handles.colors(i,:))
end

%initialize cluster menu checkboxes and cluster visualization radio buttons
handles.edit_cl = zeros(14,1);
handles.view_1st_cl = zeros(14,1);
handles.view_2nd_cl = zeros(14,1);
handles.epochs = zeros(2,10);

handles.hide_cl_scat_wires = cell(14,1);
handles.hide_spks_scat_wires = cell(14,1);
handles.hide_cl_scat_time = cell(14,1);
handles.hide_spks_scat_time = cell(14,1);
handles.polygon_indxs = [];

handles.waveforms_select = [];
handles.isi_hist_select = [];
handles.autocorr_select = [];
handles.residuals_select = [];
handles.chebyshev_select = [];
handles.silhouette_select = [];
handles.crosscorr_select = [];
handles.m_distance_select = [];
handles.lratio_select = [];

handles.output_ts = [];
handles.output_waveforms = [];
handles.output_isis = [];
handles.output_features = [];
handles.output_lratio = [];

%initialize variables for hotkey control of plots
handles.chan_combo_val = 3;
handles.feat_wires_val = 3;
handles.chan_time_val = 3;
handles.feat_time_val = 3;

handles.chan_combo_contents = {'Channel Combo', [], '1,2,3', '1,2,4', '1,3,4', '2,3,4', '1,2', '1,3', '1,4', '2,3', '2,4', '3,4'};
handles.feat_wires_contents = {'Feature', [], 'Peak Amplitude','Crest-to-Trough Amplitude','Power','PCA Scores','PCA Scores for Concatenated Waveforms','Wavelet Coefficients', 'Wavelet Coefficients for Concatenated Waveforms'};
handles.chan_time_contents = {'Channel', [], '1', '2', '3', '4'};
handles.feat_time_contents = {'Feature', [], 'Peak Amplitude','Crest-to-Trough Amplitude','Power','PCA Scores','PCA Scores for Concatenated Waveforms','Wavelet Coefficients', 'Wavelet Coefficients for Concatenated Waveforms', 'X Vs. Y'};


guidata(J_Clust_obj, handles);

function load_data_CreateFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(J_Clust_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(J_Clust_obj,'BackgroundColor','white');
end

function set_time_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to set_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in set_time.

set_time_title = 'Set Time';
set_time_prompt = {'Enter Start Time (in seconds):','Enter End Time (in seconds):'};
set_time_info = inputdlg(set_time_prompt, set_time_title, [1 40]);

if isempty(set_time_info)
    return;
end

handles.start_time = str2num(set_time_info{1}); 
handles.end_time = str2num(set_time_info{2});

%allow for adjusting the visualization of clustered units when adjusting the time display
if ~isempty(handles.unit_pts)
    if handles.preload
        last_spike_ts = max(handles.ts(handles.first_spk:handles.last_spk));
        if handles.end_time < last_spike_ts
            for i = 1:length(handles.unit_pts)
                cur_unit_pts = handles.unit_pts{i};
                [~, close_indx] = min(abs(handles.ts(cur_unit_pts) - handles.end_time));
                cur_unit_pts(close_indx-1:end) = [];
                handles.unit_pts{i} = cur_unit_pts;
            end
        end
    else
        last_spike_ts = max(handles.ts);
        if handles.end_time < last_spike_ts
            for i = 1:length(handles.unit_pts)
                cur_unit_pts = handles.unit_pts{i};
                [~, close_indx] = min(abs(handles.ts(cur_unit_pts) - handles.end_time));
                cur_unit_pts(close_indx-1:end) = [];
                handles.unit_pts{i} = cur_unit_pts;
            end
        end
    end
end

initialize_plotting; %calculate spike features and load initial plots on axes

if handles.preload
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts(handles.first_spk:handles.last_spk), handles.overlaps, handles.cluster_info);
else
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts, handles.overlaps, handles.cluster_info);
end


guidata(J_Clust_obj, handles);


function set_time_CreateFcn(J_Clust_obj, eventdata, handles, features)
% J_Clust_obj    handle to set_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.


function set_view_epoch_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to set_view_epoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on selection change in set_view_epoch.
% Hints: contents = cellstr(get(J_Clust_obj,'String')) returns set_view_epoch contents as cell array
%        contents{get(J_Clust_obj,'Value')} returns selected item from set_view_epoch

epoch_contents = get(J_Clust_obj, 'String');
epoch_val = get(J_Clust_obj, 'Value') - 2;

if epoch_val == -1 || epoch_val == 0 || epoch_val == 11
    return;
end

[handles.epochs, handles.break_set_epoch] = set_epoch(handles.epochs, handles.start_time, handles.end_time, epoch_val);
if handles.break_set_epoch
    return;
end
handles.start_time = handles.epochs(1, epoch_val);
handles.end_time = handles.epochs(2, epoch_val);

initialize_plotting;

if handles.preload
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts(handles.first_spk:handles.last_spk), handles.overlaps, handles.cluster_info);
else
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts, handles.overlaps, handles.cluster_info);
end

guidata(J_Clust_obj, handles);


function set_view_epoch_CreateFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to set_view_epoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(J_Clust_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(J_Clust_obj,'BackgroundColor','white');
end


function channel_combo_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to channel_combo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on selection change in channel_combo.

% Hints: contents = cellstr(get(J_Clust_obj,'String')) returns channel_combo contents as cell array
%        contents{get(J_Clust_obj,'Value')} returns selected item from channel_combo

handles.chan_combo_contents = get(J_Clust_obj, 'String');
handles.chan_combo_val = get(J_Clust_obj, 'Value');

if handles.chan_combo_val == 1 || handles.chan_combo_val == 2
    return;
end

handles.chan_disp1 = str2num(handles.chan_combo_contents{handles.chan_combo_val});

plot_on_channel_scatter;

handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection


guidata(J_Clust_obj, handles);


function channel_combo_CreateFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to channel_combo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(J_Clust_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(J_Clust_obj,'BackgroundColor','white');
end


function feature_wires_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to feature_wires (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on selection change in feature_wires.
% Hints: contents = cellstr(get(J_Clust_obj,'String')) returns feature_wires contents as cell array
%        contents{get(J_Clust_obj,'Value')} returns selected item from feature_wires

handles.feat_wires_contents = get(J_Clust_obj, 'String');
handles.feat_wires_val = get(J_Clust_obj, 'Value');

if handles.feat_wires_val == 1 || handles.feat_wires_val == 2
    return;
end

handles.feat_disp1 = handles.feat_wires_contents{handles.feat_wires_val};

plot_on_channel_scatter;

handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection


guidata(J_Clust_obj, handles);


function feature_wires_CreateFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to feature_wires (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(J_Clust_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(J_Clust_obj,'BackgroundColor','white');
end

function channel_time_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to channel_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on selection change in channel_time.
% Hints: contents = cellstr(get(J_Clust_obj,'String')) returns channel_time contents as cell array
%        contents{get(J_Clust_obj,'Value')} returns selected item from channel_time

handles.chan_time_contents = get(J_Clust_obj, 'String');
handles.chan_time_val = get(J_Clust_obj, 'Value');

if handles.chan_time_val == 1 || handles.chan_time_val == 2
    return;
end

handles.chan_disp2 = str2num(handles.chan_time_contents{handles.chan_time_val});

plot_on_time_scatter;

handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection


guidata(J_Clust_obj, handles);


function channel_time_CreateFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to channel_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(J_Clust_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(J_Clust_obj,'BackgroundColor','white');
end


function feature_time_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to feature_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on selection change in feature_time.

% Hints: contents = cellstr(get(J_Clust_obj,'String')) returns feature_time contents as cell array
%        contents{get(J_Clust_obj,'Value')} returns selected item from feature_time

handles.feat_time_contents = get(J_Clust_obj, 'String');
handles.feat_time_val = get(J_Clust_obj, 'Value');

if handles.feat_time_val == 1 || handles.feat_time_val == 2
    return;
end

handles.feat_disp2 = handles.feat_time_contents{handles.feat_time_val};

plot_on_time_scatter;

handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for new selection


guidata(J_Clust_obj, handles);


function feature_time_CreateFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to feature_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(J_Clust_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(J_Clust_obj,'BackgroundColor','white');
end


function Channel_Scatter_CreateFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to Channel_Scatter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.

% Hint: place code in OpeningFcn to populate Channel_Scatter


function Time_Scatter_CreateFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to Time_Scatter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.

% Hint: place code in OpeningFcn to populate Time_Scatter

% --- Executes on selection change in select_coeff_menu.
function select_coeff_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to select_coeff_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(J_Clust_obj,'String')) returns select_coeff_menu contents as cell array
%        contents{get(J_Clust_obj,'Value')} returns selected item from select_coeff_menu

select_coeff_contents = get(J_Clust_obj, 'String');
select_coeff_val = get(J_Clust_obj, 'Value');

if select_coeff_val == 1 || select_coeff_val == 2
    return;
end

cur_select_coeff = select_coeff_contents(select_coeff_val);

if strcmp(cur_select_coeff, 'PC Select')
    Select_Coefficient_PCs;
end

if strcmp(cur_select_coeff, 'Concatenated PC Select')
    Select_Coefficient_PCs_c;
end

% --- Executes during object creation, after setting all properties.
function select_coeff_menu_CreateFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to select_coeff_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(J_Clust_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(J_Clust_obj,'BackgroundColor','white');
end


function run_fcm_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to run_fcm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in run_fcm.

%acquire units from fcm
N_c = inputdlg('Enter number of clusters (1-12):', 'Cluster Number', [1 45]); %number of clusters

if isempty(N_c)
    return;
else
    N_c = str2num(cell2mat(N_c));
end

disp('Clustering...')
handles.unit_pts = run_fcm(handles.feature_wires, N_c);

for i = 1:length(handles.unit_pts)
    if ~isrow(handles.unit_pts{i})
        handles.unit_pts{i} = handles.unit_pts{i}';
    end
end

%update rest of GUI
plot_on_channel_scatter;
plot_on_time_scatter;
if handles.preload
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts(handles.first_spk:handles.last_spk), handles.overlaps, handles.cluster_info);
else
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts, handles.overlaps, handles.cluster_info);
end
disp('Done')


guidata(J_Clust_obj, handles);


function run_spc_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to run_spc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in run_spc.


function run_optics_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to run_optics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in run_optics.

min_pts = inputdlg('Enter minimum number of points in a cluster (as an integer):', 'Min Pts in Cluster', [1 45]);

if isempty(min_pts)
    return;
else
    min_pts = str2num(cell2mat(min_pts));
end

disp('Clustering...')

[handles.RD, CD, handles.optics_order] = optics_JClust(handles.feature_wires', min_pts);
% figure, bar(handles.RD(handles.optics_order));
handles.unit_pts = extract_clusters_optics(handles.RD, handles.optics_order, min_pts);
guidata(J_Clust_obj, handles);

Edit_OPTICS;
uiwait(Edit_OPTICS)

%if new cluster(s) based on Edit_OPTICS
global new_unit_pts
if ~isempty(new_unit_pts)
    handles.unit_pts = new_unit_pts;
end
clear global new_unit_pts
guidata(J_Clust_obj, handles);

for i = 1:length(handles.unit_pts)
    if ~isrow(handles.unit_pts{i})
        handles.unit_pts{i} = handles.unit_pts{i}';
    end
end

%update rest of GUI
plot_on_channel_scatter;
plot_on_time_scatter;
if handles.preload
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts(handles.first_spk:handles.last_spk), handles.overlaps, handles.cluster_info);
else
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts, handles.overlaps, handles.cluster_info);
end
disp('Done')


guidata(J_Clust_obj, handles);


function manual_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in manual.

Manual_Cluster_3D;
uiwait(Manual_Cluster_3D)

global new_unit_pts
handles.unit_pts = new_unit_pts;
clear global new_unit_pts
guidata(J_Clust_obj, handles);

for i = 1:length(handles.unit_pts)
    if ~isrow(handles.unit_pts{i})
        handles.unit_pts{i} = handles.unit_pts{i}';
    end
end

%update rest of GUI
plot_on_channel_scatter;
plot_on_time_scatter;
if handles.preload
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts(handles.first_spk:handles.last_spk), handles.overlaps, handles.cluster_info);
else
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts, handles.overlaps, handles.cluster_info);
end

guidata(J_Clust_obj, handles);


function poly_tool_2d_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to poly_tool_2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in poly_tool_2d.

if handles.preload
    handles.polygon_indxs = use_poly_tool(handles.chan_disp1, handles.chan_disp2, handles.feature_wires, handles.feature_time, handles.Channel_Scatter, handles.Time_Scatter, handles.ts(handles.first_spk:handles.last_spk));
else
    handles.polygon_indxs = use_poly_tool(handles.chan_disp1, handles.chan_disp2, handles.feature_wires, handles.feature_time, handles.Channel_Scatter, handles.Time_Scatter, handles.ts);
end

guidata(J_Clust_obj, handles);


function add_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in add.

handles.selected_pts = get_selected_pts(); %get selected points
handles.unit_pts = add_pts_to_cl(handles.unit_pts, handles.selected_pts, handles.edit_cl); %update unit points
handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for next selection

%update scatter plots and cluster table info
plot_on_channel_scatter;
plot_on_time_scatter;

if handles.preload
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts(handles.first_spk:handles.last_spk), handles.overlaps, handles.cluster_info);
else
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts, handles.overlaps, handles.cluster_info);
end

guidata(J_Clust_obj, handles);


function remove_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in remove.

handles.selected_pts = get_selected_pts(); %get selected points
handles.unit_pts = remove_pts_from_cl(handles.unit_pts, handles.selected_pts, handles.edit_cl);
handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for next selection

%update scatter plots and cluster table info
plot_on_channel_scatter;
plot_on_time_scatter;

if handles.preload
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts(handles.first_spk:handles.last_spk), handles.overlaps, handles.cluster_info);
else
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts, handles.overlaps, handles.cluster_info);
end

guidata(J_Clust_obj, handles);


function merge_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in merge.

handles.selected_pts = get_selected_pts(); %get selected points
handles.unit_pts = merge_cl(handles.unit_pts, handles.edit_cl);
handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for next selection

%update scatter plots and cluster table info
plot_on_channel_scatter;
plot_on_time_scatter;

if handles.preload
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts(handles.first_spk:handles.last_spk), handles.overlaps, handles.cluster_info);
else
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts, handles.overlaps, handles.cluster_info);
end

guidata(J_Clust_obj, handles);


function delete_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in delete.

handles.selected_pts = get_selected_pts(); %get selected points
handles.unit_pts = delete_cl(handles.unit_pts, handles.edit_cl);
handles.selected_pts = []; handles.polygon_indxs = []; %empty vector to reset for next selection

%update scatter plots and cluster table info
plot_on_channel_scatter;
plot_on_time_scatter;

if handles.preload
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts(handles.first_spk:handles.last_spk), handles.overlaps, handles.cluster_info);
else
    handles.cluster_info = update_info_table(handles.unit_pts, handles.ts, handles.overlaps, handles.cluster_info);
end

guidata(J_Clust_obj, handles);


function hide_cluster_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to hide_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in hide_cluster.

val = get(handles.hide_cluster, 'Value');

[handles.hide_cl_scat_wires] = hide_cl_wires(handles.unit_pts, handles.edit_cl, val, handles.chan_disp1, handles.feature_wires, handles.hide_cl_scat_wires, handles.Channel_Scatter);

if handles.preload
    [handles.hide_cl_scat_time] = hide_cl_time(handles.unit_pts, handles.edit_cl, val, handles.chan_disp2, handles.feature_time, handles.hide_cl_scat_time, handles.ts(handles.first_spk:handles.last_spk), handles.Time_Scatter);
else
    [handles.hide_cl_scat_time] = hide_cl_time(handles.unit_pts, handles.edit_cl, val, handles.chan_disp2, handles.feature_time, handles.hide_cl_scat_time, handles.ts, handles.Time_Scatter);
end


guidata(J_Clust_obj, handles);


function hide_spikes_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to hide_spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in hide_spikes.

val = get(handles.hide_spikes, 'Value');

[handles.hide_spks_scat_wires] = hide_spks_wires(handles.unit_pts, handles.edit_cl, val, handles.chan_disp1, handles.feature_wires, handles.hide_spks_scat_wires, handles.Channel_Scatter);

if handles.preload
    [handles.hide_spks_scat_time] = hide_spks_time(handles.unit_pts, handles.edit_cl, val, handles.chan_disp2, handles.feature_time, handles.hide_spks_scat_time, handles.ts(handles.first_spk:handles.last_spk), handles.Time_Scatter);
else
    [handles.hide_spks_scat_time] = hide_spks_time(handles.unit_pts, handles.edit_cl, val, handles.chan_disp2, handles.feature_time, handles.hide_spks_scat_time, handles.ts, handles.Time_Scatter);
end


guidata(J_Clust_obj, handles);


% --- Executes on button press in change_cl_color.
function change_cl_color_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to change_cl_color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if find(max(handles.edit_cl)) > 12
    error('Cannot change color of "all spikes" or "selected spikes"')
end

if length(find(handles.edit_cl)) > 1
    error('Cannot change colors of multiple clusters at once')
end

set_cl = find(handles.edit_cl);
new_color = uisetcolor;

if ~new_color
    return;
end

handles.colors(set_cl, :) = new_color;
var_name = ['cl_', num2str(set_cl), '_cluster_menu'];
set(handles.(var_name), 'backgroundcolor', new_color);

plot_on_channel_scatter;
plot_on_time_scatter;

guidata(J_Clust_obj, handles);


function all_spikes_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to all_spikes_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in all_spikes_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of all_spikes_cluster_menu

handles.edit_cl(14) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function brush_spikes_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to brush_spikes_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in brush_spikes_cluster_menu.

% Hint: get(J_Clust_obj,'Value') returns toggle state of brush_spikes_cluster_menu

handles.edit_cl(13) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_1_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_1_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_1_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_1_cluster_menu

handles.edit_cl(1) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_2_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_2_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_2_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_2_cluster_menu

handles.edit_cl(2) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_3_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_3_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_3_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_3_cluster_menu

handles.edit_cl(3) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_4_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_4_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_4_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_4_cluster_menu

handles.edit_cl(4) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_5_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_5_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_5_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_5_cluster_menu

handles.edit_cl(5) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_6_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_6_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_6_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_6_cluster_menu

handles.edit_cl(6) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_7_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_7_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_7_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_7_cluster_menu

handles.edit_cl(7) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_8_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_8_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_8_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_8_cluster_menu

handles.edit_cl(8) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_9_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_9_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_9_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_9_cluster_menu

handles.edit_cl(9) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_10_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_10_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_10_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_10_cluster_menu

handles.edit_cl(10) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_11_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_11_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_11_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_11_cluster_menu

handles.edit_cl(11) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cl_12_cluster_menu_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cl_12_cluster_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cl_12_cluster_menu.
% Hint: get(J_Clust_obj,'Value') returns toggle state of cl_12_cluster_menu

handles.edit_cl(12) = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function cluster_info_CellSelectionCallback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cluster_info (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selected
% handles    structure with handles and user data (see GUIDATA)
% --- Executes when selected cell(s) is changed in cluster_info.


function cluster_info_CreateFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to cluster_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.


function view_all_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_all_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_all_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_all_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(14) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_brushed_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_brushed_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_brushed_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_brushed_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(13) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_1_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_1_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_1_1st_radio.

% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_1_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(1) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_2_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_2_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_2_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_2_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(2) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_3_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_3_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_3_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_3_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(3) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_4_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_4_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_4_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_4_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(4) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_5_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_5_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_5_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_5_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(5) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_6_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_6_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_6_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_6_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(6) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_7_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_7_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_7_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_7_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(7) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_8_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_8_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_8_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_8_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(8) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_9_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_9_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_9_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_9_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(9) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_10_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_10_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_10_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_10_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(10) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_11_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_11_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_11_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_11_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(11) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_12_1st_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_12_1st_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_12_1st_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_12_1st_radio

handles.view_1st_cl = zeros(14,1); %resets variable before new selection

handles.view_1st_cl(12) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function disp_waveforms_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to disp_waveforms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in disp_waveforms.
% Hint: get(J_Clust_obj,'Value') returns toggle state of disp_waveforms

handles.waveforms_select = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function disp_isi_histogram_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to disp_isi_histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in disp_isi_histogram.
% Hint: get(J_Clust_obj,'Value') returns toggle state of disp_isi_histogram

handles.isi_hist_select = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


% --- Executes on button press in disp_residuals.
function disp_residuals_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to disp_residuals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(J_Clust_obj,'Value') returns toggle state of disp_residuals

handles.residuals_select = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function disp_autocorrelogram_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to disp_autocorrelogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in disp_autocorrelogram.
% Hint: get(J_Clust_obj,'Value') returns toggle state of disp_autocorrelogram

handles.autocorr_select = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function disp_chebyshevs_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to disp_chebyshevs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in disp_chebyshevs.
% Hint: get(J_Clust_obj,'Value') returns toggle state of disp_chebyshevs

handles.chebyshev_select = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function display_single_unit_metrics_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to display_single_unit_metrics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in display_single_unit_metrics.

Single_Unit_Metrics;


function view_all_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_all_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_all_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_all_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(14) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_brushed_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_brushed_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_brushed_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_brushed_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(13) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_1_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_1_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_1_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_1_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(1) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_2_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_2_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_2_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_2_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(2) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_3_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_3_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_3_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_3_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(3) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_4_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_4_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_4_2nd_radio.

% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_4_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(4) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_5_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_5_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_5_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_5_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(5) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_6_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_6_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_6_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_6_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(6) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_7_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_7_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_7_2nd_radio.

% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_7_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(7) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_8_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_8_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_8_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_8_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(8) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_9_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_9_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_9_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_9_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(9) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_10_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_10_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_10_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_10_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(10) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_11_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_11_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_11_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_11_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(11) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function view_cl_12_2nd_radio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to view_cl_12_2nd_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in view_cl_12_2nd_radio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of view_cl_12_2nd_radio

handles.view_2nd_cl = zeros(14,1); %resets variable before new selection

handles.view_2nd_cl(12) = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function disp_m_distance_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to disp_m_distance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in disp_m_distance.
% Hint: get(J_Clust_obj,'Value') returns toggle state of disp_m_distance

handles.m_distance_select = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function disp_crosscorrelogram_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to disp_crosscorrelogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in disp_crosscorrelogram.
% Hint: get(J_Clust_obj,'Value') returns toggle state of disp_crosscorrelogram

handles.crosscorr_select = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function disp_silhouette_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to disp_silhouette (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in disp_silhouette.
% Hint: get(J_Clust_obj,'Value') returns toggle state of disp_silhouette

handles.silhouette_select = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function disp_lratio_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to disp_lratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in disp_lratio.
% Hint: get(J_Clust_obj,'Value') returns toggle state of disp_lratio

handles.lratio_select = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function display_comparative_unit_metrics_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to display_comparative_unit_metrics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in display_comparative_unit_metrics.

Comparative_Unit_Metrics;


function unit_timestamps_output_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to unit_timestamps_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in unit_timestamps_output.
% Hint: get(J_Clust_obj,'Value') returns toggle state of unit_timestamps_output

handles.output_ts = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function unit_waveforms_output_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to unit_waveforms_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in unit_waveforms_output.
% Hint: get(J_Clust_obj,'Value') returns toggle state of unit_waveforms_output

handles.output_waveforms = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function isi_time_output_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to isi_time_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in isi_time_output.
% Hint: get(J_Clust_obj,'Value') returns toggle state of isi_time_output


handles.output_isis = get(J_Clust_obj,'Value');

guidata(J_Clust_obj, handles);


function unit_spike_features_output_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to unit_spike_features_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in unit_spike_features_output.
% Hint: get(J_Clust_obj,'Value') returns toggle state of unit_spike_features_output

handles.output_features = get(J_Clust_obj,'Value');

features_output;

guidata(J_Clust_obj, handles);


function lratio_scores_output_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to lratio_scores_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in lratio_scores_output.
% Hint: get(J_Clust_obj,'Value') returns toggle state of lratio_scores_output

handles.output_lratio = get(J_Clust_obj, 'Value');

guidata(J_Clust_obj, handles);


function save_output_Callback(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to save_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% --- Executes on button press in save_output.

save_output;

guidata(J_Clust_obj, handles);


function J_Clust_fig_WindowKeyPressFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to J_Clust_fig (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on key press with focus on J_Clust_fig or any of its controls.

keyboard_plot_control;

guidata(J_Clust_obj, handles);


function J_Clust_fig_KeyPressFcn(J_Clust_obj, eventdata, handles)
% J_Clust_obj    handle to J_Clust_fig (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on key press with focus on J_Clust_fig or any of its controls.


% --- Executes during object deletion, before destroying properties.
function J_Clust_fig_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to J_Clust_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear global cur_coeff_pcs
clear global cur_coeff_pcs_c