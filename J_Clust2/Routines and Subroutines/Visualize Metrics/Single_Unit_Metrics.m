function varargout = Single_Unit_Metrics(varargin)

% Single_Unit_Metrics MATLAB code file for Single_Unit_Metrics.fig
%
%Description: This .m file opens a GUI containing plots of metrics used to compare the quality of a single, user-selected cluster (see 'Single Unit
%Metrics' in main GUI). The metrics can be used to remove spikes (theoretical false positives) from the selected cluster, directly via this GUI.
%
% Begin initialization code - DO NOT EDIT
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Single_Unit_Metrics_OpeningFcn, ...
                   'gui_OutputFcn',  @Single_Unit_Metrics_OutputFcn, ...
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


% --- Executes just before Single_Unit_Metrics is made visible.
function Single_Unit_Metrics_OpeningFcn(single_metrics_obj, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% single_metrics_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for Single_Unit_Metrics
handles.output = single_metrics_obj;

% Update handles structure
set(single_metrics_obj,'menubar','figure')
brush 'k';

J_Clust_tag = findobj('Tag','J_Clust_fig');
data_from_main = guidata(J_Clust_tag);

waveforms_select = data_from_main.waveforms_select;
residuals_select = data_from_main.residuals_select;
isi_hist_select = data_from_main.isi_hist_select;
autocorr_select = data_from_main.autocorr_select;
chebyshev_select = data_from_main.chebyshev_select;
handles.colors = data_from_main.colors;
polygon_indxs = data_from_main.polygon_indxs;
handles.Channel_Scatter = data_from_main.Channel_Scatter;
handles.Time_Scatter = data_from_main.Time_Scatter;
handles.feature_wires = data_from_main.feature_wires;
handles.feature_time = data_from_main.feature_time;
handles.chan_disp1 = data_from_main.chan_disp1;
handles.chan_disp2 = data_from_main.chan_disp2;
handles.features = data_from_main.features;
if data_from_main.preload
    handles.ts = data_from_main.ts(data_from_main.first_spk:data_from_main.last_spk);
else
    handles.ts = data_from_main.ts;
end
handles.remove_cheb_indxs = [];
handles.cl = find(data_from_main.view_1st_cl);

if handles.cl < 13
    handles.unit_pts = data_from_main.unit_pts{handles.cl};
    unit_title = ['Unit ', num2str(handles.cl)];
elseif handles.cl == 13
    handles.unit_pts = get_selected_pts();
    unit_title = ['Selected Spikes'];
else %handles.cl == 14 (all spikes)
    handles.unit_pts = data_from_main.unit_pts;
    unit_pts_all = [];
    for j = 1:length(handles.unit_pts)
        unit_pts_all = horzcat(unit_pts_all, handles.unit_pts{j});
    end
    all_pts = [1:length(handles.ts)];
    all_pts(unit_pts_all) = [];
    handles.unit_pts = all_pts;
    unit_title = ['All Spikes'];
end

num_samples = data_from_main.num_samples;
num_spks = length(handles.unit_pts);
waveforms = data_from_main.waveforms;
waveforms_unit = waveforms(:,:,handles.unit_pts);
concat_waveforms = reshape(permute(waveforms_unit,[2,1,3]), [4 * num_samples, num_spks]);
ts_unit = sort(handles.ts(handles.unit_pts));
filt_sig = data_from_main.filt_sig;
uV_conversion = data_from_main.uV_conversion;


if waveforms_select
    axes(handles.waveforms_ax)
    handles.waveforms_plot = plot(concat_waveforms); %'Color', colors(handles.cl,:));
    hold on
    ch2_divide = plot(num_samples,([min(concat_waveforms(:)):4:max(concat_waveforms(:))]), 'k.');
    ch3_divide = plot(num_samples*2,([min(concat_waveforms(:)):4:max(concat_waveforms(:))]), 'k.');
    ch4_divide = plot(num_samples*3,([min(concat_waveforms(:)):4:max(concat_waveforms(:))]), 'k.');
    hold off
    title([unit_title, ' Waveforms'], 'FontSize', 16)
else
    set(handles.waveforms_ax, 'Visible', 'off');
    set(handles.remove_waveforms_btn, 'Visible', 'off');
end

if residuals_select
    if data_from_main.preload
        disp('Cannot display cluster residuals without continuous tetrode signal')
    else
        mean_concat_waveform = mean(concat_waveforms,2);
        std_concat_waveforms = std(concat_waveforms,0,2);
        std_noise_channels = zeros(4,1);
        std_noise_channels_samples = [];
        parfor i = 1:4
            std_noise_channels(i) = std(filt_sig(i,:));
            std_noise_channels_samples = horzcat(std_noise_channels_samples, ones(num_samples,1) * std_noise_channels(i));
        end
        axes(handles.residuals_ax)
        background_noise_plot = plot(std_noise_channels_samples(:).*uV_conversion, 'Color', [0 .447 .741], 'LineWidth', 1.25);
        hold on
        cluster_residual_plot = plot(std_concat_waveforms, 'Color', handles.colors(handles.cl,:), 'LineWidth', 1.5);
        hold off
        legend('std of background noise', 'std of unit waveforms');
        xlabel('Waveform Sample Number');
        title([unit_title, ' Waveform Residuals'], 'FontSize', 15);
    end
else
    set(handles.residuals_ax, 'Visible', 'off');
end

if isi_hist_select
    handles.ISIs = diff(ts_unit) * 1000; %convert to ms
    short_ISIs = handles.ISIs(handles.ISIs < 10);
    
    axes(handles.isi_hist_ax)
    isi_hist = histogram(short_ISIs, 40, 'FaceColor', handles.colors(handles.cl,:));
    title([unit_title, ' ISI histogram'], 'FontSize', 14);
    ylabel('Count');
    xlabel('Time(ms)');
else
    set(handles.isi_hist_ax, 'Visible', 'off');
end

if autocorr_select
    all_autocorrs = cell(length(ts_unit),1);
    parfor i = 1:length(ts_unit)
        all_autocorrs{i} = ts_unit - ts_unit(i);
    end
    handles.all_autocorrs_vec = cell2mat(all_autocorrs) *1000;
    axes(handles.autocorr_ax)
    autocorr_hist = histogram(handles.all_autocorrs_vec, -50:2.5:50, 'FaceColor', handles.colors(handles.cl,:));
    set(handles.autocorr_ax, 'YTickLabel', []);
    xlabel('Time(ms)');
    ylabel('Count');
    title([unit_title, ' AutoCorrelogram of Spike Times'], 'FontSize', 14);
else
    set(handles.autocorr_ax, 'Visible', 'off');
end

if chebyshev_select
    %...calculate using Mahal distance of feature_wires distribution space
else
    set(handles.chebyshev_ax, 'Visible', 'off');
    set(handles.chebyshev_btn, 'Visible', 'off');
end

guidata(single_metrics_obj, handles);

% UIWAIT makes Single_Unit_Metrics wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Single_Unit_Metrics_OutputFcn(single_metrics_obj, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% single_metrics_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in chebyshev_btn.
function chebyshev_btn_Callback(single_metrics_obj, eventdata, handles)
% single_metrics_obj    handle to chebyshev_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

num_spks_unit = length(handles.unit_pts);
pc_scores = handles.features{4}; pc_1_scores_unit = squeeze(pc_scores(:,handles.unit_pts,1)); 
peak_amps = handles.features{1}; peak_amps_unit = peak_amps(:,handles.unit_pts);
power = handles.features{3}; power_unit = power(:,handles.unit_pts); 
unit_feature_set = [peak_amps_unit; power_unit; pc_1_scores_unit]'; 
rep_mean_unit_feature_set = repmat(mean(unit_feature_set),num_spks_unit,1);
cov_unit_feature_set = cov(unit_feature_set);
Mahal_distance_unit_matrix = (unit_feature_set - rep_mean_unit_feature_set) * pinv(cov_unit_feature_set) * (unit_feature_set - rep_mean_unit_feature_set)'; %this is a 'number_spikes' x 'number_spikes' matrix with M-values along diagonal
Mahal_distance_unit_vec = diag(Mahal_distance_unit_matrix);

mean_Mahal_unit = mean(Mahal_distance_unit_vec);
std_Mahal_unit = std(Mahal_distance_unit_vec);
k = str2num(cell2mat(inputdlg('Enter number of standard deviations away from mean:', 'Chebyshevs Inequality for Outliers', [1 50])));
max_far_spikes = ceil(num_spks_unit * 1/k^2);
far_spikes = find(Mahal_distance_unit_vec > (mean_Mahal_unit + k*std_Mahal_unit));
if length(far_spikes) < max_far_spikes
    disp('No Violations of Chebyshevs Inequality for this unit.')
    handles.remove_cheb_indxs = [];
else
    num_spks_to_remove = length(far_spikes) - max_far_spikes;
    [~, cheb_indxs_sort] = sort(Mahal_distance_unit_vec(far_spikes), 'descend');
    handles.remove_cheb_indxs = cheb_indxs_sort(1:num_spks_to_remove);
    axes(handles.chebyshev_ax)
    scatter(handles.feature_wires(handles.chan_disp1(1), unit_pts), handles.feature_wires(handles.chan_disp1(2), unit_pts), '.')
    hold on
    scatter(handles.feature_wires(handles.chan_disp1(1), handles.remove_cheb_indxs), handles.feature_wires(handles.chan_disp1(2), handles.remove_cheb_indxs), '.')
    hold off
end

guidata(single_metrics_obj, handles);

% --- Executes on button press in remove_waveforms_btn.
function remove_waveforms_btn_Callback(single_metrics_obj, eventdata, handles)
% single_metrics_obj    handle to remove_waveforms_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

J_Clust_tag = findobj('Tag','J_Clust_fig');
data_from_main = guidata(J_Clust_tag);

for i = 1:length(handles.waveforms_plot)
    brushed_waveform_pts(i) = sum(find(handles.waveforms_plot(i).BrushData));
end

handles.remove_selected_indxs = find(brushed_waveform_pts);
delete(handles.waveforms_plot(handles.remove_selected_indxs));
handles.waveforms_plot(handles.remove_selected_indxs) = [];

handles.remove_selected_indxs = horzcat(handles.remove_selected_indxs, handles.remove_cheb_indxs); 

axes(handles.Channel_Scatter)
hold on
if length(handles.chan_disp1) > 2
    metric_scat_wires = scatter3(handles.feature_wires(handles.chan_disp1(1),handles.unit_pts(handles.remove_selected_indxs)), handles.feature_wires(handles.chan_disp1(2),handles.unit_pts(handles.remove_selected_indxs)), handles.feature_wires(handles.chan_disp1(3),handles.unit_pts(handles.remove_selected_indxs)), 90, 'k.');
else %length(chan_disp) == 2
    metric_scat_wires = scatter(handles.feature_wires(handles.chan_disp1(1),handles.unit_pts(handles.remove_selected_indxs)), handles.feature_wires(handles.chan_disp1(2),handles.unit_pts(handles.remove_selected_indxs)), 90, 'k.');
end
hold off

no_time_disp = 0;
if strcmp(data_from_main.feat_disp2, 'X vs. Y') || strcmp(data_from_main.feat_disp2, 'View All Waveforms')
    no_time_disp = 1;
end
if ~no_time_disp
    axes(handles.Time_Scatter)
    hold on
    metrics_scat_time = scatter(handles.ts(handles.unit_pts(handles.remove_selected_indxs)), handles.feature_time(handles.chan_disp2,handles.unit_pts(handles.remove_selected_indxs)), 90, 'k.');
    hold off
end

guidata(single_metrics_obj, handles);


function isi_hist_ax_ButtonDownFcn(single_metrics_obj, eventdata, handles)
% single_metrics_obj    handle to isi_hist_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on mouse press over axes background.

%if right click, user sets xlim (milliseconds) and number of bins

if strcmp(get(handles.single_metrics_fig, 'SelectionType'), 'alt')
    isi_hist_prompts = {'Enter minimum ISI time (ms):', 'Enter maximum ISI time (ms):', 'Enter number of bins for histogram:'};
    isi_hist_defaults = {'0', '10', '40'};
    isi_hist_title = 'ISI Histogram Options';
    isi_hist_user_input = inputdlg(isi_hist_prompts, isi_hist_title, [1 50], isi_hist_defaults);
    if isempty(isi_hist_user_input)
        return;
    end
    
    user_ISIs = handles.ISIs(handles.ISIs > str2num(isi_hist_user_input{1}) & handles.ISIs < str2num(isi_hist_user_input{2}));
    num_bins = str2num(isi_hist_user_input{3});
    
    axes(handles.isi_hist_ax)
    isi_hist = histogram(user_ISIs, num_bins, 'FaceColor', handles.colors(handles.cl,:));
end

%if double-left click, user can choose spikes to remove

if strcmp(get(handles.single_metrics_fig, 'SelectionType'), 'open')
    
    J_Clust_tag = findobj('Tag','J_Clust_fig');
    data_from_main = guidata(J_Clust_tag);
    
    remove_isi_prompt = 'Enter max ISI time (in ms). All spikes with ISIs less than this value will be removed.';
    remove_isi_default = {'1.5'};
    remove_isi_title = 'Remove Spikes Based on ISIs';
    remove_isi_val = inputdlg(remove_isi_prompt, remove_isi_title, [1 50], remove_isi_default);
    
    if isempty(remove_isi_val)
        return;
    end
    
    remove_isi_val = str2num(remove_isi_val{1});
    
    [ts_unit, srtd_indxs] = sort(handles.ts(handles.unit_pts));
    remove_unit_indxs = find(handles.ISIs < remove_isi_val);
    if ~isempty(remove_unit_indxs)
        [~, handles.remove_selected_indxs, ~] = intersect(srtd_indxs, remove_unit_indxs);
        
        axes(handles.Channel_Scatter)
        hold on
        if length(handles.chan_disp1) > 2
            metric_scat_wires = scatter3(handles.feature_wires(handles.chan_disp1(1),handles.unit_pts(handles.remove_selected_indxs)), handles.feature_wires(handles.chan_disp1(2),handles.unit_pts(handles.remove_selected_indxs)), handles.feature_wires(handles.chan_disp1(3),handles.unit_pts(handles.remove_selected_indxs)), 90, 'k.');
        else %length(chan_disp) == 2
            metric_scat_wires = scatter(handles.feature_wires(handles.chan_disp1(1),handles.unit_pts(handles.remove_selected_indxs)), handles.feature_wires(handles.chan_disp1(2),handles.unit_pts(handles.remove_selected_indxs)), 90, 'k.');
        end
        hold off
        
        no_time_disp = 0;
        if strcmp(data_from_main.feat_disp2, 'X vs. Y') || strcmp(data_from_main.feat_disp2, 'View All Waveforms')
            no_time_disp = 1;
        end
        if ~no_time_disp
            axes(handles.Time_Scatter)
            hold on
            metrics_scat_time = scatter(handles.ts(handles.unit_pts(handles.remove_selected_indxs)), handles.feature_time(handles.chan_disp2,handles.unit_pts(handles.remove_selected_indxs)), 90, 'k.');
            hold off
        end
    end
end

guidata(single_metrics_obj, handles);

    
function autocorr_ax_ButtonDownFcn(single_metrics_obj, eventdata, handles)
% single_metrics_obj    handle to autocorr_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on mouse press over axes background.

%if right click, user sets xlim (milliseconds) and number of bins

if strcmp(get(handles.single_metrics_fig, 'SelectionType'), 'alt')
    autocorr_prompts = {'Enter minimum time (ms):', 'Enter maximum time (ms):', 'Enter number of bins for histogram:'};
    autocorr_defaults = {'-50', '50', '40'};
    autocorr_title = 'AutoCorrelogram Options';
    autocorr_user_input = inputdlg(autocorr_prompts, autocorr_title, [1 50], autocorr_defaults);
    if isempty(autocorr_user_input)
        return;
    end
    
    user_autocorr = handles.all_autocorrs_vec(handles.all_autocorrs_vec > str2num(autocorr_user_input{1}) & handles.all_autocorrs_vec < str2num(autocorr_user_input{2}));
    num_bins = str2num(autocorr_user_input{3});
    
    axes(handles.autocorr_ax)
    autocorr_hist = histogram(user_autocorr, num_bins, 'FaceColor', handles.colors(handles.cl,:));
    set(handles.autocorr_ax, 'YTickLabel', []);
end
