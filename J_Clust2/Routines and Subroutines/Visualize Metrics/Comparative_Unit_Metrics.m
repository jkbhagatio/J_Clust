function varargout = Comparative_Unit_Metrics(varargin)
% Comparative_Unit_Metrics MATLAB code file for Comparative_Unit_Metrics.fig
%
%Description: This .m file opens a GUI containing plots of metrics used to compare two user-selected clusters together, jointly 
%(see 'Comparative Metrics' button  group in main GUI). The metrics can be used to reassign spikes to one cluster or another, or 
%remove spikes from one cluster, directly via this GUI.
%
% Begin initialization code - DO NOT EDIT
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Comparative_Unit_Metrics_OpeningFcn, ...
                   'gui_OutputFcn',  @Comparative_Unit_Metrics_OutputFcn, ...
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


% --- Executes just before comparative_metrics_fig is made visible.
function Comparative_Unit_Metrics_OpeningFcn(comparative_metrics_obj, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% comparative_metrics_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to comparative_metrics_fig (see VARARGIN)

% Choose default command line output for comparative_metrics_fig
handles.output = comparative_metrics_obj;

% Update handles structure
set(comparative_metrics_obj,'menubar','figure')

J_Clust_tag = findobj('Tag','J_Clust_fig');
data_from_main = guidata(J_Clust_tag);
if data_from_main.preload
    handles.ts = data_from_main.ts(data_from_main.first_spk:data_from_main.last_spk);
else
    handles.ts = data_from_main.ts;
end
handles.Channel_Scatter = data_from_main.Channel_Scatter;
handles.Time_Scatter = data_from_main.Time_Scatter;
handles.chan_disp1 = data_from_main.chan_disp1;
handles.feat_disp1 = data_from_main.feat_disp1;
handles.feature_wires = data_from_main.feature_wires;
handles.chan_disp2 = data_from_main.chan_disp2;
handles.feat_disp2 = data_from_main.feat_disp2;
handles.feature_time = data_from_main.feature_time;
handles.remove_selected_indxs = [];

silhouette_select = data_from_main.silhouette_select;
crosscorr_select = data_from_main.crosscorr_select;
m_distance_select = data_from_main.m_distance_select;
lratio_select = data_from_main.lratio_select;
handles.colors = data_from_main.colors;

handles.cl1 = find(data_from_main.view_1st_cl);
handles.cl2 = find(data_from_main.view_2nd_cl);

if handles.cl1 < 13
    handles.unit_pts1 = data_from_main.unit_pts{handles.cl1};
    unit_title1 = ['Unit ', num2str(handles.cl1)];
elseif handles.cl1 == 13
    handles.unit_pts1 = data_from_main.selected_pts;
    unit_title1 = ['Selected Spikes'];
else %handles.cl1 == 14 (all spikes)
    handles.unit_pts = data_from_main.unit_pts;
    unit_pts_all = [];
    for j = 1:length(handles.unit_pts)
        unit_pts_all = horzcat(unit_pts_all, handles.unit_pts{j});
    end
    all_pts = [1:length(handles.ts)];
    all_pts(unit_pts_all) = [];
    handles.unit_pts1 = all_pts;
    unit_title1 = ['All Spikes'];
end

if handles.cl2 < 13
    handles.unit_pts2 = data_from_main.unit_pts{handles.cl2};
    unit_title2 = ['Unit ', num2str(handles.cl2)];
elseif handles.cl2 == 13
    handles.unit_pts2 = get_selected_pts;
    unit_title2 = ['Selected Spikes'];
else %handles.cl2 == 14 (all spikes)
    handles.unit_pts = data_from_main.unit_pts;
    unit_pts_all = []';
    for j = 1:length(handles.unit_pts)
        unit_pts_all = horzcat(unit_pts_all, handles.unit_pts{j});
    end
    all_pts = [1:length(handles.ts)];
    all_pts(unit_pts_all) = [];
    handles.unit_pts2 = all_pts;
    unit_title2 = ['All Spikes'];
end

num_spks_unit1 = length(handles.unit_pts1);
num_spks_unit2 = length(handles.unit_pts2);
ts_unit1 = sort(handles.ts(handles.unit_pts1));
ts_unit2 = sort(handles.ts(handles.unit_pts2));

features = data_from_main.features;
pc_scores = features{4};
ts = handles.ts;
ts(handles.unit_pts1) = 0;
noise_pts = find(ts);
pc_1_scores_unit1 = squeeze(pc_scores(:,handles.unit_pts1,1)); pc_1_scores_unit2 = squeeze(pc_scores(:,handles.unit_pts2,1)); pc1_scores_noise = squeeze(pc_scores(:,noise_pts,1));
peak_amps = features{1};
peak_amps_unit1 = peak_amps(:,handles.unit_pts1); peak_amps_unit2 = peak_amps(:,handles.unit_pts2); peak_amps_noise = peak_amps(:,noise_pts); 
power = features{3};
power_unit1 = power(:,handles.unit_pts1); power_unit2 = power(:,handles.unit_pts2); power_noise = power(:,noise_pts);
unit_feature_set1 = [peak_amps_unit1; power_unit1; pc_1_scores_unit1]'; 
unit_feature_set2 = [peak_amps_unit2; power_unit2; pc_1_scores_unit2]';
noise_feature_set = [peak_amps_noise; power_noise; pc1_scores_noise]';
rep_mean_unit_feature_set1 = repmat(mean(unit_feature_set1),num_spks_unit1,1); rep_mean_unit_feature_set2 = repmat(mean(unit_feature_set1),num_spks_unit2,1); rep_mean_unit_feature_set1_noise = repmat(mean(unit_feature_set1), length(noise_pts), 1);
cov_unit_feature_set1 = cov(unit_feature_set1); cov_unit_feature_set2 = cov(unit_feature_set2); cov_noise_feature_set = cov(noise_feature_set);
Mahal_distance_unit_matrix1 = (unit_feature_set1 - rep_mean_unit_feature_set1) * pinv(cov_unit_feature_set1) * (unit_feature_set1 - rep_mean_unit_feature_set1)'; %this is a 'number_spikes' x 'number_spikes' matrix with M-values along diagonal
handles.Mahal_distance_unit_vec1 = diag(Mahal_distance_unit_matrix1);
[Mahal_distance_unit_vec_sort1, Mahal_unit_id1] = sort(handles.Mahal_distance_unit_vec1);
Mahal_distance_unit_matrix2 = (unit_feature_set2 - rep_mean_unit_feature_set2) * pinv(cov_unit_feature_set2) * (unit_feature_set2 - rep_mean_unit_feature_set2)'; %this is a 'number_spikes' x 'number_spikes' matrix with M-values along diagonal
handles.Mahal_distance_unit_vec2 = diag(Mahal_distance_unit_matrix2);
[Mahal_distance_unit_vec_sort2, Mahal_unit_id2] = sort(handles.Mahal_distance_unit_vec2);
Mahal_distance_noise_matrix = (noise_feature_set - rep_mean_unit_feature_set1_noise) * pinv(cov_unit_feature_set1) * (noise_feature_set - rep_mean_unit_feature_set1_noise)'; %this is a 'number_spikes' x 'number_spikes' matrix with M-values along diagonal
Mahal_distance_noise_vec = diag(Mahal_distance_noise_matrix);
[Mahal_distance_noise_vec_sort, Mahal_noise_id] = sort(Mahal_distance_noise_vec);

if m_distance_select
    unit_Mahal_bin_cntrs1 = [0, 0:5:num_spks_unit1]; %add 0 at beginning so first point is at 0,0
    unit_Mahal_bin_counts1 = hist(handles.Mahal_distance_unit_vec1, unit_Mahal_bin_cntrs1);
    unit_Mahal_bin_cntrs2 = [0, 0:5:num_spks_unit2]; %add 0 at  beginning so first point is at 0,0
    unit_Mahal_bin_counts2 = hist(handles.Mahal_distance_unit_vec2, unit_Mahal_bin_cntrs2);
    
    axes(handles.mahal_ax)
    handles.mahal_plot2 = plot(unit_Mahal_bin_cntrs2, unit_Mahal_bin_counts2, 'Color', handles.colors(handles.cl2,:),'LineWidth', 1.2);
    hold on, handles.mahal_plot = plot(unit_Mahal_bin_cntrs1, unit_Mahal_bin_counts1, 'Color', handles.colors(handles.cl1,:),'LineWidth', 1.2); hold off
    xlim([0 max(handles.Mahal_distance_unit_vec1)])
    legend(unit_title2, unit_title1)
    title(['Mahalanobis Distance of Spikes to ', unit_title1, ' Cluster Center'], 'FontSize', 9);
    ylabel('Count', 'FontSize', 9)
    xlabel('Mahalanobis Distance', 'FontSize', 9)
else
    set(handles.mahal_ax, 'Visible', 'off');
end

if crosscorr_select    
    all_crosscorrs = cell(length(ts_unit1),1);
    parfor i = 1:length(ts_unit1)
        all_crosscorrs{i} = ts_unit2 - ts_unit1(i);
    end
    handles.all_crosscorrs_vec = cell2mat(all_crosscorrs) *1000;
    
    axes(handles.crosscorr_ax)
    handles.crosscorr_hist = histogram(handles.all_crosscorrs_vec, -50:2.5:50, 'FaceColor', handles.colors(handles.cl1,:));
    set(handles.crosscorr_ax, 'YTickLabel', []);
    xlabel('Time(ms)', 'FontSize', 9);
    ylabel('Count', 'FontSize', 9);
    title([' CrossCorrelogram of ', unit_title1,' Spike Times with ', unit_title2, ' Spike Times'], 'FontSize', 9);
else
    set(handles.crosscorr_ax, 'Visible', 'off');
end

if silhouette_select
    sil = zeros(length(handles.unit_pts1),1);
    parfor i=1:length(ts_unit1)
        avg_unit1_dis = abs(mean(handles.Mahal_distance_unit_vec1 - handles.Mahal_distance_unit_vec1(i)));
        avg_unit2_dis = abs(mean(handles.Mahal_distance_unit_vec2 - handles.Mahal_distance_unit_vec1(i)));
        sil(i) = (avg_unit2_dis - avg_unit1_dis) / max(avg_unit1_dis, avg_unit2_dis);
    end
    
    axes(handles.silhouette_ax)
    handles.sil_bar = bar(sil);
    title(['Silhouette Scores of Each Spike in ', unit_title1, ' compared to ', unit_title2], 'FontSize', 9);
    xlabel('Spike Number', 'FontSize', 9);
    ylabel('Score', 'FontSize', 9);
else
    set(handles.silhouette_ax, 'Visible', 'off');
end

if lratio_select
    feat_to_test = handles.feature_wires([handles.chan_disp1(1) handles.chan_disp1(2)], handles.unit_pts1);
    [h, p] = kstest2(feat_to_test(1,:), feat_to_test(2,:));
    if h
        warning('This units spikes do not have a multivariate normal distribution under the current feature, so the L-Ratio will not be calculated')
        set(handles.lratio_ax, 'Visible', 'off');
    else
        L_ones_noise_vec = ones(size(noise_pts));
        L_noise_vec = L_ones_noise_vec - chi2cdf(Mahal_distance_noise_vec, 12); %12 degrees of freedom for 12-feature matrix used to calculate M-distance
        L = sum(L_noise_vec(:));
        L_ratio = L / num_spks_unit1;
        
        axes(handles.lratio_ax)
        plot([0; Mahal_distance_noise_vec_sort], [1; sort(L_noise_vec, 'descend')], 'Color', handles.colors(handles.cl1,:), 'LineWidth', 1.3)
        title(['L-ratio Plot for ', unit_title1, ', {1-cdf(x^2)}'], 'FontSize', 9);
        xlabel('Mahalanobis Distance', 'FontSize', 9);
        ylabel('|P| in unit', 'FontSize', 9);     
    end
else
    set(handles.lratio_ax, 'Visible', 'off');
end
brush 'k';

guidata(comparative_metrics_obj, handles);

% UIWAIT makes comparative_metrics_fig wait for user response (see UIRESUME)
% uiwait(handles.comparative_metrics_fig);


% --- Outputs from this function are returned to the command line.
function varargout = Comparative_Unit_Metrics_OutputFcn(comparative_metrics_obj, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% comparative_metrics_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in view_spikes.
function view_spikes_Callback(comparative_metrics_obj, eventdata, handles)
% hObject    handle to view_spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

J_Clust_tag = findobj('Tag','J_Clust_fig');
data_from_main = guidata(J_Clust_tag);

mahal_line_plot_h = handles.mahal_plot.BrushHandles;
if ~isempty(mahal_line_plot_h)
    mahal_selected_xvals = mahal_line_plot_h.Children(2).VertexData(1,:);
    mahal_brushed_indxs = find(handles.Mahal_distance_unit_vec1 > (mahal_selected_xvals(1) - 2.5) & (mahal_selected_xvals(end) + 2.5));
else
    mahal_brushed_indxs = [];
end
%mahal_selected_yvals = mahal_line_plot_vals(2).VertexData(2,:);

silhouette_brushed_indxs = find(get(handles.sil_bar, 'BrushData'));

handles.remove_selected_indxs = unique(horzcat(mahal_brushed_indxs, silhouette_brushed_indxs));

axes(handles.Channel_Scatter)
hold on
if length(handles.chan_disp1) > 2
    metrics_scat_wires = scatter3(handles.feature_wires(handles.chan_disp1(1),handles.unit_pts1(handles.remove_selected_indxs)), handles.feature_wires(handles.chan_disp1(2),handles.unit_pts1(handles.remove_selected_indxs)), handles.feature_wires(handles.chan_disp1(3),handles.unit_pts1(handles.remove_selected_indxs)), 90, 'k.');
else %length(chan_disp) == 2
    metrics_scat_wires = scatter(handles.feature_wires(handles.chan_disp1(1),handles.unit_pts1(handles.remove_selected_indxs)), handles.feature_wires(handles.chan_disp1(2),handles.unit_pts1(handles.remove_selected_indxs)), 90, 'k.');
end
hold off

no_time_disp = 0;
if strcmp(data_from_main.feat_disp2, 'X vs. Y') || strcmp(data_from_main.feat_disp2, 'View All Waveforms')
    no_time_disp = 1;
end
if ~no_time_disp
    axes(handles.Time_Scatter)
    hold on
    metrics_scat_time = scatter(handles.ts(handles.unit_pts1(handles.remove_selected_indxs)), handles.feature_time(handles.chan_disp2,handles.unit_pts1(handles.remove_selected_indxs)), 90, 'k.');
    hold off
end

guidata(comparative_metrics_obj, handles);


function crosscorr_ax_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to crosscorr_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on mouse press over axes background.

if strcmp(get(handles.comparative_metrics_fig, 'SelectionType'), 'alt')
    crosscorr_prompts = {'Enter minimum time (ms):', 'Enter maximum time (ms):', 'Enter number of bins for histogram:'};
    crosscorr_defaults = {'-50', '50', '40'};
    crosscorr_title = 'CrossCorrelogram Options';
    crosscorr_user_input = inputdlg(crosscorr_prompts, crosscorr_title, [1 50], crosscorr_defaults);
    if isempty(crosscorr_user_input)
        return;
    end
    
    user_crosscorr = handles.all_crosscorrs_vec(handles.all_crosscorrs_vec > str2num(crosscorr_user_input{1}) & handles.all_crosscorrs_vec < str2num(crosscorr_user_input{2}));
    num_bins = str2num(crosscorr_user_input{3});
    
    axes(handles.crosscorr_ax)
    crosscorr_hist = histogram(user_crosscorr, num_bins, 'FaceColor', handles.colors(handles.cl1,:));
end
