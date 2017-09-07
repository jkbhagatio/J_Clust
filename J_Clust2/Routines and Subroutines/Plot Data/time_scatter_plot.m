function feature = time_scatter_plot(chan_disp, feat_disp, features, ts, varargin)
 
%Description: This .m file is called when plotting spike features (as a scatter plot) on the 'Time Scatter' axis in the main GUI.
%
%Input: 1) single channel of tetrode to display, 2) spike feature to display, 3) cell array of all spike features cell array, 4) timestamps of all 
%spikes, 5) varargin{1} = points of any units already clustered, varargin{2} = cluster colors
%

if nargin > 4
    unit_pts = varargin{1};
    colors = varargin{2};
end

xlabname = ['Channel ', num2str(chan_disp), ', Time (s)'];
switch feat_disp
    case 'Peak Amplitude'
        feature = features{1,1};
    case 'Crest-to-Trough Amplitude'
        feature = features{2,1};
    case 'Power'
        feature = features{3,1};
    case 'PCA Scores'
        feature = features{4,1};
        global cur_coeff_pcs
        if isempty(cur_coeff_pcs)
            Select_Coefficient_PCs %add conditionals for if initial scatter or scatter after units are defined
            uiwait(Select_Coefficient_PCs)
            
            feature = squeeze(feature(:,:, cur_coeff_pcs));
        else
            feature = squeeze(feature(:,:,cur_coeff_pcs));
        end
  
        feat_disp = ['PC ', num2str(cur_coeff_pcs), ' Scores'];
        
    case 'PCA Scores for Concatenated Waveforms'
        feature = features{5,1}';
        feat_disp = ['Concatenated PC', num2str(chan_disp), ' Scores'];
        xlabname = 'Time(s)';
        
    %case 'Wavelet Coefficients'
        %...

    %case 'Wavelet Coefficients for Concatenated Waveforms'
        %...
    case 'X vs. Y'
        feature = [];
        [filename, filepath] = uigetfile;
        addpath(genpath(filepath));
        user_data = load(filename);
        data_fieldnames = fieldnames(user_data);
        pos_mtx = user_data.(data_fieldnames{1});
        pos_ts = pos_mtx(:,3);
        pos_x = pos_mtx(:,1);
        pos_y = pos_mtx(:,2);
        
        xpos = zeros(length(ts),1);
        ypos = zeros(length(ts), 1);
        parfor i = 1:length(ts)
            [~, pos_indx] = min(abs(pos_ts - ts(i)));
            xpos(i) = pos_x(pos_indx);
            ypos(i) = pos_y(pos_indx);
        end
        
        h2_scat = scatter(xpos, ypos, '.');
        if nargin > 5 %some units have already been defined
            hold on
            for i = 1:length(unit_pts)     
                h2_scat_units{i} = scatter(xpos(unit_pts{i}), ypos(unit_pts{i}), 36, colors(i,:), '.');
            end
            hold off
        end
        
        title('X vs. Y Position for Each Spike Time')
      
        return;
    case 'View All Waveforms'
        feature = [];
        all_unit_pts = [];
        J_Clust_tag = findobj('Tag','J_Clust_fig');
        data_from_main = guidata(J_Clust_tag);
        waveforms = data_from_main.waveforms;
        waveforms = squeeze(waveforms(chan_disp,:,:));
        if nargin > 5 %some units have already been defined
            for i = 1:length(unit_pts)
                plot(waveforms(:,unit_pts{i}), 'Color', colors(i,:))
                hold on
                all_unit_pts = [all_unit_pts, unit_pts{i}];
            end
        end
        all_pts = [1:length(ts)];
        all_pts(all_unit_pts) = [];
        plot(waveforms(:,all_pts), 'Color', [0 .447 .741])
        hold off
        return;

end
%% Plot Data

scatter(ts, feature(chan_disp,:),'.');
   
if nargin > 5 %some units have already been clustered
    %add new unit scatters based on number of clusters
    hold on
    for i = 1:length(unit_pts)
        scatter(ts(unit_pts{i}), feature(chan_disp,unit_pts{i}), 36, colors(i,:), '.');
    end
    hold off
end

xlabel(xlabname)
title(feat_disp)

hold off