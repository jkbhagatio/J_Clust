function feature = channel_scatter_plot(chan_disp, feat_disp, features, varargin)

%Description: This .m file is called when plotting spike features (as a scatter plot) on the 'Channel Scatter' axis in the main GUI.
%
%Input: 1) channels of tetrode to display, 2) spike feature to display, 3) cell array of all spike features cell array,
% 4) varargin{1} = points of any units already clustered, varargin{2} = cluster colors
%

if nargin > 3
    unit_pts = varargin{1};
    colors = varargin{2};
end

%Select feature based on user input

%default for axes labels
xlabname = ['Wire ', num2str(chan_disp(1))];
ylabname = ['Wire ', num2str(chan_disp(2))];
if length(chan_disp) > 2
    zlabname = ['Wire ', num2str(chan_disp(3))];
end

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
        feature = features{5,1};
        global cur_coeffs_pcs_c
        if isempty(cur_coeffs_pcs_c)
            Select_Coefficient_PCs_c
            uiwait(Select_Coefficient_PCs_c)
        
            feature = feature(:, cur_coeffs_pcs_c)';
            if isempty(feature)
                feature = feature'; %for plotting purposes, to avoid an error
            end
        else
            feature = feature(:, cur_coeffs_pcs_c)';
        end
        feat_disp = ['Concatenated PC Scores'];
        xlabname = ['PC ', num2str(chan_disp(1))];
        ylabname = ['PC ', num2str(chan_disp(2))];
        if length(chan_disp) > 2
            zlabname = ['PC ', num2str(chan_disp(3))];
        end
        
    %case 'Wavelet Coefficients'
        %...

            
    %case 'Wavelet Coefficients for Concatenated Waveforms'
        %...
end
%% Plot Data

if length(chan_disp) > 2
    scatter3(feature(chan_disp(1),:), feature(chan_disp(2),:), feature(chan_disp(3),:), '.');
    zlabel(zlabname)
else %length(chan_disp) == 2
    scatter(feature(chan_disp(1),:), feature(chan_disp(2),:),'.');
end


if nargin > 4 %some units have already been clustered
    if length(chan_disp) > 2
        hold on
        for i = 1:length(unit_pts)
            scatter3(feature(chan_disp(1), unit_pts{i}), feature(chan_disp(2), unit_pts{i}), feature(chan_disp(3), unit_pts{i}), 36, colors(i,:), '.');
        end
        hold off
        zlabel(zlabname)
    else %length(chan_disp) == 2
        hold on
        for i = 1:length(unit_pts)
            scatter(feature(chan_disp(1), unit_pts{i}), feature(chan_disp(2), unit_pts{i}), 36, colors(i,:), '.');
        end
        hold off
    end
end

xlabel(xlabname)
ylabel(ylabname)
title(feat_disp)

hold off
        