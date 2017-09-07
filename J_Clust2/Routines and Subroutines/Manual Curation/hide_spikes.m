function [hide_spks_scat_wires, hide_spks_scat_time] = hide_spikes(unit_pts, edit_cl, chan_disp1, feature_wires, chan_disp2, feature_time, hide_spks_scat_wires, hide_spks_scat_time, ts, Channel_Scatter, Time_Scatter, feat_disp1, feat_disp2, colors)

spks_to_hide = find(edit_cl);

xlabname = ['Wire ', num2str(chan_disp1(1))];
ylabname = ['Wire ', num2str(chan_disp1(2))];
if length(chan_disp1) > 2
    zlabname = ['Wire ', num2str(chan_disp1(3))];
end

if strcmp(feat_disp1, 'PCA Scores for Concatenated Waveforms')
    xlabname = ['PC ', num2str(chan_disp1(1))];
    ylabname = ['PC ', num2str(chan_disp1(2))];
    if length(chan_disp1) > 2
        zlabname = ['PC ', num2str(chan_disp1(3))];
    end
end


for i = 1:length(spks_to_hide)
    if spks_to_hide(i) > 13
        all_pts = [1:length(feature_wires)];
        unit_pts_all = [];
        for j = 1:length(unit_pts)
            unit_pts_all = horzcat(unit_pts_all, unit_pts{j});
        end
        cur_spks_to_hide = all_pts;
        cur_spks_to_hide(unit_pts_all) = []; %all points that have not been assigned to a unit
    elseif spks_to_hide(i) > 12
        cur_spks_to_hide = get_selected_pts();
    else
        cur_spks_to_hide = unit_pts{spks_to_hide(i)};
    end
    
    axes(Channel_Scatter)
    if length(chan_disp1) > 2
        h1_scat = scatter3(feature_wires(chan_disp1(1),:), feature_wires(chan_disp1(2),:), feature_wires(chan_disp1(3),:), '.');
        zlabel(zlabname)
    else %length(chan_disp) == 2
        h1_scat = scatter(feature_wires(chan_disp1(1),:), feature_wires(chan_disp1(2),:),'.');
    end
    hold on
    
    for j = 1:length(unit_pts)
        if j == spks_to_hide(i)
            if length(chan_disp1) > 2
                hide_spks_scat_wires{i} = scatter3(feature_wires(chan_disp1(1), cur_spks_to_hide), feature_wires(chan_disp1(2), cur_spks_to_hide), feature_wires(chan_disp1(3), cur_spks_to_hide), 1e-6, 'w.');
            else %length(chan_disp) == 3
                hide_spks_scat_wires{i} = scatter(feature_wires(chan_disp1(1), cur_spks_to_hide), feature_wires(chan_disp1(2), cur_spks_to_hide), 1e-6, 'w.');
            end
            continue;
        end
        if length(chan_disp1) > 2
            h1_scat_units{j} = scatter3(feature_wires(chan_disp1(1), unit_pts{j}), feature_wires(chan_disp1(2), unit_pts{j}), feature_wires(chan_disp1(3), unit_pts{j}), 36, colors(j,:), '.');
        else %length(chan_disp) == 2
            h1_scat_units{j} = scatter(feature_wires(chan_disp1(1), unit_pts{j}), feature_wires(chan_disp1(2), unit_pts{j}), 36, colors(j,:), '.');
        end
    end
    hold off
    
    xlabel(xlabname)
    ylabel(ylabname)
    title(feat_disp1)
    
    
    axes(Time_Scatter)
    h2_scat = scatter(ts, feature_time(chan_disp,:),'.');
    
    hold on
    for j = 1:length(unit_pts)
        if j == spks_to_hide(i)
            hide_spks_scat_time{i} = scatter(ts(cur_spks_to_hide), feature_time(chan_disp2, cur_spks_to_hide), 1e-6, 'w.');
            continue;
        end
        h2_scat_units{j} = scatter(ts(unit_pts{j}), feature(chan_disp2,unit_pts{j}), 36, colors(j,:), '.');
    end
    hold off
    
    xlabel(['Channel ', num2str(chan_disp2), ', Time (s)'])
    if strcmp(feat_disp2, 'PCA Scores')
        title(['PC ' num2str(coeff)])
    elseif strcmp(feat_disp2, 'PCA Scores for Concatenated Waveforms')
        title(['Concatenated PC', num2str(chan_disp), ' Scores']);
    else
        title(feat_disp2)
    end
    
end
