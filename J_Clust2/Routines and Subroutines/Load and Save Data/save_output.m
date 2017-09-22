%Description: This .m file is a script which reads in the user selections from the 'Save Options' button group in the main GUI, and saves the output
%of those selections to the MATLAB Workspace

if handles.output_ts
    for i = 1:length(handles.unit_pts)
        if isempty(handles.unit_pts{i})
            continue;
        end
        unit_ts = handles.ts(handles.unit_pts{i});
        field = ['unit', num2str(i)];
        output_struct.(field).ts = unit_ts;
    end
end

if handles.output_waveforms
    for i = 1:length(handles.unit_pts)
        if isempty(handles.unit_pts{i})
            continue;
        end
        unit_waveforms = handles.waveforms(:,:,handles.unit_pts{i});
        field = ['unit', num2str(i)];
        output_struct.(field).waveforms = unit_waveforms;
    end
end


if handles.output_isis
    for i = 1:length(handles.unit_pts)
        if isempty(handles.unit_pts{i})
            continue;
        end
        unit_ts = sort(handles.ts(handles.unit_pts{i}));
        unit_isis = diff(unit_ts) * 1000; %converted to ms
        field = ['unit', num2str(i)];
        output_struct.(field).ISIs = unit_isis;
    end
end

if handles.output_features
    global features_vector_out
    for i = 1:length(handles.unit_pts)
        if isempty(handles.unit_pts{i})
            continue;
        end
        for j = 1:length(features_vector_out)
            if features_vector_out(j)
                cur_feature = handles.features{j,1};
                cur_feature_name = handles.features{j,2};
                if strcmp(cur_feature_name, 'PC_Scores')
                    unit_feature = cur_feature(:,handles.unit_pts{i},:);
                elseif strcmp(cur_feature_name, 'Concatenated_PC_Scores')
                    unit_feature = cur_feature(handles.unit_pts{i},:);
                else
                    unit_feature = cur_feature(:,handles.unit_pts{i});
                end
                field = ['unit', num2str(i)];
                output_struct.(field).(cur_feature_name) = unit_feature;
            end
        end
    end
    clear global features_vector_out
end

if handles.output_lratio
    for i = 1:length(handles.unit_pts)
        if isempty(handles.unit_pts{i})
            continue;
        end
        field = ['unit', num2str(i)];
        feat_to_test = handles.feature_wires([handles.chan_disp1(1) handles.chan_disp1(2)], handles.unit_pts{i});
        [h, p] = kstest2(feat_to_test(1,:), feat_to_test(2,:));
        if h %if distribution of spikes in current feature space is non-multivariate normal
            output_struct.(field).L_ratio = 'undefined';
            continue;
        else %calculate L-ratio
            if handles.preload
                ts = handles.ts(handles.first_spk:handles.last_spk);
            else
                ts = handles.ts;
            end
            ts(handles.unit_pts{i}) = 0;
            noise_pts = find(ts);
            peak_amps = handles.features{1};
            power = handles.features{3};
            pc_scores = handles.features{4};
            peak_amps_noise = peak_amps(:,noise_pts);
            power_noise = power(:,noise_pts);
            pc1_scores_noise = squeeze(pc_scores(:,noise_pts,1));
            noise_feature_set = [peak_amps_noise; power_noise; pc1_scores_noise]';
            rep_mean_unit_feature_set1_noise = repmat(mean_unit_feature_set1, length(noise_pts), 1);
            cov_noise_feature_set = cov(noise_feature_set);
            Mahal_distance_noise_matrix = (noise_feature_set - rep_mean_unit_feature_set1_noise) * pinv(cov_unit_feature_set1) * (noise_feature_set - rep_mean_unit_feature_set1_noise)'; %this is a 'number_spikes' x 'number_spikes' matrix with M-values along diagonal
            Mahal_distance_noise_vec = diag(Mahal_distance_noise_matrix);
            [Mahal_distance_noise_vec_sort, Mahal_noise_id] = sort(Mahal_distance_noise_vec);
            
            L_ones_noise_vec = ones(size(noise_pts));
            L_noise_vec = L_ones_noise_vec - chi2cdf(Mahal_distance_noise_vec, 12);
            L = sum(L_noise_vec(:));
            L_ratio = L / length(unit_pts{i});
            output_struct.(field).L_ratio = L_ratio;
        end
    end
end


clear global cur_coeff_pcs
clear global cur_coeff_pcs_c
                
user_output_var = inputdlg('Enter name of variable you would like to save output structure as: ', 'Save Output', [1 60]);
assignin('base', user_output_var{1}, output_struct);