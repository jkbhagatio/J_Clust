function unit_pts = run_fcm(feature, N_c)

%Description: This .m file runs Fuzzy C-means clustering to assign all detected spikes to the cluster the algorithm deems most probable
%
%Input: 'feature' = spike feature used to cluster on (across all four wires), 'N_c' = user-selected number of clusters
%
%Output: 'unit_pts' = cell array, where each cell contains the points assigned to that cluster

[c, U] = fcm(feature', N_c, [2 100 .00001 0]); %see the doc for 'fcm' to understand these parameters
[P_member, cluster_member] = max(U,[],1);

unit_pts = cell(N_c, 1);

for i = 1:N_c
    unit_pts{i} = find(cluster_member == i);
end
