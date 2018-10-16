function cluster_info = update_info_table(unit_pts, ts, overlaps, cluster_info)

%Description: This .m file updates the 'cluster_info' table in the main GUI
%
%Input: 'unit_pts' = cell array, where each cell contains the points assigned to that cluster; 'ts' = spike timestamps; 'overlaps' = spikes with
%temporally (< 1.5 ms apart) overlapped waveforms; 'cluster_info' = table in main GUI showing number of spikes in each cluster, number of overlapped 
%waveforms in each cluster, shared spike assignments between clusters, and the respective clusters of those shared spike assignments 

%Output: 'cluster_info' = updates 'cluster_info'

updated_cluster_info = cell(13,4);

updated_cluster_info{1,1} = length(ts);
updated_cluster_info{1,2} = length(overlaps);

shared_spikes = cell(length(unit_pts), 1);
shared_cluster = cell(length(unit_pts), 1);

parfor i = 1:length(unit_pts)
    for j = 1:length(unit_pts)
        if i == j
            continue;
        else
            cur_shared_spikes = intersect(unit_pts{i}, unit_pts{j});
            if ~isempty(cur_shared_spikes)
                shared_spikes{i} = horzcat(shared_spikes{i}, cur_shared_spikes);
                shared_cluster{i} = horzcat(shared_cluster{i}, j);
            end
        end
    end
end


for i = 2:length(unit_pts)+1 %entries in table corresponding to units
    updated_cluster_info{i,1} = length(unit_pts{i-1});
    updated_cluster_info{i,2} = length(intersect(unit_pts{i-1}, overlaps)); 
    updated_cluster_info{i,3} = length(shared_spikes{i-1});
    updated_cluster_info{i,4} = num2str(shared_cluster{i-1}); 
end

cluster_info.Data = updated_cluster_info;
end