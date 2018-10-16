function [RD, CD, order] = optics_JClust(feature, minpts)

% Description:
%
% Ordering Points of a data set To Identify the Clustering Structure. To understand algorithm, see References below.
% 
% Input: 
% 'feature' = dataset(m,n): m = observations, n = variables
% 'minpts' = number of points in a neighborhood of the selected point
% (minimal number of points considered as a cluster)
% 
% Output: 
% 'RD' - vector with reachability distances for each point
% 'CD' - vector with core distances for each point
% 'order' - vector specifying the order of points as determined by minimum spanning tree
%
% References: 
% [1] M. Ankrest, M. Breunig, H. Kriegel, J. Sander (1999) OPTICS: Ordering Points To Identify the Clustering Structure
%

[m, n] = size(feature);
cov_x = inv(cov(feature));
CD = zeros(1,m);
RD = ones(1,m) * 10e9; %set all Reachability Distances very high b/c we'll eventually be taking min distances at each point
cov_flag = 0;

% Calculate Core Distances

try
    CD_all_vec = pdist(feature, 'mahalanobis', cov_x);
catch ME
    if strcmp(ME.identifier, 'stats:pdist:SingularCov')
        cov_flag = 1;
    end
end

if cov_flag
    CD_all_vec = pdist(feature);
end
CD_mtx = squareform(CD_all_vec);
sort_CD_mtx = sort(CD_mtx);
CD = sort_CD_mtx(minpts+1,:);


order = zeros(m,1);
seeds = [1:m];

min_RD_indx = 1;
ord_count = 1;

while ~isempty(seeds)
    cur_seed = seeds(min_RD_indx); %set 'cur_seed' to previously found index 'ind'       
    seeds(min_RD_indx) = []; %remove this index from 'seeds'
    order(ord_count) = cur_seed; %concatenate order in terms of closest reachability distance
    
    CDs_remaining_pts = ones(1,length(seeds)) * CD(cur_seed);
    if ~cov_flag
        dist2remaining_pts = pdist2(feature(cur_seed,:), feature(seeds,:), 'mahalanobis', cov_x); 
    else
        dist2remaining_pts = pdist2(feature(cur_seed,:), feature(seeds,:));
    end
    
    cur_RDs = max(CDs_remaining_pts, dist2remaining_pts); %Reachability distances for 'cur_seed' (max b/w core dist for 'cur_seed', & dist b/w 'cur_seed' and all other pts in 'x' that have not yet been visited) 
    old_RD_indxs = (RD(seeds)) > cur_RDs; %checks which previously set and previously minimum calculated reachability distances are greater than current reachability distances for 'cur_seed'
    RD(seeds(old_RD_indxs)) = cur_RDs(old_RD_indxs); %for all those reachability distances that are greater than current 'cur_seed' reachability distances, set them equal to current
    [min_RD_indx_val, min_RD_indx] = min(RD(seeds)); %find index of minimum reachability distance from 'cur_seed'
    ord_count = ord_count + 1;
end

RD(1) = max(RD(2:m)) + .1*max(RD(2:m)); %set RD of 1st point to be 10% higher than previous max RD since it is otherwise undefined

end
