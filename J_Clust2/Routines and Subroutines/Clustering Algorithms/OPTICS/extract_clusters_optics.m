function clusters = extract_clusters_optics(RD, order, minpts)
% Description: Extract Flat Clusters from a Reachability-Distance Plot
%
%This code finds a Global Maximum, splits all points into nodes left and right of the global maximum, and finds local maxima within these two nodes. 
%We assume each local maxima is a potential split point between two clusters, and test each maxima to see if it passes to be considered a split point. 
%If the maxima is deemed significant, we create clusters to the left and right of the maxima, as long as there are 'minpts' between the maxima and the 
%next maxima, or between the maxima and the first or last index.
%
% Inputs:
%
% 'RD' - Reachability-Distance of each point to the next nearest point
% 'order' - A vector containing all points, where each subsequent point in
%vector is the index of the next closest point
% 'minpts' - minimum number of points to create a cluster
%
% Output:
% 'clusters' - a cell array, with each cell containing the indices of a
%cluster
%
%References: 
% [1] J. Sander et al. Automatic Extraction of Clusters from Heirarchical Clustering Representations (2003)
%

clusters = [];

cl_right = cell(1);
cl_left = cell(1);

ord_RD = RD(order); %put RD in order of RD plot
m = length(RD);

[maxima_vals, maxima_indxs] = findpeaks(ord_RD);

[g_max, indx] = max(maxima_vals); %global maxima
g_max_indx = maxima_indxs(indx);

if isempty(g_max)
    return;
end

sig_pct = .85; %percentage of maxima which must be > than surrounding areas

if (g_max * sig_pct < mean(ord_RD(1:g_max_indx))) & (g_max * sig_pct < mean(ord_RD(g_max_indx+1:g_max_indx+1+minpts))) %return 0 clusters if global maxima is not significant
    return;
end

%Organize into clustering on right side of global maxima, then left side of global maxima

%right-side clusters

cur_pts = [g_max_indx:m];
[right_maximas, right_maximas_indxs] = findpeaks(ord_RD(cur_pts));
if isempty(right_maximas) %if empty, rest of R-plot is a cluster
    [min_valley, min_valley_indx] = findpeaks(-ord_RD(cur_pts), 'NPeaks', 1); %find first local minima
    min_valley = min_valley * -1;
    cl_indxL = g_max_indx + find(ord_RD(cur_pts) < ((g_max+min_valley)* .6), 1);
    cl_indxR = cl_indxL + find(ord_RD(cl_indxL:m) > g_max,1) - 1;
    if (cl_indxR - cl_indxL) >= minpts 
        cl_right{1} = order(cl_indxL:cl_indxR);
    end
else
    right_maximas_indxs = right_maximas_indxs + g_max_indx -1;
    cur_max = g_max;
    cur_max_indx = g_max_indx;
    %for each right maxima, check to see if a possible cluster to the right is significant, if so, then cluster to the left of current maxima to previous maxima
    for i = 1:length(right_maximas)
        prev_max = cur_max;
        prev_max_indx = cur_max_indx;
        cur_max = right_maximas(i);
        cur_max_indx = right_maximas_indxs(i);
        cur_pts = [prev_max_indx:cur_max_indx];
        
        if (cur_max * sig_pct < mean(ord_RD(cur_max_indx:cur_max_indx+minpts+1))) %if current maxima is not signifcant to the right of maxima
            cur_max = prev_max;
            cur_max_indx = prev_max_indx;
            continue
        else %else create a cluster from previous maxima to current maxima (ends at next, previous maxima value, or at height of current maxima)
            %find start and end of cluster
            [min_valley, min_valley_indx] = findpeaks(-ord_RD(cur_pts), 'NPeaks', 1); %find first local minima
            min_valley = min_valley * -1;
            cl_indxL = prev_max_indx + find(ord_RD(cur_pts) < ((prev_max+min_valley)* .6), 1);
            cl_indxR = cl_indxL + find(ord_RD(cl_indxL:cur_max_indx) > prev_max,1) - 1;
            if isempty(cl_indxR)
                cl_indxR = cur_max_indx;
            end
            if (cl_indxR - cl_indxL) >= minpts
                cl_right{i} = order(cl_indxL:cl_indxR);
            end
        end
    end
    
    if isempty(cl_right) %if none of the right-maxima were significant
        %create a cluster from global maxima
        cur_pts = [g_max_indx:m];
        [min_valley, min_valley_indx] = findpeaks(-ord_RD(cur_pts), 'NPeaks', 1); %find first local minima
        min_valley = min_valley * -1;
        cl_indxL = g_max_indx + find(ord_RD(cur_pts) < ((g_max+min_valley)* .6), 1);
        cl_indxR = cl_indxL + find(ord_RD(cl_indxL:m) > g_max,1) - 1;
        if (cl_indxR - cl_indxL) >= minpts
            cl_right{1} = order(cl_indxL:cl_indxR);
        end
    else %add far right cluster from final right maxima to end of RD plot
        cur_pts = [cur_max_indx:m];
        [min_valley, min_valley_indx] = findpeaks(-ord_RD(cur_pts), 'NPeaks', 1); %find first local minima
        min_valley = min_valley * -1;
        cl_indxL = cur_max_indx + find(ord_RD(cur_pts) < ((cur_max+min_valley) *.6),1);
        cl_indxR = cl_indxL + find(ord_RD(cl_indxL:m) > cur_max,1) - 1;
        if (cl_indxR - cl_indxL) >= minpts
            cl_right{length(cl_right) + 1} = order(cl_indxL:cl_indxR);
        end
    end
end

%left-side clusters

cur_pts = [1:g_max_indx];
[left_maximas, left_maximas_indxs] = findpeaks(ord_RD(cur_pts));
if isempty(left_maximas) %if empty, rest of R-plot is a cluster
    [min_valley, min_valley_indx] = findpeaks(-ord_RD(cur_pts), 'NPeaks', 1); %find first local minima
    min_valley = min_valley * -1;
    cl_indxL = find(ord_RD(cur_pts) < ((g_max+min_valley)* .6), 1);
    cl_indxR = cl_indxL + find(ord_RD(cl_indxL:g_max_indx) > ((g_max+min_valley)/2),1);
    if (cl_indxR - cl_indxL) >= minpts 
        cl_left{1} = order(cl_indxL:cl_indxR);
    end
else
    cur_max = g_max;
    cur_max_indx = g_max_indx;
    left_maximas = flip(left_maximas); left_maximas_indxs = flip(left_maximas_indxs); %flip maximas and indxs to start clustering from g_max down to index 1
    %for each left maxima, check to see if a possible cluster to the left is significant, if so, then cluster to the right of current maxima from previous maxima
    for i = 1:length(left_maximas)
        prev_max = cur_max;
        prev_max_indx = cur_max_indx;
        cur_max = left_maximas(i);
        cur_max_indx = left_maximas_indxs(i);
        cur_pts = [cur_max_indx:prev_max_indx];
        
        if (cur_max * sig_pct < mean(ord_RD(cur_max_indx-minpts-1:cur_max_indx))) %if current maxima is not signifcant to the left of maxima
            cur_max = prev_max;
            cur_max_indx = prev_max_indx;
            continue
        else %else create a cluster from current maxima to previous maxima (ends at previous maxima value, or at height of current maxima)
            %find start and end of cluster
            [min_valley, min_valley_indx] = findpeaks(-ord_RD(cur_pts), 'NPeaks', 1); %find first local minima
            min_valley = min_valley * -1;
            cl_indxL = cur_max_indx + find(ord_RD(cur_pts) < ((cur_max+min_valley)* .6), 1);
            cl_indxR = cl_indxL + find(ord_RD(cl_indxL:prev_max_indx) > cur_max,1) - 1;
            if isempty(cl_indxR)
                cl_indxR = prev_max_indx;
            end
            if (cl_indxR - cl_indxL) >= minpts
                cl_left{i} = order(cl_indxL:cl_indxR);
                prev_max = cur_max;
                prev_max_indx = cur_max_indx;
            end
        end
    end
    
    if isempty(cl_left) %if none of the left-maxima were significant
        %create a cluster from global maxima
        cur_pts = [1:g_max_indx];
        [min_valley, min_valley_indx] = findpeaks(-ord_RD(cur_pts), 'NPeaks', 1); %find first local minima
        min_valley = min_valley * -1;
        cl_indxL = find(ord_RD(cur_pts) < ((g_max+min_valley)* .6), 1);
        cl_indxR = cl_indxL + find(ord_RD(cl_indxL:g_max_indx) > ((g_max+min_valley)/2),1);
        if (cl_indxR - cl_indxL) >= minpts
            cl_left{1} = order(cl_indxL:cl_indxR);
        end
    else %add initial cluster
        for i = 2:prev_max_indx
            if ord_RD(i)/ord_RD(i+1) < 1.75
                cl_indxL = i+1;
                cur_max = ord_RD(i+1);
                break
            end
        end
        
        cl_indxR = cl_indxL + find(ord_RD(cl_indxL:prev_max_indx) > cur_max*1.5, 1);
        if isempty(cl_indxR)
            cl_indxR = prev_max_indx;
        end
        if (cl_indxR - cl_indxL) >= minpts
            cl_left{length(cl_left) + 1} = order(cl_indxL:cl_indxR);
        end
    end
end

%remove empties for insignificant maxima
cl_left = cl_left(~cellfun('isempty', cl_left));
cl_right = cl_right(~cellfun('isempty', cl_right));

clusters = [cl_left, cl_right];

if length(clusters) > 12 %prune clusters if >12
    for i = 13:length(clusters)
        clusters{i} = [];
    end
    clusters = clusters(~cellfun('isempty', clusters));
end

end
    
    
