function hide_clusters_or_spikes(edit_cl, Channel_Scatter, Time_Scatter, colors, val, val2)

%Description: This .m file is called via the 'Hide Cluster' or 'Hide Spikes' toggle button in the main GUI, and depending on the toggle state, either 
%1) hides spikes(s) belonging to a paritcular cluster(s) on the 'Time_Scatter' axis by setting the visibility of the corresponding scatter plot to
%'off' (via 'Hide Cluster' button) or sets the color of the coresponding scatter plot to white (via 'Hide Spikes') OR 
%2) reveals cluster(s) that have been hidden
%
%Input: 'edit_cl' = vector containing info on which cluster(s) to hide; 'Channel_Scatter' = axis; 'val' = toggle state of 'Hide Cluster' button; 
%'val2' = toggle state of 'Hide Spikes' button; 'Time_Scatter' = axis; 'colors' =  matrix of RGB-triplet colors
%

to_hide = find(edit_cl); %cluster indices to hide
if to_hide(to_hide == 14) %special case for "All Spikes"
    to_hide(to_hide == 14) = 0;
end
if to_hide(to_hide == 13) %special case for "Brushed Data"
    error('Cannot hide brushed spikes. Suggestion: set "brush" tool to white, and manually select those spikes you wish to hide.')
end

all_chan_scats = get(Channel_Scatter, 'children');
all_time_scats = get(Time_Scatter, 'children');

if max(to_hide) > (length(all_chan_scats) - 1)
    error('You have selected cluster(s) that do not exist')
end
if isempty(edit_cl)
    error('No cluster(s) selected')
end

%% Hide Clusters

if ~isempty(val)
    
    if val
        for i = 1:length(to_hide)
            if to_hide(i) == 0
                set(all_chan_scats(end), 'Visible', 'off')
                set(all_time_scats(end), 'Visible', 'off')
            else
                set(all_chan_scats(end-to_hide(i)), 'Visible', 'off')
                set(all_time_scats(end-to_hide(i)), 'Visible', 'off')
            end
        end
    end
    
    if ~val
        for i = 1:length(to_hide)
            if to_hide(i) == 0
                set(all_chan_scats(end), 'Visible', 'on')
                set(all_time_scats(end), 'Visible', 'on')
            else
                set(all_chan_scats(end-to_hide(i)), 'Visible', 'on')
                set(all_time_scats(end-to_hide(i)), 'Visible', 'on')
            end
        end
    end
    
end

%% Hide Spikes

if ~isempty(val2)
    
    if val2
        for i = 1:length(to_hide)
            if to_hide(i) == 0
                set(all_chan_scats(end), 'Visible', 'off')
                set(all_time_scats(end), 'Visible', 'off')
            else
                set(all_chan_scats(end-to_hide(i)), 'CData', [1 1 1])
                set(all_time_scats(end-to_hide(i)), 'CData', [1 1 1])
            end
        end
    end
    
    if ~val2
        for i = 1:length(to_hide)
            if to_hide(i) == 0
                set(all_chan_scats(end), 'Visible', 'on')
                set(all_time_scats(end), 'Visible', 'on')
            else
                set(all_chan_scats(end-to_hide(i)), 'CData', colors(to_hide(i),:))
                set(all_time_scats(end-to_hide(i)), 'CData', colors(to_hide(i),:))
            end
        end
    end
    
end

