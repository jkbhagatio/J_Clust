function unit_pts = remove_pts_from_cl(unit_pts, selected_pts, edit_cl)

%Description: This .m file removes spikes to a user selected cluster
%
%Input: 'unit_pts' = points belonging to all current units, 'selected_pts' = user selected points from main GUI, 
%'edit_cl' = vector containing info on which cluster to add points to
%
%Output: 'unit_pts' = updated 'unit_pts'
%

if isempty(selected_pts)
    error('No points selected')
end

if ~sum(edit_cl(:))
    error('No cluster selected')
end

units_to_remove_from = find(edit_cl);

if max(units_to_remove_from) > 12
    error('Cannot remove spikes from "Selected Spikes" or "All Non-Clustered Spikes"');
end

for i = 1:length(units_to_remove_from)
    unit_removed_from = unit_pts{units_to_remove_from(i)};
    if ~isrow(unit_removed_from)
        unit_removed_from = unit_removed_from';
    end
    if ~isrow(selected_pts)
        selected_pts = selected_pts';
    end
    [found_pts, remove_pts, ~] = intersect(unit_removed_from, selected_pts);
    unit_removed_from(remove_pts) = [];
    unit_pts{units_to_remove_from(i)} = unit_removed_from;
end
    
end