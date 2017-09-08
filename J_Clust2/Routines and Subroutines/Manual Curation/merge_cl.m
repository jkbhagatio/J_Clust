function unit_pts = merge_cl(unit_pts, edit_cl)

%Description: This .m file merges two (or more) user selected clusters
%
%Input: 'unit_pts' = points belonging to all current units, 'edit_cl' = vector containing info on which clusters to merge, 
%
%Output: 'unit_pts' = updated 'unit_pts'
%

cl_to_merge = find(edit_cl);

if max(cl_to_merge) > 12
    error('Cannot Merge selected spikes because they are not a cluster. Use "add" instead.');
end

if length(find(edit_cl)) < 2
    error('Not enough clusters selected')
end

% merge units

for i = 2:length(cl_to_merge) %start at 2 b/c we already have the first unit pts
    unit_pts{min(cl_to_merge)} = [unit_pts{min(cl_to_merge)}, unit_pts{cl_to_merge(i)}];
    unit_pts{cl_to_merge(i)} = [];
end



end