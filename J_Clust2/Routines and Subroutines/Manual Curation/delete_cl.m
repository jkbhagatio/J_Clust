function unit_pts = delete_cl(unit_pts, edit_cl)

%Description: This .m file deletes a user selected cluster(s)
%
%Input: 'unit_pts' = points belonging to all current units, 'edit_cl' = vector containing info on which clusters to merge, 
%
%Output: 'unit_pts' = updated 'unit_pts'
%

cl_to_delete = find(edit_cl);

if max(cl_to_delete) > 12
    error('Cannot Delete selected spikes because they are not a cluster. Use "remove" instead.');
end

if ~sum(edit_cl(:))
    error('No cluster selected')
end

for i = 1:length(cl_to_delete) %start at 2 b/c we already have the first unit pts
    unit_pts{cl_to_delete(i)} = [];
end

end