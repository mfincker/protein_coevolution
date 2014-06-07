function subgroupList = getSubgroup(sectorDB)
%GETSUBGROUP takes a cell array containing a database of sector as input.
%It creates a set (cell array with no duplicate) containing the Subgroup 
%properties of the sectors in the sector database. 
%   GETSUBROUP returns an empty array if no sector has the subgroup
%   argument as Subgroup property.

    subgroupList = {};
    for i = 1:numel(sectorDB)
        if sum(strcmp(sectorDB{i}.Subgroup, subgroupList)) == 0
            subgroupList = [subgroupList sectorDB{i}.Subgroup];
        end
    end
end
    
    