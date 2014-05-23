function groupList = getGroup(sectorDB)
%GETGROUP takes a cell array containing a database of sector as input.
%It creates a set (cell array with no duplicate) containing the Group 
%properties of the sectors in the sector database. 
%   GETGROUP returns an empty array if no sector has the subgroup
%   argument as Subgroup property.

    groupList = {};
    for i = 1:numel(sectorDB)
        if sum(strcmp(sectorDB{i}.Group, groupList)) == 0
            groupList = [groupList sectorDB{i}.Group];
        end
    end
end
    
    