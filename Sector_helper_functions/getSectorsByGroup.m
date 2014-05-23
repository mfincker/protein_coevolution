function sectorList = getSectorsByGroup(sectorDB, group)
%GETSECTORSBYGROUP takes a cell array containing a database of sector 
%and a group name. It returns a cell array containing sectors 
%whose Group attribute is group.
%   GETSECTORSBYGROUP returns an empty array if no sector has the group
%   argument as Group property.

    sectorList = {};
    indexInSectorDB = [];
    for i = 1:numel(sectorDB)
        if strcmp(sectorDB{i}.Group, group) == 1
            sectorList = [sectorList sectorDB(i)];
            indexInSectorDB = [indexInSectorDB i];
        end
    end
end
    
    
    
