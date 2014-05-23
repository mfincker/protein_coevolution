function sectorList = getSectorsBySubgroup(sectorDB, subgroup)
%GETSECTORSBYSUBGROUP takes a cell array containing a database of sector 
%and a subgroup name. It returns a cell array containing sectors 
%whose Subgroup attribute is subgroup.
%   GETSECTORSBYSUBGROUP returns an empty array if no sector has the 
%   subgroup argument as Subgroup property.

    sectorList = {};
    indexInSectorDB = [];
    for i = 1:numel(sectorDB)
        if strcmp(sectorDB{i}.Subgroup, subgroup) == 1
            sectorList = [sectorList sectorDB(i)];
            indexInSectorDB = [indexInSectorDB i];
        end
    end
end
    
    