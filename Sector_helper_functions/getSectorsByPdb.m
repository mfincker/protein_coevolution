function [sectorList, indexInSectorDB] = getSectorsByPdb(sectorDB, pdbId)
%GETSECTORSBYPDB takes a cell array containing a database of sector 
%and a pdbId . It returns a cell array containing sectors 
%whose Pdb attribute is pdbId.
%   GETSECTORSBYPDB returns an empty array if no sector has the 
%   pdbId argument as Pdb property.

    sectorList = {};
    indexInSectorDB = [];
    for i = 1:numel(sectorDB)
        if strcmp(sectorDB{i}.Pdb, pdbId) == 1
            sectorList = [sectorList sectorDB(i)];
            indexInSectorDB = [indexInSectorDB i];
        end
    end
end