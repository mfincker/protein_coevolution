function pdbList = getPdb(sectorDB)
%This is NOT the Matlab getpdb() function !!!!
%See 'help getpdb' to learn more about the Matlab function.
%GETPDB takes a cell array containing a database of sector as input.
%It creates a set (cell array with no duplicate) containing the Pdb 
%properties of the sectors in the sector database. 
%   GETPDB returns an empty array if no sector has the pdbId
%   argument as Pdb property.

    pdbList = {};
    for i = 1:numel(sectorDB)
        if sum(strcmp(sectorDB{i}.Pdb, pdbList)) == 0
            pdbList = [pdbList sectorDB{i}.Pdb];
        end
    end
end