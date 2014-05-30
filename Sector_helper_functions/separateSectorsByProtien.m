function [ proteinsectors] = separateSectorsByProtien( sectorDB )
% This function will separate protein sectors in a sector database by
% protein and outputs a cell array of sector DBs.
%   Detailed explanation goes here

%create open array to store list of redundant pdbs
pdbList = {};

% get redundant pdbs
for i = 1:numel(sectorDB)
    pdbList{end+1}= sectorDB{i}.Pdb;
end

[pdbList, ~ ,pdbindex]= unique(pdbList);

proteinsectors= cell(1,length(pdbList));

% create empty cell of sectors
for i = 1:length(proteinsectors)
        proteinsectors{i}= {};
end


for i = 1:length(sectorDB);
    
    currentindex = pdbindex(i);
    
    proteinsectors{currentindex}(end+1) = sectorDB(i);
end


end

