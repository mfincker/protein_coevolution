directory = './prokaryote_Sectors/' ;
d = dir( './prokaryote_Sectors/prokaryote_MSA_SectorDatabase*.mat');
sectorDatabaseCopy = {};
for i = 1:numel(d)
    sectorDatabaseCopy = [sectorDatabaseCopy importdata( [directory d(i).name] ) ];
end 

sectorDatabase = {};
for j = 1:numel(sectorDatabaseCopy)
    if strcmp(sectorDatabaseCopy{j}.Sequence, 'undefined') == 0
        sectorDatabase = [sectorDatabase {sectorDatabaseCopy{j}}];
    end
end
save('./prokaryote_Sectors/prokaryoticSectorDB.mat','sectorDatabase');