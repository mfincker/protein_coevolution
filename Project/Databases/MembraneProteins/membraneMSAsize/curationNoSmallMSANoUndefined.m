%% Curation of the sectorDB

% Remove sectors from small MSA (< 50 seauences)
uMembraneSectorDB_12_curated = importdata('/Users/Maeva/Documents/Stanford/Stanford courses/Spring 2014/BIOE 115/proteinCoevolution/protein_coevolution/MembraneProteins/uMembraneSectorDB/uMembraneSectorDB_11GroupSubMPRP.mat');
disp(['initial size of the database : ' num2str(numel(uMembraneSectorDB_12_curated))]);

smallMSA = importdata('/Users/Maeva/Documents/Stanford/Stanford courses/Spring 2014/BIOE 115/proteinCoevolution/protein_coevolution/MembraneProteins/uMembraneSmallMSA.mat');
disp(['number of small MSA : ' num2str(numel(smallMSA))]);

for i = 1:numel(smallMSA)
    disp(['Small MSA : ' smallMSA{i}]);
    [sectors, index] = getSectorsByPdb(uMembraneSectorDB_12_curated, smallMSA{i});
    disp([char(9) 'number of sectors with that Pdb id : ' num2str(numel(sectors))]);
    uMembraneSectorDB_12_curated(index) = [];
end

disp (['Size of the database after removing small MSAs : ' num2str(numel(uMembraneSectorDB_12_curated))]);

% Remove undefined sectors

uMembraneSectorDB_12_undefined = {};

for i = numel(uMembraneSectorDB_12_curated)
    if strcmp(uMembraneSectorDB_12_curated{i}.Sequence, 'undefined') == 1
        uMembraneSectorDB_12_undefined = [uMembraneSectorDB_12_undefined uMembraneSectorDB_12_curated(i)];
        uMembraneSectorDB_12_curated(i) = [];
        
    end
end

disp (['Final of the database after removing undefined sectors : ' num2str(numel(uMembraneSectorDB_12_curated))]);


save('/Users/Maeva/Documents/Stanford/Stanford courses/Spring 2014/BIOE 115/proteinCoevolution/protein_coevolution/MembraneProteins/uMembraneSectorDB/uMembraneSectorDB_12_noSmallMSAnoUndefined.mat', 'uMembraneSectorDB_12_curated');
save('/Users/Maeva/Documents/Stanford/Stanford courses/Spring 2014/BIOE 115/proteinCoevolution/protein_coevolution/MembraneProteins/uMembraneSectorDB/uMembraneSectorDB_12_undefined.mat', 'uMembraneSectorDB_12_undefined');




