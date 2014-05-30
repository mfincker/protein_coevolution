% Script to consoidate all membraneSectorDB files from Leah and me

membraneSectorDB = {};

% Consolidate in a single cell array all sectors:
fileList_Leah = dir('membraneSectorDB_Leah_*');
fileList_Maeva = dir('uMembraneSectorDB_Maeva_*');

for i = 1:numel(fileList_Leah)
    partialMembraneSectorDB = importdata(fileList_Leah(i).name);
    membraneSectorDB = [membraneSectorDB partialMembraneSectorDB];
end

for j = 1:numel(fileList_Maeva)
    partialMembraneSectorDB = importdata(fileList_Maeva(j).name);
    for sector = 1:numel(partialMembraneSectorDB)
        membraneSectorDB = [membraneSectorDB {partialMembraneSectorDB(sector)}];
    end
end

% Remove undefined sectors:
membraneSectorDB_undefined = {};
undefined_index = [];
for k = 1:numel(membraneSectorDB)
    if strcmp(membraneSectorDB{k}.Sequence, 'undefined') == 1
        membraneSectorDB_undefined = [membraneSectorDB_undefined membraneSectorDB(k)];
        undefined_index = [undefined_index k]
    end
end
membraneSectorDB(undefined_index) = [];


% Remove sectors from small MSA:
small_L = importdata('../membraneMSAsize/new_uMembraneSmallMSA_Leah.mat');
small_M = importdata('../membraneMSAsize/new_uMembraneSmallMSA_Maeva.mat');

disp(['number of small MSA : ' num2str(numel(small_L) + numel(small_M))]);

for l = 1:numel(small_L)
    disp(['Small MSA : ' small_L{l}]);
    [sectors, index] = getSectorsByPdb(membraneSectorDB, small_L{l});
    disp([char(9) 'number of sectors with that Pdb id : ' num2str(numel(sectors))]);
    membraneSectorDB(index) = [];
end

for m = 1:numel(small_M)
    disp(['Small MSA : ' small_M{m}]);
    [sectors, index] = getSectorsByPdb(membraneSectorDB, small_M{m});
    disp([char(9) 'number of sectors with that Pdb id : ' num2str(numel(sectors))]);
    membraneSectorDB(index) = [];
end



save('membraneSectorDB_052914.mat', 'membraneSectorDB'); 



