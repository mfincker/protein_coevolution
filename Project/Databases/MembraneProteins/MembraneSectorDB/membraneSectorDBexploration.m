% Playing with the new membrane sector DB:

sectorLength = zeros(numel(membraneSectorDB),1);
for i = 1:numel(membraneSectorDB)
    sectorLength(i) = membraneSectorDB{i}.Length;
end

% Distribution of the sector lengths:
figure;
hist(sectorLength,50);
title('Distribution of the sector length');

% Length delimited sector DB :
membraneSectorDB_over10 = {};
membraneSectorDB_over100 = {};
membraneSectorDB_under100 = {};

for i = 1:numel(membraneSectorDB)
    if membraneSectorDB{i}.Length >= 10 
        membraneSectorDB_over10 = [membraneSectorDB_over10 {membraneSectorDB{i}}];
    end
    if membraneSectorDB{i}.Length >= 100 
        membraneSectorDB_over100 = [membraneSectorDB_over100 {membraneSectorDB{i}}];
    end
    if membraneSectorDB{i}.Length < 100 
        membraneSectorDB_under100 = [membraneSectorDB_under100 {membraneSectorDB{i}}];
    end
end

save('membraneSectorDB_over10_Maeva_052714.mat', 'membraneSectorDB_over10');
save('membraneSectorDB_over100_Maeva_052714.mat', 'membraneSectorDB_over100');
save('membraneSectorDB_under100_Maeva_052714.mat', 'membraneSectorDB_under100');

