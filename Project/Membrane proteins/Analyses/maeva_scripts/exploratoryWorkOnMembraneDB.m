%% Exploring the membrane sector DB

clc; close all;

% Importing the database
DB = importdata('../MembraneSectorDB/membraneSectorDB_060114_forCoord.mat'); 

%% Some statistics about the database:

%%  Protein statistics

% Number of proteins:
pdb = getPdb(DB);
numPdb = numel(pdb);

% Number of groups:
groups = getGroup(DB);
numGroup = numel(groups);

% Number of proteins per group:
proteinPerGroup = zeros(1, numGroup);
for group = 1:numGroup
    sectors = getSectorsByGroup(DB, groups{group});
    proteinPerGroup(1,group) = numel(getPdb(sectors));
end

%%
figure('Name','Distribution of the proteins per group');
bar(proteinPerGroup);
set(gca,'XTickLabel',char(groups));

%%
% Number of subgroups:
subgroups = getSubgroup(DB);
numSubgroup = numel(subgroups);

% Number of proteins per subgroup:
proteinPerSubgroup = zeros(1, numSubgroup);
for subgroup = 1:numSubgroup
    sectors = getSectorsBySubgroup(DB, subgroups{subgroup});
    proteinPerSubgroup(1,subgroup) = numel(getPdb(sectors));
    [proteinPerSubgroup_sorted, permutations] = sort(proteinPerSubgroup);
end

%%
figure('Name','Distribution of the proteins per subgroup');
bar(proteinPerSubgroup_sorted);
set(gca,'XTickLabel',char(subgroups(permutations)));

%% Sector statistics

% Number of sectors:
numSector = numel(DB);

% Number of sector per protein:
sectorPerProtein = zeros(1, numPdb);
for prot = 1:numPdb
    sectors = getSectorsByPdb(DB, pdb{prot});
    sectorPerProtein(1, prot) = numel(sectors);
end

meanSectorPerProtein = mean(sectorPerProtein);

%%
figure('Name','Distribution of the number of sector per protein');
subplot(1,2,1);
boxplot(sectorPerProtein);
title('Distribution of the number of sector per protein');
ylabel('Number of sectors');
set(gca,'XTickLabel','Membrane sector database');
subplot(1,2,2);
hist(sectorPerProtein);
title('Distribution of the number of sector per protein');
xlabel('Nuumber of sectors');
ylabel('Number of proteins');

%%
% Statistics on the sector length

lengths = zeros(1,numSector);
for sec = 1:numSector
    lengths(1,sec) = DB{sec}.Length;
end

meanLength = mean(lengths);
medianLength = median(lengths);

%%
figure('Name','Distribution of the length of sectors');
boxplot(lengths);
ylabel('Number of reisdues', 'FontSize',11,'FontName', 'arial');
set(gca,'XTickLabel','Membrane sector database','FontSize',11,'FontName', 'arial');



%%
% Aminoc acid frequency in the sector DB:
numResidue = sum(lengths);
countAAinDB = aaSectorDBCount(DB);
countAAinDB = sum(countAAinDB, 2);
freqAAinDB = countAAinDB / numResidue;

% %%
% % Amino acid frequency in the membrane proteins used to make the DB:
% countAAinMembraneProteins = zeros(20,numPdb);
% for prot = 1:numPdb
%     seq = getgenpept(pdb{prot});
%     seq = seq.Sequence;
%     count = aacount(seq);
% 	count = struct2cell(count);
% 	count = cell2mat(count);
%     countAAinMembraneProteins(:, prot) = count;
% end
% 
% countAAinMembraneProteins = sum(countAAinMembraneProteins,2);
% save('countAAinMembraneProteins.mat', 'countAAinMembraneProteins');

countAAinMembraneProteins = importdata('countAAinMembraneProteins.mat');
freqAAinMembraneProteins = countAAinMembraneProteins ...
                                / sum(sum(countAAinMembraneProteins,1),2);
       

% Amino acid frequency in UniprotKB/trEMBl
freqAAinUniprot = [8.73 5.38 4.11 5.34 1.19 4 6.2 7.1 2.18 6.09 ...
                10 5.29 2.50 4.03 4.55 6.52 5.52 1.29 3.06 6.81]' / 100;
            
%%
figure('Name', 'Amino acid frequency in the membrane sectors');
subplot(3,3,1);
bar(freqAAinDB);
title('Amino acid frequency in the membrane sectors')
ylabel('Frequency');
ylim([0 0.12]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);

subplot(3,3,4);
bar(freqAAinUniprot);
title('Amino acid frequency in Uniprot')
ylabel('Frequency');
ylim([0 0.12]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);

subplot(3,3,7);
bar(freqAAinDB - freqAAinUniprot);
title(['Difference in amino acid frequency:' char(10) 'sector - uniprot'])
ylabel('Frequency');
ylim([-0.02 0.02]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);

subplot(3,3,2);
bar(freqAAinDB);
ylim([0 0.12]);
title('Amino acid frequency in the membrane sectors')
ylabel('Frequency');
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);

subplot(3,3,5);
bar(freqAAinMembraneProteins);
title('Amino acid frequency in the membrane proteins')
ylabel('Frequency');
ylim([0 0.12]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);

subplot(3,3,8);
bar(freqAAinDB - freqAAinMembraneProteins);
title(['Difference in amino acid frequency:' char(10) 'sector - membrane proteins'])
ylabel('Frequency');
ylim([-0.02 0.02]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);


subplot(3,3,3);
bar(freqAAinMembraneProteins);
title('Amino acid frequency in the membrane proteins')
ylabel('Frequency');
ylim([0 0.12]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);

subplot(3,3,6);
bar(freqAAinUniprot);
title('Amino acid frequency in Uniprot')
ylabel('Frequency');
ylim([0 0.12]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);

subplot(3,3,9);
bar(freqAAinMembraneProteins - freqAAinUniprot);
title(['Difference in amino acid frequency:' char(10) 'membrane proteins - uniprot'])
ylabel('Frequency');
ylim([-0.02 0.02]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);
            

%%
figure('Name', 'Amino acid frequency in the membrane sectors');

bar([(freqAAinDB - freqAAinMembraneProteins) ./ freqAAinMembraneProteins, ...
    (freqAAinDB - freqAAinUniprot) ./ freqAAinUniprot]);
title(['Relative enrichment in the membrane sectors'], 'FontSize',20,...
       'FontWeight','bold')
ylabel('Frequency', 'FontSize',18,...
       'FontWeight','bold');
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20, 'FontSize',18,...
       'FontWeight','bold');
a = legend('relative enrichment from the membrane proteins', 'relative enrichment from Uniprot ');
set(a, 'Location', 'NorthWest', 'FontSize',18,...
       'FontWeight','bold');


 


%%
% Analysis: 
% The difference in distribution between the aa freq in the
% sectors and in uniprot is very close to the difference betwee the
% membrane proteins and uniprot. All differences go in the same direction
% except for the leucine and isoleucine : it looks like the sectors are
% depleted is branched aa instead of being enriched as they are in membrane
% proteins.
% In genral, membrane proteins are highly enriched in glycine, phenylalanine
% and tyrosine and more modestly in valine, tryptophan, threonine, leucine,
% asparagine and iso-leucine. All are non-polar (apart for Asparagine) and
% either aliphatic or aromatic. Sectors are even more enriched in 
% asparagine, glycine, tryptophan and tyrosine.

%%
% Cluster of aa freq in sectors:
% freqAAperSector = aaSectorDBFreq(DB);
% freqCluster = clustergram(freqAAperSector, 'Cluster', 'Row', 'RowPdist', 'jaccard');

% Nothing here ...

% varFreqAAperSector = freqAAperSector - repmat(freqAAinMembraneProteins,1,size(freqAAperSector,2));
% freqCluster = clustergram(varFreqAAperSector, 'Cluster', 'Row', 'RowPdist', 'jaccard');

% Nothing here ...

%% Exploring the coordinates of the sectors


[group1, index1] = getSectorsByGroup(DB, groups{1});
[group2, index2] = getSectorsByGroup(DB, groups{2});

subgroupSectors = {};
subgroupIndex = [];
for i = 1:numel(subgroups)
    [subgroupSector, subgroupI] = getSectorsBySubgroup(DB, subgroups{i});
    subgroupSectors = [subgroupSectors subgroupSector];
    subgroupIndex = [subgroupIndex {subgroupI}];
end

% Eigenvalue of the coordinates
% Run PCA, returns a 3 x num sector matrix of eigenValues 1, 2 and 3
eigenValueDB = sectorDBcoordPCA(DB, 0);

% Normalizing eigenvalue:
sumEigenValueDB = sum(eigenValueDB,1);
normEigenValueDB = eigenValueDB ./ repmat(sumEigenValueDB,3,1);

%%
% Plotting the eigenValues and noramlized ones:
figure('Name', 'Coordinate eigenValues');
subplot(4,2,1);
scatter3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:));
title('Scatter plot of the eigen values of each sector in the database')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');
subplot(4,2,3);
scatter(eigenValueDB(1,:),eigenValueDB(2,:));
title('Projection of the eigen values of each sector in the database')
xlabel('eigenValue 1');
ylabel('eigenValue 2');
subplot(4,2,5);
scatter(eigenValueDB(1,:),eigenValueDB(3,:));
title('Projection of the eigen values of each sector in the database')
xlabel('eigenValue 1');
ylabel('eigenValue 3');
subplot(4,2,7);
scatter(eigenValueDB(2,:),eigenValueDB(3,:));
title('Projection of the eigen values of each sector in the database')
xlabel('eigenValue 2');
ylabel('eigenValue 3');

subplot(4,2,2);
scatter3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:));
title('Scatter plot of the  normalized eigen values of each sector in the database')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');
subplot(4,2,4);
scatter(normEigenValueDB(1,:),normEigenValueDB(2,:));
title('Projection of the  normalized eigen values of each sector in the database')
xlabel('eigenValue 1');
ylabel('eigenValue 2');
subplot(4,2,6);
scatter(normEigenValueDB(1,:),normEigenValueDB(3,:));
title('Projection of the  normalized eigen values of each sector in the database')
xlabel('eigenValue 1');
ylabel('eigenValue 3');
subplot(4,2,8);
scatter(normEigenValueDB(2,:),normEigenValueDB(3,:));
title('Projection of the  normalized eigen values of each sector in the database')
xlabel('eigenValue 2');
ylabel('eigenValue 3');

%%
% Subgroup of interest: index in subgroups

% Electron transport:
eTransport = [1 2 4 14];

% ATPase
ATPase = [19 53 69];

% Oxygenase
oxygenase = [35 38 47];
cyclooxygenase = 32;

% Dehydrogenase
dehydrogenase = [33 36];

% Oxydoreductase
oxydoreductase = [37 70];

%%
c = [0.75 0.25 0.25 ; 0.25 0.75 0.25 ; 0.25 0.25 0.75 ; 0.75 0.75 0.25];
figure('Name', 'EigenValues of the subgroups of interest')
subplot(2,5,1);
for i =1:numel(eTransport)
    plot3(eigenValueDB(3,subgroupIndex{eTransport(i)}), ...
        eigenValueDB(1,subgroupIndex{eTransport(i)}), ...
        eigenValueDB(2,subgroupIndex{eTransport(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 15);
    hold on;
end
% plot3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:),'.', ...
%     'Color',[0.7 0.7 0.7], 'MarkerSize', 8);
ylim([0 5000]);
xlim([0 250]);
zlim([0 500]);
grid on;
title('Electron transport sector')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

subplot(2,5,6);
for i =1:numel(eTransport)
    plot3(normEigenValueDB(3,subgroupIndex{eTransport(i)}), ...
        normEigenValueDB(1,subgroupIndex{eTransport(i)}), ...
        normEigenValueDB(2,subgroupIndex{eTransport(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 15);
    hold on;
end
% plot3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:),'.', ...
%     'Color',[0.3 0.3 0.3], 'MarkerSize', 4);
xlim([0 0.3]);
ylim([.4 1]);
zlim([0 0.5]);
grid on;
title('Electron transport sector')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

%%%%%%

subplot(2,5,2);
for i =1:numel(ATPase)
    plot3(eigenValueDB(3,subgroupIndex{ATPase(i)}), ...
        eigenValueDB(1,subgroupIndex{ATPase(i)}), ...
        eigenValueDB(2,subgroupIndex{ATPase(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 15);
    hold on;
end
% plot3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:),'.', ...
%     'Color',[0.7 0.7 0.7], 'MarkerSize', 8);
ylim([0 5000]);
xlim([0 250]);
zlim([0 500]);
grid on;
title('ATPase')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

subplot(2,5,7);
for i =1:numel(ATPase)
    plot3(normEigenValueDB(3,subgroupIndex{ATPase(i)}), ...
        normEigenValueDB(1,subgroupIndex{ATPase(i)}), ...
        normEigenValueDB(2,subgroupIndex{ATPase(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 15);
    hold on;
end
% plot3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:),'.', ...
%     'Color',[0.3 0.3 0.3], 'MarkerSize', 4);
xlim([0 0.3]);
ylim([.4 1]);
zlim([0 0.5]);
grid on;
title('ATPase')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

%%%%%%

subplot(2,5,3);
for i =1:numel(ATPase)
    plot3(eigenValueDB(3,subgroupIndex{oxygenase(i)}), ...
        eigenValueDB(1,subgroupIndex{oxygenase(i)}), ...
        eigenValueDB(2,subgroupIndex{oxygenase(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 15);
    hold on;
end
% plot3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:),'.', ...
%     'Color',[0.7 0.7 0.7], 'MarkerSize', 8);
ylim([0 5000]);
xlim([0 250]);
zlim([0 500]);
grid on;
title('oxygenase')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

subplot(2,5,8);
for i =1:numel(oxygenase)
    plot3(normEigenValueDB(3,subgroupIndex{oxygenase(i)}), ...
        normEigenValueDB(1,subgroupIndex{oxygenase(i)}), ...
        normEigenValueDB(2,subgroupIndex{oxygenase(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 15);
    hold on;
end
% plot3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:),'.', ...
%     'Color',[0.3 0.3 0.3], 'MarkerSize', 4);
xlim([0 0.3]);
ylim([.4 1]);
zlim([0 0.5]);
grid on;
title('Oxygenase')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

%%%%%%

subplot(2,5,4);
for i =1:numel(dehydrogenase)
    plot3(eigenValueDB(3,subgroupIndex{dehydrogenase(i)}), ...
        eigenValueDB(1,subgroupIndex{dehydrogenase(i)}), ...
        eigenValueDB(2,subgroupIndex{dehydrogenase(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 15);
    hold on;
end
% plot3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:),'.', ...
%     'Color',[0.7 0.7 0.7], 'MarkerSize', 8);
ylim([0 5000]);
xlim([0 250]);
zlim([0 500]);
grid on;
title('oxygenase')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

subplot(2,5,9);
for i =1:numel(dehydrogenase)
    plot3(normEigenValueDB(3,subgroupIndex{dehydrogenase(i)}), ...
        normEigenValueDB(1,subgroupIndex{dehydrogenase(i)}), ...
        normEigenValueDB(2,subgroupIndex{dehydrogenase(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 15);
    hold on;
end
% plot3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:),'.', ...
%     'Color',[0.3 0.3 0.3], 'MarkerSize', 4);
xlim([0 0.3]);
ylim([.4 1]);
zlim([0 0.5]);
grid on;
title('Dehydrogenase')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

%%%%%%

subplot(2,5,5);
for i =1:numel(oxydoreductase)
    plot3(eigenValueDB(3,subgroupIndex{oxydoreductase(i)}), ...
        eigenValueDB(1,subgroupIndex{oxydoreductase(i)}), ...
        eigenValueDB(2,subgroupIndex{oxydoreductase(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 15);
    hold on;
end
% plot3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:),'.', ...
%     'Color',[0.7 0.7 0.7], 'MarkerSize', 8);
ylim([0 5000]);
xlim([0 250]);
zlim([0 500]);
grid on;
title('Oxydoreductase')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

subplot(2,5,10);
for i =1:numel(oxydoreductase)
    plot3(normEigenValueDB(3,subgroupIndex{oxydoreductase(i)}), ...
        normEigenValueDB(1,subgroupIndex{oxydoreductase(i)}), ...
        normEigenValueDB(2,subgroupIndex{oxydoreductase(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 15);
    hold on;
end
% plot3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:),'.', ...
%     'Color',[0.3 0.3 0.3], 'MarkerSize', 4);
xlim([0 0.3]);
ylim([.4 1]);
zlim([0 0.5]);
grid on;
title('Oxydoreductase')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

% %%
% h = clustergram(normEigenValueDB,'Cluster' , 'Row', 'Rowpdist', 'correlation', 'Colormap', jet);
% 
% %%
% hATPase = clustergram(eigenValueDB,'Cluster' , 'Row', 'Rowpdist', 'correlation', 'Colormap', jet);
% 
% %%
% hATPase = clustergram(normEigenValueDB(:,[subgroupIndex{ATPase}]), 'Colormap', jet);
% 
% %%
% heTransport = clustergram(normEigenValueDB(:,[subgroupIndex{eTransport}]), 'Colormap', jet);
% 
% %%
% hOxygenase = clustergram(normEigenValueDB(:,[subgroupIndex{oxygenase}]), 'Colormap', jet);
% %Maybe something interesting:
% oxygenaseSector = DB([subgroupIndex{oxygenase}]);
% % Sector 22, 35 and 25 share the same eigen value distribution.
% %%
% hDehydrogenase = clustergram(normEigenValueDB(:,[subgroupIndex{dehydrogenase}]), 'Colormap', jet);
% 
% %%
% hOxydoreductase = clustergram(normEigenValueDB(:,[subgroupIndex{oxydoreductase}]), 'Colormap', jet);

%% Analysis of the sectors dtrectched in a single direction

%%
% Identification of sectors where normEigenVal1 > 0.9 and sector length > 5
stretchedIndex = find(normEigenValueDB(1,:) >= 0.9 );
% 264 sectors picked
notLongEnough = [];
for i = 1:numel(stretchedIndex)
    if DB{stretchedIndex(i)}.Length < 5
        notLongEnough = [notLongEnough i];
    end
end

stretchedIndex(notLongEnough) = [];

%DB stretched:
DBstretched = DB(stretchedIndex);
% 122 sectors are left.

%%
% Groups in DBstretched:
stretchedGroup = getGroup(DBstretched);
% All 3 groups are represented

%%
% Subgroups
stretchedSubgroup = getSubgroup(DBstretched);
% Only 32 out of the 70 subgroups are represented

%%
% Amino acid representation:
stretchedAACount = sum(aaSectorDBCount(DBstretched),2);
stretchedAAFreq = stretchedAACount / sum(stretchedAACount);

%%
figure('Name','Amino acid frequency in the stretched sectors')
subplot(2,1,1);
bar(stretchedAAFreq);
title('Amino acid frequency in the stretchted sectors')
ylabel('Frequency');
ylim([0 0.12]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);

subplot(2,1,2);
bar([(stretchedAAFreq-freqAAinDB)./freqAAinDB, (stretchedAAFreq-freqAAinMembraneProteins)./freqAAinMembraneProteins,(stretchedAAFreq-freqAAinUniprot)./freqAAinUniprot ]);
title('Relative amino acid frequency in the stretched sectors compared to different protein database')
ylabel('Relative enrichement');
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);
legend('relative enrichment from the sectors db','relative enrichment from the membrane protein db', 'relative enrichment from the uniprot db', 'Location', 'SouthEast');

% Compared to the whole sector db, the stretched sectors are enriched in
% Cisteine (20% more), Tryptophan (10%) and Tyrosine (20%) and are depleted
% in Arginine (15%), Histidine (22%) and Methionine (40%);

%%
% Residue PCA:
freqAAstretchedDB = aaSectorDBFreq(DBstretched);
[vect, val, eigenBase] = residuePCA(freqAAstretchedDB, 3, 1);

%%
h = clustergram(freqAAstretchedDB, 'Colormap', jet);

%%
% 3 lowest sectors on mode 20:
[~, perm20] = sort(eigenBase(end,:));
sector_lowest20 = perm20(end-2:end);

% 3 lowest sectors on mode 19:
[~, perm19] = sort(eigenBase(end-1,:));
sector_lowest19 = perm19(end-2:end);

% 3 highest sectors on mode 18:
[~, perm18] = sort(eigenBase(end-2,:));
sector_lowest18 = perm18(1:3);

% %%
% % Eigen vectors with the highest spread on C, H, M and Y
% % C should be enriched (aa number 5)
% [~, permC] = sortrows(abs(vect)', 5);
% modeC = permC(end);
% 
% % H should be depleted (aa number 9)
% [~, permH] = sortrows(abs(vect)', 9);
% modeH = permH(end);
% 
% % M should be depleted (aa number 13)
% [~, permM] = sortrows(abs(vect)', 13);
% modeM = permM(end);
% 
% % Y should be enriched (aa number 19)
% [~, permY] = sortrows(abs(vect)', 19);
% modeY = permY(end);
% 
% %%
% figure('Name', 'Scatter plot on the mode or higher enrichment/depletion of aa')
% scatter3(eigenBase(modeH,:), eigenBase(modeM,:) , eigenBase(modeY,:));
% xlabel('Histidine');
% ylabel('Methionine');
% zlabel('Tyrosine');


%% Exploring transmembrane
transmembraneSectorCell = importdata('transmembraneSectorCell.mat');
%%
% Percentage of the aa in sectors in the transmembrane region:
aaInTrans = 0;
aaInSector = 0;
percentInTransPerSector = zeros(1, numel(DB));
for i = 1:numel(transmembraneSectorCell)
    aaInTrans = sum(transmembraneSectorCell{i}) + aaInTrans;
    aaInSector = aaInSector + numel(transmembraneSectorCell{i});
    percentInTransPerSector(i) = sum(transmembraneSectorCell{i}) / numel(transmembraneSectorCell{i});
end

percentInTrans = aaInTrans / aaInSector;

%%
% Amino acid frequency in the transmembrane part of the sectors
seq = '';
for i = 1:numel(transmembraneSectorCell)
    sequence = DB{i}.Sequence;
    seq = [seq sequence(transmembraneSectorCell{i} == 1)];
end

aaFreqTransmembrane = aacount(seq);
aaFreqTransmembrane = struct2cell(aaFreqTransmembrane);
aaFreqTransmembrane = cell2mat(aaFreqTransmembrane);
aaFreqTransmembrane = aaFreqTransmembrane / sum(aaFreqTransmembrane);

%%
figure('Name', 'Percentage of amino acid in transmembrane region per sector');
hist(percentInTransPerSector, 20);
ylabel('Number of sectors', 'FontSize',18,...
       'FontWeight','bold');
xlabel('Percentage of transmembrane amino acids', 'FontSize',18,...
       'FontWeight','bold');
title('Distribution of the transmembrane residues in sectors', 'FontSize',20,...
       'FontWeight','bold');
set(gca, 'FontSize',18,...
       'FontWeight','bold');
%%
figure('Name','Amino acid frequency in the transmembrane aa in sectors')
subplot(2,1,1);
bar(stretchedAAFreq);
title('Amino acid frequency')
ylabel('Frequency');
ylim([0 0.12]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);

subplot(2,1,2);
bar([(aaFreqTransmembrane-freqAAinDB)./freqAAinDB, (aaFreqTransmembrane-freqAAinMembraneProteins)./freqAAinMembraneProteins,(aaFreqTransmembrane-freqAAinUniprot)./freqAAinUniprot ]);
title('Relative amino acid frequency in the stretched sectors compared to different protein database')
ylabel('Relative enrichement');
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);
legend('relative enrichment from the sectors db','relative enrichment from the membrane protein db', 'relative enrichment from the uniprot db', 'Location', 'SouthEast');



%%
% Let's isolate the sectors which have more than 80 % of their aa in the
% transmembrane region.

trans80 = find(percentInTransPerSector >= 0.8);

notLongEnough = [];
for i = 1:numel(trans80)
    if DB{trans80(i)}.Length < 5
        notLongEnough = [notLongEnough i];
    end
end

trans80(notLongEnough) = [];

DBtrans80 = DB(trans80);

%%e
figure('Name', 'EigenValues of the sectors with more the 80% in transmembrane regions')
subplot(2,1,1);

plot3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:),'.', ...
    'Color',[0.7 0.7 0.7], 'MarkerSize', 8);
hold on;
plot3(normEigenValueDB(3,trans80), ...
    normEigenValueDB(1,trans80), ...
    normEigenValueDB(2,trans80),'.', ...
    'MarkerSize', 15, 'Color', 'b');




grid on;
title('Transmembrane sectors')
xlabel('normalized eigenValue 3');
ylabel('normalized eigenValue 1');
zlabel('normalized eigenValue 2');

subplot(2,1,2);

plot3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:),'.', ...
    'Color',[0.7 0.7 0.7], 'MarkerSize', 8);
hold on;
plot3(eigenValueDB(3,trans80), ...
    eigenValueDB(1,trans80), ...
    eigenValueDB(2,trans80),'.', ...
    'MarkerSize', 15, 'Color', 'b','Markersize', 16);



grid on;
title('Transmembrane sectors')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');


%%
% Let's isolate the sectors which have more than 90 % of their aa in the
% transmembrane region.

trans90 = find(percentInTransPerSector >= 0.9);

notLongEnough = [];
for i = 1:numel(trans90)
    if DB{trans90(i)}.Length < 5
        notLongEnough = [notLongEnough i];
    end
end

trans90(notLongEnough) = [];

DBtrans90 = DB(trans90);

%%e
figure('Name', 'EigenValues of the sectors with more the 90% in transmembrane regions')
subplot(2,1,1);

% plot3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:),'.', ...
%     'Color',[0.7 0.7 0.7], 'MarkerSize', 8);
%hold on;
plot3(normEigenValueDB(3,trans90), ...
    normEigenValueDB(1,trans90), ...
    normEigenValueDB(2,trans90),'.', ...
    'MarkerSize', 15, 'Color', 'b');




grid on;
title('Transmembrane sectors')
xlabel('normalized eigenValue 3');
ylabel('normalized eigenValue 1');
zlabel('normalized eigenValue 2');

subplot(2,1,2);
% 
% plot3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:),'.', ...
%     'Color',[0.7 0.7 0.7], 'MarkerSize', 8);
%hold on;
plot3(eigenValueDB(3,trans90), ...
    eigenValueDB(1,trans90), ...
    eigenValueDB(2,trans90),'.', ...
    'MarkerSize', 15, 'Color', 'b');



grid on;
title('Transmembrane sectors')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');



%%
% Let's isolate the sectors which have more than 50 % of their aa in the
% transmembrane region.

trans50 = find(percentInTransPerSector >= 0.5);

notLongEnough = [];
for i = 1:numel(trans50)
    if DB{trans50(i)}.Length < 5
        notLongEnough = [notLongEnough i];
    end
end

trans50(notLongEnough) = [];

DBtrans50 = DB(trans50);

%%
% Let's isolate the sectors which have more than 75 % of their aa in the
% transmembrane region.

trans75 = find(percentInTransPerSector >= 0.75);

notLongEnough = [];
for i = 1:numel(trans75)
    if DB{trans75(i)}.Length < 5
        notLongEnough = [notLongEnough i];
    end
end

trans75(notLongEnough) = [];

DBtrans75 = DB(trans75);

%%
numTransSector = numel(DBtrans75);
numTransGroup = size(getGroup(DBtrans75),2);
numTransSubgroup = size(getSubgroup(DBtrans75),2);
numTransPdb = numel(getPdb(DBtrans75));


aaFreqDB75 = aaSectorDBCount(DBtrans75);
aaFreqDB75 = sum(aaFreqDB75,2)/ sum(sum(aaFreqDB75,2));
aaEnrichment = aaFreqDB75 ./ freqAAinMembraneProteins;

%%
figure('Name', 'EigenValues of the sectors with more than 75% aa in transmembrane regions')
subplot(2,1,1);

plot3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:),'.', ...
    'Color',[0.5 0.5 0.5], 'MarkerSize', 7);
hold on;
plot3(normEigenValueDB(3,trans75), ...
    normEigenValueDB(1,trans75), ...
    normEigenValueDB(2,trans75),'.', ...
    'MarkerSize', 18, 'Color', 'b');




grid on;
% title('Transmembrane sectors principal component projection', 'FontSize',20,...
%        'FontWeight','bold');
xlabel('normalized eigenValue 3', 'FontSize',11);
ylabel('normalized eigenValue 1', 'FontSize',11);
zlabel('normalized eigenValue 2', 'FontSize', 11);
 set(gca, 'FontSize',11);
 ylim([0.3 1])
 zlim([0 0.5])
 xlim([0 0.3])

   legend('All membrane sectors','Transmembrane sectors', 'FontSize',11);
subplot(2,1,2);
% 
plot3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:),'.', ...
    'Color',[0.5 0.5 0.5], 'MarkerSize', 7);
hold on;
plot3(eigenValueDB(3,trans75), ...
    eigenValueDB(1,trans75), ...
    eigenValueDB(2,trans75),'.', ...
    'MarkerSize', 18, 'Color', 'b');



grid on;
zlim([0 500]);
% title('Transmembrane sectors principal component projection', 'FontSize',20,...
%        'FontWeight','bold');
xlabel('eigenValue 3','FontSize',11);
ylabel('eigenValue 1', 'FontSize',11);
zlabel('eigenValue 2', 'FontSize',11);
set(gca, 'FontSize',11);
legend('All membrane sectors','Transmembrane sectors', 'FontSize',11);
xlim([0 250]);
ylim([0 5000]);
zlim([0 500]);
%    legend('All sectors','Transmembrane sectors', 'FontSize',18,...
%        'FontWeight','bold');

%%
% Clustering of the trans sector normalized eigenVal:
h = clustergram([normEigenValueDB(1,trans75); ...
    normEigenValueDB(2,trans75); ...
    normEigenValueDB(3,trans75)], 'Colormap', jet);

%%
% Clustering of the trans sector  eigenVal:
h = clustergram(zscore([eigenValueDB(1,trans75); ...
    eigenValueDB(2,trans75); ...
    eigenValueDB(3,trans75)]), 'Colormap', jet);

%%
% K mean of the normalized eigenVal:
opts = statset('Display','final');
X = [normEigenValueDB(1,trans75); ...
    normEigenValueDB(2,trans75); ...
    normEigenValueDB(3,trans75)]';

[idx,ctrs] = kmeans(X,5,...
              'Replicates',10000,'Options',opts);

          %%
figure('Name','K mean clustering of the normalized eigen value of trans75');
 plot3(X(idx==1,3),X(idx==1,1),X(idx==1,2),'b.','MarkerSize',18)
 hold on
 grid on
 plot3(X(idx==2,3),X(idx==2,1),X(idx==2,2),'g.','MarkerSize',18)
 plot3(X(idx==3,3),X(idx==3,1),X(idx==3,2),'y.','MarkerSize',18)
 plot3(X(idx==4,3),X(idx==4,1),X(idx==4,2),'m.','MarkerSize',18)
 plot3(X(idx==5,3),X(idx==5,1),X(idx==5,2),'r.','MarkerSize',18)
%  title('Scatter plot of the normalized eigen value of the transmembrane sectors');
 xlabel('Normalized eigenValue 3');
ylabel('Normalized eigenValue 1');
zlabel('Normalized eigenValue 2');
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5', 'FontSize',18,...
       'FontWeight','bold');
 
 %%
 Y = [eigenValueDB(1,trans75); ...
    eigenValueDB(2,trans75); ...
    eigenValueDB(3,trans75)]';

figure('Name','Clustering of the eigen value of trans75');
 plot3(Y(idx==1,1),Y(idx==1,2),Y(idx==1,3),'b.','MarkerSize',12)
 hold on
 grid on
 plot3(Y(idx==2,1),Y(idx==2,2),Y(idx==2,3),'g.','MarkerSize',12)
 plot3(Y(idx==3,1),Y(idx==3,2),Y(idx==3,3),'y.','MarkerSize',12)
 plot3(Y(idx==4,1),Y(idx==4,2),Y(idx==4,3),'m.','MarkerSize',12)
 plot3(Y(idx==5,1),Y(idx==5,2),Y(idx==5,3),'r.','MarkerSize',12)
 title('Scatter plot of the eigen values of the transmembrane sectors');
 xlabel('EigenValue 1');
ylabel('EigenValue 2');
zlabel('EigenValue 3');

%%
% PLotting by group
groups75 = getGroup(DBtrans75);
[~, groups75_1] = getSectorsByGroup(DBtrans75, groups75{1});
[~, groups75_2] = getSectorsByGroup(DBtrans75, groups75{2});


figure('Name','Clustering of the eigen value of trans75');
subplot(2,1,1)
 plot3(Y(groups75_1,1),Y(groups75_1,2),Y(groups75_1,3),'m.','MarkerSize',12)
 hold on
 grid on
 plot3(Y(groups75_2,1),Y(groups75_2,2),Y(groups75_2,3),'b.','MarkerSize',12)
 title('Scatter plot of the eigen values of the transmembrane sectors');
 xlabel('EigenValue 1');
ylabel('EigenValue 2');
zlabel('EigenValue 3');
legend('Alpha helical', 'Beta barrel');

subplot(2,1,2)
 plot3(X(groups75_1,1),X(groups75_1,2),X(groups75_1,3),'m.','MarkerSize',12)
 hold on
 grid on
 plot3(X(groups75_2,1),X(groups75_2,2),X(groups75_2,3),'b.','MarkerSize',12)
 title('Scatter plot of the normalized eigen values of the transmembrane sectors');
  xlabel('Normalized eigenValue 1');
ylabel('Normalized eigenValue 2');
zlabel('Normalized eigenValue 3');

%%
%PLotting by subgroup
subgroups75 = getSubgroup(DBtrans75);
subgroups75_indexes = {};
for i = 1:numel(subgroups75)
    [~, indexesTemp] = getSectorsBySubgroup(DBtrans75, subgroups75{i});
    subgroups75_indexes = [subgroups75_indexes {indexesTemp}];
end

figure
for i = 1:numel(subgroups75)
    subplot(5,5,i)
    plot3(X(subgroups75_indexes{i},1),X(subgroups75_indexes{i},2),X(subgroups75_indexes{i},3),'.','MarkerSize',12)
     title(subgroups75{i});
     grid on
     zlim([0 0.3]);
     xlim([.4 1]);
     ylim([0 0.5]);
end

%%
% Plotting by subgroup and cluster
subgroups75_indexes_cluster = cell(5, numel(subgroups75));
for j = 1:numel(subgroups75)
    [~, indexesTemp] = getSectorsBySubgroup(DBtrans75, subgroups75{j});
    for i = 1:5
        indexesTempi = intersect(idx(idx == i), indexesTemp);
        subgroups75_indexes_cluster{i,j} = indexesTempi ;
    end
end
c = [1 0 0 ; 0 1 0 ; 0 0 1 ; 0.5 0.5 0; 0.5 0 0.5];
figure
for i = 1:numel(subgroups75)
    subplot(5,5,i)
    for j = 1:5
        if ~ isempty(subgroups75_indexes_cluster{j,i})
        plot3(X(subgroups75_indexes_cluster{j,i},1),X(subgroups75_indexes_cluster{j,i},2),X(subgroups75_indexes_cluster{j,i},3),'.','Color', c(j,:),'MarkerSize',12)
        hold on;
        end
    end
     grid on
     zlim([0 0.3]);
     xlim([.4 1]);
     ylim([0 0.5]);
     hold off;
     title(subgroups75{i});
end

%%
DBtrans75_1 = DBtrans75((idx==1));
DBtrans75_2 = DBtrans75((idx==2));
DBtrans75_3 = DBtrans75((idx==3));
DBtrans75_4 = DBtrans75((idx==4));
DBtrans75_5 = DBtrans75((idx==5));

save('DBtrans75_1.mat', 'DBtrans75_1');
save('DBtrans75_2.mat', 'DBtrans75_2');
save('DBtrans75_3.mat', 'DBtrans75_3');
save('DBtrans75_4.mat', 'DBtrans75_4');
save('DBtrans75_5.mat', 'DBtrans75_5');
