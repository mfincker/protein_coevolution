%% Exploring the membrane sector DB

clc; close all;

% Importing the database
DB = importdata('../MembraneSectorDB/membraneSectorDB_053014_forCoord.mat'); 

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

%%
figure('Name','Distribution of the length of sectors');
subplot(1,2,1);
boxplot(lengths);
title('Distribution of the length of sectors');
ylabel('Number of reisdue');
set(gca,'XTickLabel','Membrane sector database');
subplot(1,2,2);
hist(lengths,25);
title('Distribution of the length of sectors');
xlabel('Nuumber of residue');
ylabel('Number of sectors');

%%
% Aminoc acid frequency in the sector DB:
numResidue = sum(lengths);
countAAinDB = aaSectorDBCount(DB);
countAAinDB = sum(countAAinDB, 2);
freqAAinDB = countAAinDB / numResidue;

% %%
% % Amino acid frequency in the membrane proteins used to make the DB:
countAAinMembraneProteins = zeros(20,numPdb);
for prot = 1:numPdb
    seq = getgenpept(pdb{prot});
    seq = seq.Sequence;
    count = aacount(seq);
	count = struct2cell(count);
	count = cell2mat(count);
    countAAinMembraneProteins(:, prot) = count;
end

countAAinMembraneProteins = sum(countAAinMembraneProteins,2);
save('countAAinMembraneProteins.mat', 'countAAinMembraneProteins');

% countAAinMembraneProteins = importdata('countAAinMembraneProteins.mat');
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
subgroupIndex = {};
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
eTransport = [1 40];

% ATPase
ATPase = [8, 48, 49];

% Oxygenase
oxygenase = [17, 52, 32];
cyclooxygenase = 32;

% Dehydrogenase
dehydrogenase = [30, 45];

% Oxydoreductase
oxydoreductase = [31, 50];

%%
c = [0.75 0.25 0.25 ; 0.25 0.75 0.25 ; 0.25 0.25 0.75];
figure('Name', 'EigenValues of the subgroups of interest')
subplot(2,5,1);
for i =1:numel(eTransport)
    plot3(eigenValueDB(3,subgroupIndex{eTransport(i)}), ...
        eigenValueDB(1,subgroupIndex{eTransport(i)}), ...
        eigenValueDB(2,subgroupIndex{eTransport(i)}),'.', ...
        'Color',c(i,:), 'MarkerSize', 12);
    hold on;
end
plot3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:),'.', ...
    'Color',[0.3 0.3 0.3], 'MarkerSize', 4);
ylim([0 20000]);
xlim([0 200]);
zlim([0 200]);
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
        'Color',c(i,:), 'MarkerSize', 12);
    hold on;
end
plot3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:),'.', ...
    'Color',[0.3 0.3 0.3], 'MarkerSize', 4);
ylim([0 0.5]);
xlim([0 0.5]);
zlim([0 0.5]);
grid on;
title('Electron transport sector')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

