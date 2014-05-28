%% Coordinate PCA trials - 27/05/2014

% Keep track of different things I tried on the CoordinateS PCA.

% % Load database
% filename = './MembraneProteins/MembraneSectorDB/membraneSectorDB_052714.mat';
% sectorDatabase = importdata(filename);

% % Extract membrane sectors indexes:
% membraneGroup = zeros(numel(sectorDatabase),1);
% for i = 1:numel(sectorDatabase)
% 	if sectorDatabase{i}.Membrane == 1
% 		membraneGroup(i) = i;
% 	end 
% end
% membraneGroup = nonzeros(membraneGroup); 
sectorDatabase = membraneSectorDB;
groups = getGroup(sectorDatabase);
subgroups = getSubgroup(sectorDatabase);

[group1, index1] = getSectorsByGroup(sectorDatabase, groups{1});
[group2, index2] = getSectorsByGroup(sectorDatabase, groups{2

subgroupSectors = {};
subgroupIndex = {};
for i = 1:numel(subgroups)
    [subgroupSector, subgroupI] = getSectorsBySubgroup(sectorDatabase, subgroups{i});
    subgroupSectors = [subgroupSectors subgroupSector];
    subgroupIndex = [subgroupIndex {subgroupI}];
end
        
    
% Run PCA, returns a 3 x num sector matrix of eigenValues 1, 2 and 3
eigenValueDB = sectorDBcoordPCA(sectorDatabase, 0);

%%
% Scatter plot of the eigenValues:
figure
scatter3(eigenValueDB(3,:),eigenValueDB(1,:),eigenValueDB(2,:));
title('Scatter plot of the eigen values of each sector in the database')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');

% Can't see much, except that we do have some spread.

%% Eigenvalue distribution and stats

% EigenValue 1 varies along a much wider range than
% eigenValue 2 or 3.
figure
boxplot(eigenValueDB','labels',{'eigenValue 1', ...
					'eigenValue 2','eigenValue 3'});

%%
% Stats on eigenvalues
eigMean = mean (eigenValueDB,2);
eigStd = std(eigenValueDB')';
eigMedian = median(eigenValueDB')';
disp({'Mean', 'Std', 'Median'});
disp([eigMean eigStd eigMedian]);

%%
% Histograms of the distribution
figure
subplot(2,3,1); hist(eigenValueDB(1,:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r');
title('Histogram of eigenValue 1 in the sectors');
xlabel('Value'); ylabel('Count');

subplot(2,3,2); hist(eigenValueDB(2,:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
title('Histogram of eigenValue 2 in the sectors');
xlabel('Value'); ylabel('Count');

subplot(2,3,3); hist(eigenValueDB(3,:));
title('Histogram of eigenValue 3 in the sectors');
xlabel('Value'); ylabel('Count');

subplot(2,3,4); hist(eigenValueDB(1,:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r');
hold on
hist(eigenValueDB(2,:));
title('Overlayed histograms of eigenValue 1 and 2');
xlabel('Value'); ylabel('Count');

subplot(2,3,5); hist(eigenValueDB(1,:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r');
hold on
hist(eigenValueDB(3,:));
title('Overlayed histograms of eigenValue 1 and 3');
xlabel('Value'); ylabel('Count');

subplot(2,3,6); hist(eigenValueDB(2,:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
hold on
hist(eigenValueDB(3,:));
title('Overlayed histograms of eigenValue 2 and 3');
xlabel('Value'); ylabel('Count');



%%
% Cut off to create groups:
% Group 1 : sectors where eig1~ >= 90%
% Group 2 : sectors not in group 1 and where eig1~ + eig2~ >= 90%
% Group 3 : the rest

sumEigenValueDB = sum(eigenValueDB,1);
sumEigenValueDB = repmat(sumEigenValueDB,3,1);
eigenValueDB_contributionRaw = eigenValueDB ./ sumEigenValueDB;

eigenValueDB_cumulContribution = [eigenValueDB_contributionRaw(1,:) ; ...
			eigenValueDB_contributionRaw(1,:)+eigenValueDB_contributionRaw(2,:) ; ...
			sum(eigenValueDB_contributionRaw,1)];

group1 = find(eigenValueDB_contributionRaw(1,:) >= 0.9);

group2 = find(eigenValueDB_cumulContribution(2,:) >= 0.9) ;
group2 = setdiff(group2, group1);

group3 = setdiff([1:numel(eigenValueDB(1,:))],union(group2,group1));

sectorDatabase_group1 = {sectorDatabase{group1}};
sectorDatabase_group2 = {sectorDatabase{group2}};
sectorDatabase_group3 = {sectorDatabase{group3}};


%%
% Normalized eigenValues
figure
scatter3(eigenValueDB_contributionRaw(3,:),eigenValueDB_contributionRaw(1,:),eigenValueDB_contributionRaw(2,:));
title('Scatter plot of the eigen values of each sector in the database')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');


figure
subplot(1,3,1);
plot(eigenValueDB_contributionRaw(1,:),eigenValueDB_contributionRaw(2,:),'g*');
hold on;
plot(eigenValueDB_contributionRaw(1,index2),eigenValueDB_contributionRaw(2,index2),'r*');
plot(eigenValueDB_contributionRaw(1,index1),eigenValueDB_contributionRaw(2,index1),'b*');
subplot(1,3,2);
plot(eigenValueDB_contributionRaw(1,:),eigenValueDB_contributionRaw(3,:),'*');
subplot(1,3,3);
plot(eigenValueDB_contributionRaw(2,:),eigenValueDB_contributionRaw(3,:),'*');





























