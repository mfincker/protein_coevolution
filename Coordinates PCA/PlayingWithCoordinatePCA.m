%% Coordinate PCA trials

% Keep track of different things I tried on the CoordinateS PCA.

% Load database
filename = './prokaryote_Sectors/prokaryoticSectorDB.mat';
sectorDatabase = importdata(filename);

% Extract membrane sectors indexes:
membraneGroup = zeros(numel(sectorDatabase),1);
for i = 1:numel(sectorDatabase)
	if sectorDatabase{i}.Membrane == 1
		membraneGroup(i) = i;
	end 
end
membraneGroup = nonzeros(membraneGroup); 

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
% Normalizing the database
normEigenValueDB = zscore(eigenValueDB')';

%%
% Scatter plot of the normalized eigenValues:
figure
scatter3(normEigenValueDB(3,:),normEigenValueDB(1,:),normEigenValueDB(2,:));
title('Scatter plot of the normalized eigen values of each sector in the database')
xlabel('normalized eigenValue 3');
ylabel('normalized eigenValue 1');
zlabel('normalized eigenValue 2');

%%
% Histograms of the distribution
figure
subplot(2,3,1); hist(normEigenValueDB(1,:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r');
title('Histogram of normalized eigenValue 1 in the sectors');
xlabel('Value'); ylabel('Count');

subplot(2,3,2); hist(normEigenValueDB(2,:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
title('Histogram of normalized eigenValue 2 in the sectors');
xlabel('Value'); ylabel('Count');

subplot(2,3,3); hist(normEigenValueDB(3,:));
title('Histogram of normalized eigenValue 3 in the sectors');
xlabel('Value'); ylabel('Count');

subplot(2,3,4); hist(normEigenValueDB(1,:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r');
hold on
hist(normEigenValueDB(2,:));
title('Overlayed histograms of normalized eigenValue 1 and 2');
xlabel('Value'); ylabel('Count');

subplot(2,3,5); hist(normEigenValueDB(1,:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r');
hold on
hist(normEigenValueDB(3,:));
title('Overlayed histograms of normalized eigenValue 1 and 3');
xlabel('Value'); ylabel('Count');

subplot(2,3,6); hist(normEigenValueDB(2,:));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b');
hold on
hist(normEigenValueDB(3,:));
title('Overlayed histograms of normalized eigenValue 2 and 3');
xlabel('Value'); ylabel('Count');


%%
% Cut off to create groups:
% Group 1 : sectors where eig1~ >= 90%
% Group 2 : sectors not in group 1 and where eig1~ + eig2~ >= 90%
% Group 3 : the rest

sumEigenValueDB = sum(eigenValueDB,1);
sumEigenValueDB = repmat(sumEigenValueDB,3,1);
eigenValueDB_contribution = eigenValueDB ./ sumEigenValueDB;

eigenValueDB_cumulContribution = [eigenValueDB_contribution(1,:) ; ...
			eigenValueDB_contribution(1,:)+eigenValueDB_contribution(2,:) ; ...
			sum(eigenValueDB_contribution,1)];

group1 = find(eigenValueDB_contribution(1,:) >= 0.9);

group2 = find(eigenValueDB_cumulContribution(2,:) >= 0.9) ;
group2 = setdiff(group2, group1);

group3 = setdiff([1:numel(eigenValueDB(1,:))],union(group2,group1));

sectorDatabase_group1 = {sectorDatabase{group1}};
sectorDatabase_group2 = {sectorDatabase{group2}};
sectorDatabase_group3 = {sectorDatabase{group3}};

%%
% Displaying the color coded scatter plot
figure
c = zeros(size(eigenValueDB,2),3);
c(group1,:) = repmat([0 0.75 0.75],numel(group1),1);
c(group2,:) = repmat([0.75 0 0.75],numel(group2),1);
c(group3,:) = repmat([0.75 0.75 0],numel(group3),1);
s = scatter3(eigenValueDB(3,:),eigenValueDB(1,:), ...
				eigenValueDB(2,:),40, c,'MarkerFaceColor',[0.4 0.4 0.4]);

title('Scatter plot of the eigen values of each sector in the database')
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');


%%
% Membrane proteins

% Extract membrane sectors indexes:
membraneGroup = zeros(numel(sectorDatabase),1);
for i = 1:numel(sectorDatabase)
	if sectorDatabase{i}.Membrane == 1
		membraneGroup(i) = i;
	end 
end
membraneGroup = nonzeros(membraneGroup); 

sectorDatabase_membrane = {sectorDatabase{membraneGroup}};
figure
c = repmat([0.5 0.5 0.5],size(eigenValueDB,2),1);
c(membraneGroup,:) = repmat([0.75 0 0.75],numel(membraneGroup),1);
s = scatter3(eigenValueDB(3,:),eigenValueDB(1,:), ...
				eigenValueDB(2,:),40, c);

title(['Scatter plot of the eigen values of each sector in the database' ...
 			char(10) 'Membrane protein are in pink.']);
xlabel('eigenValue 3');
ylabel('eigenValue 1');
zlabel('eigenValue 2');


























