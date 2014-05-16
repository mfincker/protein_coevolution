function [ eigValuesDB ] = sectorDBcoordPCA( sectorDatabase, showplots)
%This function will be able to perform PCA over all sectors in the database
    %pca_loop will return the pca results, and plot and compare differnet
    %pca results for all sectors. 

	% to keep track of the  lambdas for each sector.
	eigValuesDB = zeros(3, numel(sectorDatabase)) ;

	for i = 1:numel(sectorDatabase)
	    eigValuesDB( :, i) = coordinatePCA( sectorDatabase{i}, ...
	    									showplots );
	end

	% visualize the data in a meaningful way.

% >> eigValuesDB = sectorDBcoordPCA(sectorDatabase, 0);
% >> zScore = zscore(eigValuesDB')';
% >> figure, scatter3(zScore(3,:),zScore(1,:),zScore(2,:))

end

