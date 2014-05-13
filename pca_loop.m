function [ PCA_loop_results ] = pca_loop( sectorDatabase, showplots)
%This function will be able to perform PCA over all sectors in the database
    %pca_loop will return the pca results, and plot and compare differnet
    %pca results for all sectors. 

% to keep track of the  lambdas for each sector.
PCAresults = zeros(numel(sectorDatabase), 3) 
size(PCAresults)

for i = 1:numel(sectorDatabase)
    PCAresults( i, :) = sector_pca( sectorDatabase{1, i}.Coordinates, 0 )
end

% visualize the data in a meaningful way.

end

