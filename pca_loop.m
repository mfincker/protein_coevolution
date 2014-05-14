function [ eigenvector_sector_database ] = pca_loop( sectorDatabase, showplots)
%This function will be able to perform PCA over all sectors in the database
    %pca_loop will return the pca results, and plot and compare differnet
    %pca results for all sectors. 
    

% to keep track of the  lambdas for each sector.

eigenvector_sector_database = cell(1,length(sectorDatabase));
% get eigenvectors from sector_PCA on all sectors

    for i = 1:numel(sectorDatabase)
    [sect_eigvectors] = sector_pca(sectorDatabase{1,i}, showplots )
        for i = 1:length(sectorDatabase)
            eigenvector_sector_database{:,i} = sect_eigvectors
            
        end
    end

% visualize the data in a meaningful way.

end

