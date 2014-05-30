function [sectorIndex] = is_in_sector(clusters, index)

%% Is the Index in Sectors?
% This function determines whether an index is within the given sectors.
% 
% Parameters
%       clusters: a cell array of double arrays containing the clusters


sectorIndex = 0;
for i=1:length(clusters)
    cluster = cell2mat(clusters(1,i));
    if ismember(index, cluster)
        sectorIndex = i;
        break
    end
end