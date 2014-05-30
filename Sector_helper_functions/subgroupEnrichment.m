function [ enrichment,medianenrichment ] = subgroupEnrichment( sectorDB )
% This function will calculate the enrichement over an entire subgroup.  It
% will average the enrichment over all the proteins in the subgroup which
% is the output of the function proteinEnrichment
%   Detailed explanation goes here

% separate subgroup sectors by protein
sectorbyProt = separateSectorsByProtien(sectorDB);

% make 20xn matrix where n is the number of protiens in the subgroup
enrichment = zeros(20, length(sectorbyProt));
allrescount = zeros(1, length(sectorbyProt));

for i= 1:length(sectorbyProt)
    [currentenrichment, numresidues] = proteinEnrichment(sectorbyProt{i});
    allrescount(i) = numresidues;
    enrichment(:,i) = currentenrichment;
    
end

enrichment(:, allrescount<0) = [ ];

% avg enrichment per subgroup
medianenrichment = median(enrichment,2);
% plot


end

