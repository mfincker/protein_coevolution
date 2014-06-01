%% Load Data 
%%
% Get Cluster Data
%[TP53_clusters, TP53_extra] = blast2clust('BAC16799');
load('clusters.mat');
%
%[TP53_clusters, TP53_extra] = msa2clust(msa);
TP53_clusters = clusters;
%%
% Read Functional Assessment Data
function_assess_raw = importdata('functionalAssessmentIARC TP53 Database, R17.txt');
function_assess_raw = function_assess_raw.textdata;
%%
% Read Somatic Mutation Data
[num, txt, somatic_mutation_raw] = xlsread('somatic.xlsx');

%% Trim Neutral Mutations
% Column 20 is the WT_AA, column 21 is the Mutant_AA.
mut_count = 1;
for i=1:length(somatic_mutation_raw)
    % nucleotide mutation without AA change
    if ~strcmp(somatic_mutation_raw(i,20),somatic_mutation_raw(i,21)) 
        
        if ~strcmp(somatic_mutation_raw(i,20),'NA') || ~strcmp(somatic_mutation_raw(i,21),'NA')
            if ~ischar(somatic_mutation_raw{i,7}) && somatic_mutation_raw{i,7} <= 391
            somatic_mutation_trim(mut_count,:) = somatic_mutation_raw(i,:);
            mut_count = mut_count + 1;
            end
        end
    end
end

%% Trim Columns not going to be used
% Column 1-4, 6, 8-16, 22, 27-30, 36-65
somatic_mutation_trim = somatic_mutation_trim(:,[5 7 17:21 23:26 31:35 65]);
 
%% Overall Sector Enrichment

inSectorCount= zeros(1,7);
outSectorCount = 0;

for i=1:length(somatic_mutation_trim)
    %I changed the somatic index accession from 7 to 2 to reflect our new trimmed
    %matrix. This is because the mutations are in a new column
    cluster_index = is_in_sector(TP53_clusters, cell2mat(somatic_mutation_trim(i,2))); 
    %changed clusters to TP53_clusters
    if cluster_index > 0 %added this logical
        mutation_in_sector(inSectorCount(1,cluster_index) + 1,:,cluster_index) = somatic_mutation_trim(i,:);
        %added a multidimensional component to mutation_in_sector to allow
        %easy access to different sectors.
        inSectorCount(1,cluster_index) = inSectorCount(1,cluster_index) + 1;
    else
        mutation_out_sector(outSectorCount + 1,:) = somatic_mutation_trim(i,:);
        outSectorCount = outSectorCount + 1;
    end
    somatic_mutation_trim{i, 17} = cluster_index;
end

%% Relative Sector Enrichment
% For each sector
Sector_enrichment = zeros(length(TP53_clusters),1);
for j = 1:length(TP53_clusters)
Sector_enrichment(j,1) = (inSectorCount(1,j)./length(somatic_mutation_trim(:,2)))./(length(cell2mat(TP53_clusters(1,j)))./391);
end
bar(1:length(Sector_enrichment),Sector_enrichment,'r');
title('Mutation enrichment per sector', 'FontSize', 20);
xlabel('Sector number', 'FontSize', 18);
ylabel('Enrichment (% of (mutants in sector/totalmutants)./(residues in sector/totalresidues)');

%% Compensatory mutations
% Extract the residues in the MSA corresponding to the reference sequence.
map = find(msa(1,:) ~= 25);
seqmsa = msa(:,map);

%%
% Find the maximum size of cluster to facilitate setting up the MSA of
% clusters.
maxResidue_clust = 0;
for i=1:length(TP53_clusters)
    if (length(cell2mat(TP53_clusters(1,i))) > maxResidue_clust)
        maxResidue_clust = length(cell2mat(TP53_clusters(1,i)));
    end
end
%%
% Store the msa of the reference sequence into each cluster. The data is
% stored in a 3D matrix, with each layer containing the indices of the 
% clustered residues in the reference sequence and the MSA corresponding to
% those indices.
clustMSA = zeros(length(seqmsa)+1, maxResidue_clust, length(TP53_clusters));
for i=1:length(TP53_clusters)
    clust_index = cell2mat(TP53_clusters(1,i)); % indices of ref sequence
    clustMSA(1,1:length(clust_index),i) = clust_index;
    clustMSA(2:length(clustMSA),1:length(clust_index),i) = seqmsa(:,clust_index);
end

%%
%% Map to pdb
%pdb is 3Q05
