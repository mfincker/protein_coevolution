%% Load Data
%%
% Get Cluster Data
[TP53_clusters, TP53_extra] = blast2clust('BAC16799');

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
            if somatic_mutation_raw(i,7) <= 391
            somatic_mutation_trim(mut_count,:) = somatic_mutation_raw(i,:);
            mut_count = mut_count + 1;
        end
    end
end

%% Trim Columns not going to be used
% Column 1-4, 6, 8-16, 22, 27-30, 36-65
somatic_mutation_trim = somatic_mutation_trim(:,[5 7 17:21 23:26 31:35]);
%% Overall Sector Enrichment

inSectorCount = 0;
outSectorCount = 0;

for i=1:length(somatic_mutation_trim)
    %I changed the somatic index accession from 7 to 2 to reflect our new trimmed
    %matrix. This is because the mutations are in a new column
    cluster_index = is_in_sector(TP53_clusters, cell2mat(somatic_mutation_trim(i,2))); 
    %changed clusters to TP53_clusters
    if cluster_index > 0 %added this logical
        mutation_in_sector(inSectorCount + 1,:,cluster_index) = somatic_mutation_trim(i,:);
        %added a multidimensional component to mutation_in_sector to allow
        %easy access to different sectors.
        inSectorCount(1,cluster_index) = inSectorCount(1,cluster_index) + 1;
    else
        mutation_out_sector(outSectorCount + 1,:) = somatic_mutation_trim(i,:);
        outSectorCount = outSectorCount + 1;
    end
end

%% Relative Sector Enrichment
% For each sector
Sector_enrichment = zeros(length(TP53_clusters),1);
for j = 1:length(TP53_clusters)
Sector_enrichment(j) = (inSectorCount(1,j)./length(somatic_mutation_trim(:,2)))./(length(cell2mat(TP53_clusters(1,i)))./391);
end

%% Compensatory mutations
% 

