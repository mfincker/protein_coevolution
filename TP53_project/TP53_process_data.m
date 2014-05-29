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
            somatic_mutation_trim(mut_count,:) = somatic_mutation_raw(i,:);
            mut_count = mut_count + 1;
        end
    end
end

%% Trim Columns not going to be used
% Column 1-4, 6, 8-16, 22, 27-30, 36-65

%% Overall Sector Enrichment

inSectorCount = 0;
outSectorCount = 0;

for i=1:length(somatic_mutation_trim)
    cluster_index = is_in_sector(clusters, cell2mat(somatic_mutation_trim(i,7)));
    if cluster_index
        mutation_in_sector(inSectorCount + 1,:) = somatic_mutation_trim(i,:);
        
        inSectorCount = inSectorCount + 1;
    else
        mutation_out_sector(outSectorCount + 1,:) = somatic_mutation_trim(i,:);
        outSectorCount = outSectorCount + 1;
    end
end

%% Relative Sector Enrichment
% For each sector

%% Compensatory mutations
% 