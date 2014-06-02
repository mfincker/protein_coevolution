%% Load Data 
%%
% Get Cluster Data
%[TP53_clusters, TP53_extra] = blast2clust('BAC16799');
load('clusters.mat');
%
%[TP53_clusters, TP53_extra] = msa2clust(msa);
TP53_clusters = clusters;

%% Visualize MSA
imagesc(msa);
title('MSA of Homo Sapien p53 Sequences', 'FontSize', 20);
ylabel('No. of Sequences', 'FontSize', 18);
xlabel('Residue Index', 'FontSize', 18);
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
            if ~ischar(somatic_mutation_raw{i,7}) && somatic_mutation_raw{i,7} <= 393
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
mutations_and_what_sectors_theyarein = find(cell2mat(somatic_mutation_trim(2:end,17))>0);
somatic_mutations_only_in_sectors = somatic_mutation_trim(mutations_and_what_sectors_theyarein + 1,:);

%% Relative Sector Enrichment
% For each sector
figure
Sector_enrichment = zeros(length(TP53_clusters),1);
for j = 1:length(TP53_clusters)
Sector_enrichment(j,1) = (inSectorCount(1,j)./length(somatic_mutation_trim(:,2)))./(length(cell2mat(TP53_clusters(1,j)))./393);
end
bar(1:length(Sector_enrichment),Sector_enrichment,'b');
title('Mutation enrichment per sector', 'FontSize', 20);
xlabel('Sector number', 'FontSize', 18);
ylabel('Enrichment (% of (mutants in sector/totalmutants)./(residues in sector/totalresidues)');

%% Compensatory mutations
%
% First we will look at each mutation at each residue within a sector
% across species to find "fixed" mutations. In order to do this we must
% make a map of all the residues in p53 to compare all of the other
% sequences to.
%%
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
% Extract the mutations in the clusters. For somatic_mutation_trim, the
% columns are
%%
% 2: indices of reference sequence
%%
% 6: wildtype AA
%%
% 7: mutant AA
%%
% 17: the index of the cluster the mutation is in
%%
% We also have to convert the three letter amino acid codes to the integer
% representation in MATLAB..
clust_mut_map = find(cell2mat(somatic_mutation_trim(:,17)));
clust_mutation_cell = somatic_mutation_trim(clust_mut_map,[2 6 7 17]);

mut_count = 1;
for i=1:length(clust_mutation_cell)
    % only consider deletion, substitution and premature stops, since there
    % are mutation from 1 AA to 2 AA or other peculiar ones in the database
    if (length(clust_mutation_cell{i,2}) <=4 && length(clust_mutation_cell{i,3}) <=4)
        clust_mutation(mut_count,1) = cell2mat(clust_mutation_cell(i,1));
        clust_mutation(mut_count,4) = cell2mat(clust_mutation_cell(i,4));
        for aa = 2:3 % WT AA and MT AA
            if (strcmp(clust_mutation_cell(i,aa), 'NA'))
                 % deletion causing gaps
                clust_mutation(mut_count,aa) = aa2int('-');
            elseif (strcmp(clust_mutation_cell(i,aa), 'STOP'))
                 % translation stop, * in MATLAB
                clust_mutation(mut_count,aa) = aa2int('*');
            else % substitution
                AminoAcid = aminolookup('Abbreviation', clust_mutation_cell{i,aa});
                oneLetterCode = AminoAcid(1,1);
                clust_mutation(mut_count,aa) = aa2int(AminoAcid(1,1));
            end
        end
        mut_count = mut_count + 1;
    end
end
clust_mutation = unique(clust_mutation, 'rows');



%% Create trimmed database of unique mutations for each residue
unique_mutations_cell = somatic_mutation_trim(:,[2,6,7,17]);
mut_count = 1;
for i = 1:length(unique_mutations_cell(:,1))
    for j = 1:length(unique_mutations(1,:))
        if isa(unique_mutations{i,j},'double')
            unique_mutations{i,j} = num2str(unique_mutations{i,j});
        end
        if (length(unique_mutations_cell{i,2}) <=4 && length(unique_mutations_cell{i,3}) <=4)
        unique_mutations(mut_count,1) = cell2mat(unique_mutations_cell(i,1));
        unique_mutations(mut_count,4) = cell2mat(unique_mutations_cell(i,4));
            for aa = 2:3 % WT AA and MT AA
                if (strcmp(unique_mutations_cell(i,aa), 'NA'))
                     % deletion causing gaps
                    unique_mutations(mut_count,aa) = aa2int('-');
                elseif (strcmp(unique_mutations_cell(i,aa), 'STOP'))
                     % translation stop, * in MATLAB
                    unique_mutations(mut_count,aa) = aa2int('*');
                else % substitution
                    AminoAcid = aminolookup('Abbreviation', unique_mutations_cell{i,aa});
                    oneLetterCode = AminoAcid(1,1);
                    unique_mutations(mut_count,aa) = aa2int(AminoAcid(1,1));
                end
            end
        mut_count = mut_count + 1;
        end
        unique_mutations(i,j) = cellstr(unique_mutations(i,j));
    end
end
unique_mutations_new = unique(unique_mutations,'rows');




% figure
% New_Sector_enrichment = zeros(length(TP53_clusters),1);
% for j = 1:length(TP53_clusters)
% New_Sector_enrichment(j,1) = (length(find(clust_mutation(:,4)==j))./length(clust_mutation(:,4)))./(length(cell2mat(TP53_clusters(1,j)))./393);
% end

%%
% The actual part of finding compensatory mutation. First, we look at each
% mutation in the mutation database, and count the frequencies of each
% carcinogenic mutation in species other than human beings. The fifth
% column ot clust_mutation is added in this sections, and it contains the
% sequences in clustMSA with the given mutation.
clust_mutation_cell = num2cell(clust_mutation);
%%

for i=1:length(clust_mutation)
    mut_index = clust_mutation(i,1); % index of mutation residue
    mut_aa = clust_mutation(i,3); % mutant AA
    clust_index = clust_mutation(i,4); % index of cluster the mutation is in
    clust_msa = clustMSA(:,:,clust_index);
    % index of the current mutation in the cluster's MSA
    msa_mut_index = find(clust_msa(1,:) == mut_index); 
    % the residues in each sequence at the mutation index
    mut_msa = clust_msa(:,msa_mut_index);
    % find the sequences with the current mutation 
    mut_msa_map = find(mut_msa(:,1) == mut_aa);
    mut_seqs = clust_msa(mut_msa_map,:);
    % store the information as cells
    clust_mutation_cell{i,5} = length(mut_seqs(:,1));
    clust_mutation_cell{i,6} = mut_seqs;
end

%%
% Plot the number of carcinogenic mutations present in other species in the
% MSA.
for i=1:length(TP53_clusters)
    subplot(2,3,i);
    clust_mut_i = find(cell2mat(clust_mutation_cell(:,4)) == i);
    x_mutIndex = cell2mat(clust_mutation_cell(clust_mut_i,1));
    y_mutCount = cell2mat(clust_mutation_cell(clust_mut_i,5));
    scatter(x_mutIndex, y_mutCount, 35, 'filled');
    xlabel('Index of Carcinogenic Mutation', 'FontSize', 18);
    titleStr = sprintf('Sector %d', i);
    title(titleStr, 'FontSize', 20);
    ylabel('No. of Seqs with Mutation', 'FontSize', 18);
end

%%
% From the above graph, we can see that there is one mutation in cluster 3
% that exists in most of the sequences in the MSA. We extract the MSA of
% the sequences with this mutation, and visualize it using imagesc.
max_mut_count_index = find(cell2mat(clust_mutation_cell(:,5))==max(cell2mat(clust_mutation_cell(:,5))));
max_mut_clust_index = cell2mat(clust_mutation_cell(max_mut_count_index,4));
max_mut_msa = vertcat(clustMSA(:,:,max_mut_clust_index));
% set up MSA so that top 1/3 of the image would be the reference sequence
max_mut_msa = max_mut_msa(2:length(max_mut_msa),1:length(cell2mat(TP53_clusters(1,max_mut_clust_index))));
max_mut_msa = vertcat(repmat(max_mut_msa(1,:),length(max_mut_msa)/2,1),max_mut_msa);
imagesc(max_mut_msa);
xlabel('Residues in Cluster 3', 'FontSize', 18);
ylabel('Sequences with Mutation and Reference Sequence', 'FontSize', 18);
title('Cluster MSA of Sequences', 'FontSize', 20);
xticks = linspace(1,length(max_mut_msa(1,:)), length(max_mut_msa(1,:)));
set(gca, 'XTick', xticks, 'XTickLabel', xticks);
set(gca, 'YTick', [], 'YTickLabel', []);

%% Phenotypic enrichment for each sector
%This section is meant to find the number of different phenotypes expressed
%in each sector, and the relative enrichment of each phenotype within each
%sector.







%% Map to pdb
%pdb is 3Q05
