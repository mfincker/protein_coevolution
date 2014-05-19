%% Script to make a database of sectors from a
% text file containing pdbId, group, subgroup, memberProtein and relatedProtein
% on 5 different lines.


membraneSectorDB_Maeva = {};
PDBnotWorking_Maeva = struct([]);
PDBworking_Maeva = {};

file = fopen('./allMembraneProteinsDB_Maeva.txt');

pdbId = fgetl(file);
while ischar(pdbId)
    group = fgetl(file);
    subgroup = fgetl(file);

    memberProtein = fgetl(file);
    if isempty(memberProtein)
        memberProteins = regexp(memberProtein, '\t', 'split');
    else
        memberProteins = {};
    end

    relatedProtein = fgetl(file);
    if isempty(relatedProtein)
        relatedProteins = regexp(relatedProtein, '\t', 'split');
    else
        relatedProteins = {};
    end

    % Generate msa and clusters
    % try
        [clusters, data] = blast2clust_maeva(pdbId);
        disp(['For' pdbId ': found ', num2str(length(clusters)), ' clusters!']);

        for j = 1:numel(clusters)
            sector = Sector(data, clusters{j});
            sector.Group = group;
            sector.Subgroup = subgroup;
            sector.MemberProteins = memberProteins;
            sector.RelatedProteins = relatedProteins;
            membraneSectorDB_Maeva = [membraneSectorDB_Maeva {Sector(data, clusters{j})}];
        end
        PDBworking_Maeva = [PDBworking_Maeva; {pdbId}];

    catch err
        PDBnotWorking_Maeva = [PDBnotWorking_Maeva, struct('file',{pdbId},'error',{err})];
    end

    pdbId = fgetl(file);
end

save(['./membraneSectorDB_Maeva.mat'], 'membraneSectorDB_Maeva');
save(['./membranePDBnotWorking_Maeva.mat'], 'PDBnotWorking_Maeva');
save(['./membranePDBnotWorking_Maeva.mat'], 'PDBworking_Maeva');
clear membraneSectorDB PDBnotWorking PDBworking;
clc;









