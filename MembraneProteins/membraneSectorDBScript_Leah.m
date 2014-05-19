%% Script to make a database of sectors from a
% text file containing pdbId, group, subgroup, memberProtein and relatedProtein
% on 5 different lines.


membraneSectorDB_Leah = {};
PDBnotWorking_leah = struct([]);
PDBworking_Leah = {};

file = fopen('./allMembraneProteinsDB_Leah.txt');

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
            membraneSectorDB_Leah = [membraneSectorDB_Leah {Sector(data, clusters{j})}];
        end
        PDBworking_Leah = [PDBworking_Leah; {pdbId}];

    catch err
        PDBnotWorking_leah = [PDBnotWorking_leah, struct('file',{pdbId},'error',{err})];
    end

    pdbId = fgetl(file);
end

save(['./membraneSectorDB_Leah.mat'], 'membraneSectorDB_Leah');
save(['./membranePDBnotWorking_Leah.mat'], 'PDBnotWorking_leah');
save(['./membranePDBnotWorking_Leah.mat'], 'PDBworking_Leah');
clear membraneSectorDB PDBnotWorking PDBworking;
clc;









