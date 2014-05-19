%% Script to make a database of sectors from a
% text file containing pdbId, group, subgroup, memberProtein and relatedProtein
% on 5 different lines.


membraneSectorDB = {};
PDBnotWorking = struct([]);
PDBworking = {};

file = fopen('./dbTest.txt');

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
            membraneSectorDB = [membraneSectorDB {Sector(data, clusters{j})}];
        end
        PDBworking = [PDBworking; {pdbId}];

    % catch err
    %     PDBnotWorking = [PDBnotWorking, struct('file',{pdbId},'error',{err})];
    % end

    pdbId = fgetl(file);
end

save(['./membraneSectorDB.mat'], 'membraneSectorDB');
save(['./membranePDBnotWorking.mat'], 'PDBnotWorking');
save(['./membranePDBworking.mat'], 'PDBworking');
clear membraneSectorDB PDBnotWorking PDBworking;
clc;









