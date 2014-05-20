%% Script to make a database of sectors from a
% text file containing pdbId, group, subgroup, memberProtein and relatedProtein
% on 5 different lines.

fileList = dir('./uniqueMembraneProteinsDB_Leah*');

disp(['Number of files: ' num2str(numel(fileList))]);
for i = 1:numel(fileList)
    disp([char(9) 'Starting file ' num2str(i)]);
    num = regexp(fileList(i).name, 'uniqueMembraneProteinsDB_Leah(?<nb>\d+).txt','names');
    num = num.nb;

    membraneSectorDB_Leah = {};
    PDBnotWorking_Leah = struct([]);
    PDBworking_Leah = {};


    file = fopen(fileList(i).name);

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
        try
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
             PDBnotWorking_Leah = [PDBnotWorking_Leah, struct('file',{pdbId},'error',{err})];
         end

        pdbId = fgetl(file);
    end

    save(['./membraneSectorDB_Leah' num2str(num) '.mat'], 'membraneSectorDB_Leah');
    save(['./membranePDBnotWorking_Leah' num2str(num) '.mat'], 'PDBnotWorking_Leah');
    save(['./membranePDBWorking_Leah' num2str(num) '.mat'], 'PDBworking_Leah');
    clear membraneSectorDB_Leah PDBnotWorking_Leah PDBworking_Leah;
    clc;
end









