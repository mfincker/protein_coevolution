%% Script to make a database of sectors from a
% text file containing pdbId, group, subgroup, memberProtein and relatedProtein
% on 5 different lines.

fileList = dir('./splitMembrane_*');

disp(['Number of files: ' num2str(numel(fileList))]);
for i = 1:numel(fileList)
    disp([char(9) 'Starting file ' num2str(i)]);
    num = regexp(fileList(i).name, 'splitMembrane_(?<nb>.+)','names');
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
                membraneSectorDB_Leah = [membraneSectorDB_Leah {sector}];
            end
            PDBworking_Leah = [PDBworking_Leah; {pdbId}];

         catch err
             PDBnotWorking_Leah = [PDBnotWorking_Leah, struct('file',{pdbId},'error',{err})];
         end

        pdbId = fgetl(file);
    end

    save(['./membraneSectorDB_Leah_' num '.mat'], 'membraneSectorDB_Leah');
    save(['./membranePDBnotWorking_Leah_' num '.mat'], 'PDBnotWorking_Leah');
    save(['./membranePDBWorking_Leah_' num '.mat'], 'PDBworking_Leah');
    clear membraneSectorDB_Leah PDBnotWorking_Leah PDBworking_Leah;
    clc;
end









