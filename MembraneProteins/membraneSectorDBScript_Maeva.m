%% Script to make a database of sectors from a
% text file containing pdbId, group, subgroup, memberProtein and relatedProtein
% on 5 different lines.

fileList = dir('./UniqueProteinsDB_Maeva/uniqueMembraneProt_Maeva_*');
% fileList = 'UniqueProteinsDB_Maeva/uniqueMembraneProt_Maeva_ak';

disp(['Number of files: ' num2str(numel(fileList))]);
% disp(['File name : ' fileList]);

for i = 1:numel(fileList)
    disp([char(9) 'Starting file ' num2str(i)]);
    num = regexp(fileList(i).name, 'uniqueMembraneProt_Maeva_(?<nb>.+)','names');
%     num = regexp(fileList, 'uniqueMembraneProt_Maeva_(?<nb>.+)','names');
    num = num.nb;

    membraneSectorDB_Maeva = {};
    PDBnotWorking_Maeva = struct([]);
    PDBworking_Maeva = {};


    file = fopen(fileList(i).name);
%     file = fopen(fileList);
    
    pdbId = fgetl(file);
    while (strcmp(pdbId,'') == 0 && ischar(pdbId))
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
                membraneSectorDB_Maeva = [membraneSectorDB_Maeva sector];
            end
            PDBworking_Maeva = [PDBworking_Maeva; {pdbId}];

         catch err
             PDBnotWorking_Maeva = [PDBnotWorking_Maeva, struct('file',{pdbId},'error',{err})];
         end

        pdbId = fgetl(file);
    end

    save(['./uMembraneSectorDB_Maeva' num '.mat'], 'membraneSectorDB_Maeva');
    save(['./uMembranePDBnotWorking_Maeva' num '.mat'], 'PDBnotWorking_Maeva');
    save(['./uMembranePDBWorking_Maeva' num '.mat'], 'PDBworking_Maeva');
    clear membraneSectorDB_Maeva PDBnotWorking_Maeva PDBworking_Maeva;
    clc;
end









