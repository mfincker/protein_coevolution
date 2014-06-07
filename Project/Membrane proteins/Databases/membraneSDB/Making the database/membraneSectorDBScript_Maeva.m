%% Defining the sectors from proteins in the pdb
% 05/30/14 - Maeva
% Script to make a database of sectors from a
% text file containing pdbId, group, subgroup, memberProtein and relatedProtein
% on 5 different lines.

% Paths have to be updated before running again.

% List all files containing the proteins information
fileList = dir('./UniqueProteinsDB_Maeva/uniqueMembraneProt_Maeva_*');

disp(['Number of files: ' num2str(numel(fileList))]);

% For each file:
for i = 1:numel(fileList)
    disp([char(9) 'Starting file ' num2str(i)]);
    num = regexp(fileList(i).name, 'uniqueMembraneProt_Maeva_(?<nb>.+)','names');
    num = num.nb;

    % Create a SDB and keep track of the protein working and not working
    membraneSectorDB_Maeva = {};
    PDBnotWorking_Maeva = struct([]);
    PDBworking_Maeva = {};

    % Open the file and start reading it
    file = fopen(fileList(i).name);
    
    % Get the first pdb ID
    pdbId = fgetl(file);
    % As long as you haven't reached the EOF, repeat:
    while (strcmp(pdbId,'') == 0 && ischar(pdbId))
        % Get protein info
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
                % Update the SDB with the new sectors
                membraneSectorDB_Maeva = [membraneSectorDB_Maeva sector];
            end
            % Update the working PDB
            PDBworking_Maeva = [PDBworking_Maeva; {pdbId}];

         catch err
            % In case of an error, keep track of the PDB and the error generated
            PDBnotWorking_Maeva = [PDBnotWorking_Maeva, struct('file',{pdbId},'error',{err})];
         end
        % Process next protein
        pdbId = fgetl(file);
    end

    % Save variables
    save(['./uMembraneSectorDB_Maeva_' num '.mat'], 'membraneSectorDB_Maeva');
    save(['./uMembranePDBnotWorking_Maeva_' num '.mat'], 'PDBnotWorking_Maeva');
    save(['./uMembranePDBWorking_Maeva_' num '.mat'], 'PDBworking_Maeva');
    clear membraneSectorDB_Maeva PDBnotWorking_Maeva PDBworking_Maeva;
    clc;
    % Loop until EOF
end
% Loop until no more file









