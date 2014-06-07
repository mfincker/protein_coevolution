%% Defining the sectors from proteins in the pdb
% 05/30/14 - Maeva
% Script to make a database of sectors from a
% text file containing pdbId, group, subgroup, memberProtein and relatedProtein
% on 5 different lines.

% Paths have to be updated before running again.

% List all files containing the proteins information
fileList = dir('./splitMembrane_*');

disp(['Number of files: ' num2str(numel(fileList))]);

% For each file:
for i = 1:numel(fileList)
    disp([char(9) 'Starting file ' num2str(i)]);
    num = regexp(fileList(i).name, 'splitMembrane_(?<nb>.+)','names');
    num = num.nb;

    % Create a SDB and keep track of the protein working and not working
    membraneSectorDB_Leah = {};
    PDBnotWorking_Leah = struct([]);
    PDBworking_Leah = {};

    % Open the file and start reading it
    file = fopen(fileList(i).name);

    % Get the first pdb ID
    pdbId = fgetl(file);

    % As long as you haven't reached the EOF, repeat:
    while ischar(pdbId)
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
                membraneSectorDB_Leah = [membraneSectorDB_Leah {sector}];
            end
            % Update the working PDB
            PDBworking_Leah = [PDBworking_Leah; {pdbId}];

         catch err
            % In case of an error, keep track of the PDB and the error generated
             PDBnotWorking_Leah = [PDBnotWorking_Leah, struct('file',{pdbId},'error',{err})];
         end

        % Process next protein
        pdbId = fgetl(file);
    end

    % Save variables
    save(['./membraneSectorDB_Leah_' num '.mat'], 'membraneSectorDB_Leah');
    save(['./membranePDBnotWorking_Leah_' num '.mat'], 'PDBnotWorking_Leah');
    save(['./membranePDBWorking_Leah_' num '.mat'], 'PDBworking_Leah');
    clear membraneSectorDB_Leah PDBnotWorking_Leah PDBworking_Leah;
   % Loop until EOF
end
% Loop until no more file









