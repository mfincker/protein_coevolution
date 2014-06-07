%% Making SDB from MSA in fasta files called *-aligned.fasta
% where * is a pdbID.
% 05/27/14 - Maeva
% Script to make a database of sectors from a
% folder containing MSA called '*-aligned.fasta' in
% the '../prokaryote_MSA/' directory.


% List all the MSA to process
fileList = dir(['*-aligned.fasta']);

% Create a SDB and the control files for working
% and non working pdbID.
membraneSectorDB_Maeva = {};
PDBnotWorking_Maeva = struct([]);
PDBworking_Maeva = {};

% Loop over the MSA files
for i = 1:numel(fileList)
    name = fileList(i).name;
    disp([char(9) 'File : ' name]);
    
    % Extract pdb ID
    pdbId = regexp(name, '(?<id>.*)-aligned.fasta','names');
    pdbId = pdbId.id;
    
    % Generate msa and clusters
    try
        % Blast2clust_maeva will use the MSA file
        [clusters, data] = blast2clust_maeva_prok(pdbId);
        disp(['For' pdbId ': found ', num2str(length(clusters)), ' clusters!']);

        for j = 1:numel(clusters)
            % Update the SDB
            membraneSectorDB_Maeva = [membraneSectorDB_Maeva Sector(data, clusters{j})];
        end
    
        PDBworking_Maeva = [PDBworking_Maeva; {pdbId}];
    
    catch err
        % Catch an error if it occurs
        PDBnotWorking_Maeva = [PDBnotWorking_Maeva, struct('file',{pdbId},'error',{err})];
    end

    
end

% Save variables
save(['../prokaryote_Sectors/new_prokaryote_SectorDatabase.mat'], 'membraneSectorDB_Maeva');
save(['../prokaryote_Sectors/new_prokaryote_MSAnotWorking.mat'], 'PDBnotWorking_Maeva');
save(['../prokaryote_Sectors/new_prokaryote_MSAworking.mat'], 'PDBworking_Maeva');
clear membraneSectorDB_Maeva PDBnotWorking_Maeva PDBworking_Maeva;




