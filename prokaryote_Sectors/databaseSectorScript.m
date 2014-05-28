%% Script to make a database of sectors from a
% folder containing MSA called '*-aligned.fasta' in
% the '../prokaryote_MSA/' directory.


% MSAdirectories = dir('../prokaryote_MSA_*');
% disp(['Number of directories: ' num2str(numel(MSAdirectories))]);
% 
% for i = 1:numel(MSAdirectories)
%     disp(['Into :' MSAdirectories(i).name]);
%     num = regexp(MSAdirectories(i).name, 'prokaryote_MSA_(?<nb>\d+)','names');
%     num = num.nb;
%     MSAdirectory = ['../' MSAdirectories(i).name '/'];
    fileList = dir(['*-aligned.fasta']);

    membraneSectorDB_Maeva = {};
    PDBnotWorking_Maeva = struct([]);
    PDBworking_Maeva = {};
    
    for i = 1:numel(fileList)
        name = fileList(i).name;
        disp([char(9) 'File : ' name]);
        %         filePath = [MSAdirectory name];
        
        pdbId = regexp(name, '(?<id>.*)-aligned.fasta','names');
        pdbId = pdbId.id;
        
        membraneSectorDB_Maeva = {};
        PDBnotWorking_Maeva = struct([]);
        PDBworking_Maeva = {};
        
        % Generate msa and clusters
        try
            [clusters, data] = blast2clust_maeva_prok(pdbId);
            disp(['For' pdbId ': found ', num2str(length(clusters)), ' clusters!']);

            for j = 1:numel(clusters)
                membraneSectorDB_Maeva = [membraneSectorDB_Maeva Sector(data, clusters{j})];
            end
        
            PDBworking_Maeva = [PDBworking_Maeva; {pdbId}];
        
        catch err
            PDBnotWorking_Maeva = [PDBnotWorking_Maeva, struct('file',{pdbId},'error',{err})];
        end

        
    end

    save(['../prokaryote_Sectors/new_prokaryote_SectorDatabase.mat'], 'membraneSectorDB_Maeva');
    save(['../prokaryote_Sectors/new_prokaryote_MSAnotWorking.mat'], 'PDBnotWorking_Maeva');
    save(['../prokaryote_Sectors/new_prokaryote_MSAworking.mat'], 'PDBworking_Maeva');
    clear membraneSectorDB_Maeva PDBnotWorking_Maeva PDBworking_Maeva;




