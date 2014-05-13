%% Script to make a database of sectors from a
% folder containing MSA called '*-aligned.fasta' in
% the '../prokaryote_MSA/' directory.


MSAdirectory = '../prokaryote_MSA/';

fileList = dir([MSAdirectory '*-aligned.fasta']);

databaseSector = {};

for i = 1:numel(fileList)
	name = fileList(i).name;
	filePath = [MSAdirectory name];

	%% Importing MSA into MATLAB
	rawdata = fastaread(filePath);
	numseq = length(rawdata);
	seqlength = length(rawdata(1).Sequence);

	%% Convert MSA to numbers
	msa = char(zeros(numseq,seqlength));
	for curseq = 1:numseq
	    % this line converts characters to their number representation
	    msa(curseq,:) = rawdata(curseq).Sequence;

	end
	msa = aa2int(msa);

	% get reference sequence
	pdbDirectory = '../PDBfiles/';
	pdbId = regexp(name, '(?<id>.*)-aligned.fasta','names');
	pdbId = pdbId.id;
    if (exist(pdbDirectory, 'file') == 7 & exist([pdbDirectory pdbId '.pdb'], 'file') == 2)
    	data = pdbread([pdbDirectory pdbId '.pdb']);

    else
      	% Try to get information from PDB online
        % and create a local copy of the pdb file for next time
        data = getpdb(pdbId);
        pdbwrite([pdbDirectory pdbId '.pdb'], data);
    end

	sequence = data.Sequence.Sequence;
	residue_numbers = [1:length(sequence)] + data.DBReferences.seqBegin - 1 ;
	%% get clusters!!
	[clusters,extra] = miscwrapper(msa,sequence,residue_numbers);
	disp(['For' pdbId ': found ', num2str(length(clusters)), ' clusters!']);

	for j = 1:numel(clusters)
		databaseSector = [databaseSector {Sector(data, clusters{j})}];
	end

end

save('./smallSectorDatabase.mat', 'databaseSector')




