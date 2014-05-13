%% Script to make a database of sectors from a
% folder containing MSA called '*-aligned.fasta' in
% the '../prokaryote_MSA/' directory.


MSAdirectory = '../prokaryote_MSA_copy/';

fileList = dir([MSAdirectory '*-aligned.fasta']);

databaseSector = {};
MSAnotWorking = {};
MSAworking = {};

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
	pdbId = regexp(name, '(?<id>.*)-aligned.fasta','names');
	data = getPdbData(pdbId.id);

	sequence = data.Sequence.Sequence;
	residue_numbers = [1:length(sequence)] + data.DBReferences.seqBegin - 1 ;
	%% get clusters!!
    try
        [clusters,extra] = miscwrapper(msa,sequence,residue_numbers);
        disp(['For' pdbId.id ': found ', num2str(length(clusters)), ' clusters!']);
        MSAworking = [MSAworking {filePath}];
    catch err
        MSAnotWorking = [MSAnotWorking {filePath}];
    end

	for j = 1:numel(clusters)
		databaseSector = [databaseSector {Sector(data, clusters{j})}];
	end

end

save('./prokaryote_MSA_SectorDatabase.mat', 'databaseSector');
save('./prokaryote_MSAnotWorking.mat', 'MSAnotWorking');
save('./prokaryote_MSAworking.mat', 'MSAworking');




