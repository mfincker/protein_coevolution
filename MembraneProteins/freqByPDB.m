%% this gets the total sequence for the entire protein!
% get pdbsequences for toxin proteins
seq_7HAL = getpdb('7AHL', 'SequenceOnly', 'True');
%%
seq_1PFO = getpdb('1PFO', 'SequenceOnly', 'True');
%%
seq_4HSC = getpdb('4HSC', 'SequenceOnly', 'True');
%%
seq_3B07 = getpdb('3B07', 'SequenceOnly', 'True');
%%
seq_3O44 = getpdb('3O44', 'SequenceOnly', 'True');
%%
% These are aminoacid counts per amino acid in the whole protein
counts_7HAL = aacount(seq_7HAL{1});
counts_1PF0 = aacount(seq_1PFO);


%%
% get sequences from sectors
sectorSeq7HAL={}
for i = 1:numel(tox1_Sectors)
    sectorSeq7HAL{i}= tox1_Sectors{1,i}.Sequence
end

% join into one string to count later
totSectorSeq7HAL = strjoin(sectorSeq7HAL);
%% trying to get correct frequencies for sectors

% count AA for amino acids in all sectors in a given protein here we use
% 7HAL
aacountSect7HAL = aacount(totSectorSeq7HAL);

% get counts per amino acid into cell array
cellSect7HAL = struct2cell(aacountSect7HAL);

% total number of residues in sectors
totalsectorresi7HAL = sum(cell2mat(cellSect7HAL));

% get frequencies of sector residues as a fraction of total sector residues

sumsectorFreq_7HAL = sum(aaSectorDBCount(tox1_Sectors),2)

% THIS IS THE FREQUENCY WE WANT. Frequency of given AA in sectors in all
% the residues in all the sectors in a given protein.
freqAAinSect = sumsectorFreq_7HAL/totalsectorresi7HAL

% measure relative enrichment
% for whole protein
% get counts per amino acid into cell array
cell7HAL = cell2mat(struct2cell(counts_7HAL));

% get total number of residues in protein
residuecount = sum(cell7HAL)

% frequency of residues in whole protein
totProtFreq = cell7HAL./residuecount

% total number of residues in sectors
totalsectorresi7HAL = sum(cell2mat(cellSect7HAL))

    
%% get relative frequency
relativeFreq = freqAAinSect./totProtFreq
relativeFreq(isnan(relativeFreq)) = 1;
relativeFreq7HAL= relativeFreq




