function [seqDistPDB_sort, I, nCloseSeq] = msa_sort(pdbID)

%% Sorting MSA by Similarity to PDB Sequence
% This function takes in a string of PDB ID, sorts the sequences in the
% corresponding MSA by similarity to the PDB sequence.

%%
% Input:
%%
% 1. pdbID: the PDB database ID (I thought about adding an option to
% extract the sequence from the local database, but thought it may be
% easier to just extract from online.)
%%
% Output:
%%
% 1. A histogram of the distances from the sequences in the MSA to the
% PDB sequence. The distance is measured with scoring matrix BLOSUM60.
%%
% 2. seqDistPDB_sort: the sorted array of distances between MSA sequences
% and the PDB sequence.
%%
% 3. I: the array containing the original index of MSA sequences before 
% sorting.

%% Get PDB Sequence into a cell from online
protein_pdb = getpdb(pdbID);
pdb_sequence = cellstr(protein_pdb.Sequence.Sequence);

%% Import MSA
dirName = sprintf('./prokaryote_MSA/%s-aligned.fasta',pdbID);
rawData = fastaread(dirName);
nSeq = length(rawData);

%% Convert MSA into Cell Array of strings of characters
msa = cell(nSeq,1);
for curr = 1:nSeq
    msa(curr,:) = cellstr(rawData(curr).Sequence);
end

%% Compare each sequence in MSA to the PDB sequence
seqDistPDB = zeros(nSeq,1);
nCloseSeq = 1;
for curr = 1:nSeq
    currSeq = msa(curr,:);
    % Prune the gaps from the MSA sequences
    currSeq = strrep(currSeq, '-', '');
    dist = seqpdist(vertcat(pdb_sequence,currSeq),'SCORINGMATRIX','BLOSUM60');
    if dist <= 10
        seqDistPDB(nCloseSeq,:) = dist;
        nCloseSeq = nCloseSeq + 1;
    end
end
[seqDistPDB_sort, I] = sort(seqDistPDB);

%% The number of close sequences to PDB sequence
% Only counting the sequences with distances smaller than 10, because the
% value may be as large as 30 and dwarf all the other values, making the
% scatterplot mostly an horizontal line.

nCloseSeq = nSeq;
for curr = 1:length(seqDistPDB_sort)
    if seqDistPDB_sort(curr,:) >= 10
        nCloseSeq = nCloseSeq - 1;
    end
end