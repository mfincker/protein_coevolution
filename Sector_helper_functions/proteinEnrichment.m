function [ enrichment, totalsectorresi] = proteinEnrichment( proteinsectorDB )
% takes protein sector DB and returns amino acid enrichment
%   Detailed explanation goes here

% get the sequences from all sectors in a given protein
sectorSeq={}
for i = 1:numel(proteinsectorDB)
    sectorSeq{i}= proteinsectorDB{i}.Sequence;
end
% join into one string to count later
totSectorSeq = strjoin(sectorSeq);

% count AA for amino acids in all sectors in a given protein here we use
% 7HAL
aacountSect = aacount(totSectorSeq);

% get counts per amino acid into cell array
cellSect = struct2cell(aacountSect);

% total number of residues in sectors
totalsectorresi = sum(cell2mat(cellSect));


% get frequencies of sector residues as a fraction of total sector residues

sumsectorFreq = sum(aaSectorDBCount(proteinsectorDB),2);

% THIS IS THE FREQUENCY WE WANT. Frequency of given AA in sectors in all
% the residues in all the sectors in a given protein.
freqAAinSect = sumsectorFreq/totalsectorresi;

pdbID = proteinsectorDB{1}.Pdb;

%this is the entire protein sequence
totalSeq = getpdb(pdbID, 'SequenceOnly', 'True');


if strcmp(class(totalSeq), 'cell')
    totalSeq = totalSeq{1};
end


% get length of total seq
lengthtotalSeq = length(totalSeq);

aaCountProtein = cell2mat(struct2cell(aacount(totalSeq)));

% sequence frequency in protein

proteinFreq = aaCountProtein/lengthtotalSeq;

% define enrichement
enrichment = freqAAinSect./proteinFreq;
enrichment(isnan(enrichment)) = 1;


end

