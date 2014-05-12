%% clear all
clear all
%% Importing MSA into MATLAB
rawdata = fastaread('2BH9-500-aligned-sequences-P11413.fasta');
numseq = length(rawdata);
seqlength = length(rawdata(1).Sequence);

%% Convert MSA to numbers
msa = char(zeros(numseq,seqlength));
for curseq = 1:numseq
    % this line converts characters to their number representation
    msa(curseq,:) = rawdata(curseq).Sequence;
end
msa = aa2int(msa);

msa = uint32(msa);
%% get reference sequence
refID = 'P11413';
gensearch = getgenpept(refID);
sequence = gensearch.Sequence;
residue_numbers = 1:length(sequence);
%% get clusters!!
[clusters,extra] = miscwrapper(msa,sequence,residue_numbers);
disp(['Found ', num2str(length(clusters)), ' clusters!']);
