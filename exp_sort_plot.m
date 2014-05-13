%% 
% This script plots the distances from the first 10 MSA sequences to the
% PDB sequence.

files = dir('./prokaryote_MSA/*-aligned.fasta');
names = cell(10,1);
for curr = 1:10
    currName = files(curr,1).name;
    currName = strtok(currName, '-');
    names(curr,1) = cellstr(currName);
    [distSort, I, nSeq] = msa_sort(currName);
    scatter(linspace(1, nSeq, nSeq), distSort);
    hold on
end
xLabel = xlabel('Number of Sequences', 'Fontsize', 16);
yLabel = ylabel('Distance to PDB Sequence', 'Fontsize', 16);
Title = title('Distribution of Distances from MSA Sequences to PDB Sequence', 'Fontsize', 18);
Legend = legend(names, 'Fontsize', 14); legendTitle = get(Legend, 'title');
set(legendTitle, 'string', 'PDB ID');