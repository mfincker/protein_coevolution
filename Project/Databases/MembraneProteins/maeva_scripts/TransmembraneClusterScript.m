%% Analysis of the 5 transmembrane sector clusters

DBtrans75_1 = importdata('DBtrans75_1.mat');
DBtrans75_2 = importdata('DBtrans75_2.mat');
DBtrans75_3 = importdata('DBtrans75_3.mat');
DBtrans75_4 = importdata('DBtrans75_4.mat');
DBtrans75_5 = importdata('DBtrans75_5.mat');

%%
DBtrans75 = [DBtrans75_1 DBtrans75_2 DBtrans75_3 DBtrans75_4 DBtrans75_5];

%% Amino acid frequency and enrichment
aaFreqDB75_1 = aaSectorDBCount(DBtrans75_1);
aaFreqDB75_1 = sum(aaFreqDB75_1,2)/ sum(sum(aaFreqDB75_1,2));

aaFreqDB75_2 = aaSectorDBCount(DBtrans75_2);
aaFreqDB75_2 = sum(aaFreqDB75_2,2)/ sum(sum(aaFreqDB75_2,2));

aaFreqDB75_3 = aaSectorDBCount(DBtrans75_3);
aaFreqDB75_3 = sum(aaFreqDB75_3,2)/ sum(sum(aaFreqDB75_3,2));

aaFreqDB75_4 = aaSectorDBCount(DBtrans75_4);
aaFreqDB75_4 = sum(aaFreqDB75_4,2)/ sum(sum(aaFreqDB75_4,2));

aaFreqDB75_5 = aaSectorDBCount(DBtrans75_5);
aaFreqDB75_5 = sum(aaFreqDB75_5,2)/ sum(sum(aaFreqDB75_5,2));

aaFreqDB75 = aaSectorDBCount(DBtrans75);
aaFreqDB75 = sum(aaFreqDB75,2)/ sum(sum(aaFreqDB75,2));

%%


%% Amino acid count;
aaCount75_1 = sum(aaSectorDBCount(DBtrans75_1),2);
aaCount75_2 = sum(aaSectorDBCount(DBtrans75_2),2);
aaCount75_3 = sum(aaSectorDBCount(DBtrans75_3),2);
aaCount75_4 = sum(aaSectorDBCount(DBtrans75_4),2);
aaCount75_5 = sum(aaSectorDBCount(DBtrans75_5),2);


%% Transmembrane residues
transmembraneSectorCell = importdata('transmembraneSectorCell.mat');
DB = importdata('../MembraneSectorDB/membraneSectorDB_060114_forCoord.mat'); 
%Amino acid frequency in the transmembrane part of the sectors
seq = '';
for i = 1:numel(transmembraneSectorCell)
    sequence = DB{i}.Sequence;
    seq = [seq sequence(transmembraneSectorCell{i} == 1)];
end

aaFreqTransmembrane = aacount(seq);
aaFreqTransmembrane = struct2cell(aaFreqTransmembrane);
aaFreqTransmembrane = cell2mat(aaFreqTransmembrane);
aaFreqTransmembrane = aaFreqTransmembrane / sum(aaFreqTransmembrane);

%%
figure('Name', 'Amino acid frequency in the 5 clusters')
subplot(3,1,1)
bar([aaFreqDB75_3,aaFreqDB75_5]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);
title('Amino acid frequency');

subplot(3,1,2);
bar([ (aaFreqDB75_3-aaFreqTransmembrane) ./ aaFreqTransmembrane,(aaFreqDB75_5-aaFreqTransmembrane) ./ aaFreqTransmembrane]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);
title('Relative enrichment to the transmembrane residues');

subplot(3,1,3);
% bar([ (aaFreqDB75_3-aaFreqDB75_5) ./ aaFreqDB75_5,(aaFreqDB75_5-aaFreqDB75_3) ./ aaFreqDB75_3]);
% set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
%     'M','F','P','S','T','W','Y','V'},'XTick',1:20);
% legend('Cluster 3 vs Cluster 5', 'Cluster 5 vs Cluster 3');
% title('Relative enrichment ');
chance = ones(20,1)/20;
bar([ (aaFreqDB75_3-chance) ./ chance,(aaFreqDB75_5-chance) ./ chance]);
set(gca,'XTickLabel',{'A','R','N','D','C','Q','E','G','H','I','L','K', ...
    'M','F','P','S','T','W','Y','V'},'XTick',1:20);
title('Relative enrichment to uniform distribution');



