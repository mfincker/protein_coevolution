function [clusters, data, extra] = blast2clust_maeva(seqacc) 
%% Modified version of blast2clust that uses msa_sort_maeva
% to constrain the phylogenetic distance covered in the MSA
% and blast2msa_maeva to check for raw and aligned MSA fasta 
% files before blasting on NCBI. 
% See blast2clust, blast2msa amd msa_sort_maeva for more help

[msa,seedseq,seednum] = blast2msa_maeva(seqacc,500,1);

% seqacc must be a pdbId in that case !!!!
msa = msa_sort_maeva(seqacc, msa, seedseq);
% BRAD PRUNING FUNCTION GOES HERE
% e.g. msa2 = phyloprunemsa(msa,cutoff)
data = getPdbData(seqacc);

% Hack that sequence interest is seq number 1 !!!
residue_numbers = [1:numel(seedseq)] + data.DBReferences(1).seqBegin - 1 ;

disp([char(9) 'Into miscwrapper']);
[clusters,extra] = miscwrapper(msa,seedseq,residue_numbers);

end