function [clusters, data, extra] = blast2clust_maeva(seqacc) 

[msa,seedseq,seednum] = blast2msa_maeva(seqacc,500,1);

% seqacc must be a pdbId in that case !!!!
msa = msa_sort_maeva(seqacc, msa, seedseq);
% BRAD PRUNING FUNCTION GOES HERE
% e.g. msa2 = phyloprunemsa(msa,cutoff)
data = getPdbData(seqacc);

% Hack that sequence interest is seq number 1 !!!
residue_numbers = [1:numel(seedseq)] + data.DBReferences(1).seqBegin - 1 ;

[clusters,extra] = miscwrapper(msa,seedseq,residue_numbers);

end