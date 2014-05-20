function resmap = msaSeqRef(msa,refseq,resnum)
% MSASEQREF      aligns a sequence to an MSA (rows=sequence), and returns
%                the mapping of column indices in the MSA to the residue
%                positions of seq. 
%
% "msa" is a numseq x numpos integer matrix corresponding to a multiple
% sequence alignment of the protein of interest
%
% "seq" is the sequence you would like to map to the MSA. It can be a
% partial sequence, but make sure to correctly specify resnum. It can be a
% character array or an array of integers returned by aa2int().
%
% "resnum" is a vector of the same size as "seq", with each index
% corresponding to the corresponding residue number in "seq". This is
% useful when you only have a partial sequence available.
%
% resmap = MSASEQREF(msa, seq, resnum) returns a vector of size "numpos" -
% the number of positions in the msa. Positions in the MSA that do not have 
% a corresponding position in the sequence of interest are NaN.

% msa concensus sequence:
s1 = int2aa(mode(double(msa),1));

% reference sequence:
if ischar(refseq)
    s2 = refseq;
else
    s2 = int2aa(refseq);
end

% global align:
[~,al] = nwalign(s2,s1);


% count up the matches (: or |), throw away gaps
topseq = cumsum(aa2int(al(1,:))<=20);
topseq(find((aa2int(al(1,:))>20)+(al(2,:)==' '))) = 0;
% convert map according to resnum
resmap = zeros(1,size(msa,2));
resmap(topseq~=0) = resnum(topseq(topseq~=0));
end