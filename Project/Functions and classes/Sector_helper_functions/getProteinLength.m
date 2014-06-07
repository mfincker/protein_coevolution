function [ proteinLength ] = getProteinLength(data, chainIndex)
%% GETPROTEINLENGHT returns the length of the protein characterized
% by the 'chainIndex' in the pdb structure 'data' entered as input arg.
% Input :
%	- data: pdb structure
%	- chainIndex: index of the chain in the pdb record
	if chainIndex == -1
		proteinLength = -1;
	else
		proteinLength = data.DBReferences(chainIndex).seqEnd - ...
                       		data.DBReferences(chainIndex).seqBegin +1;
    end
end