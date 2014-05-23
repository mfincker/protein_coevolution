function [ chainID, chainIndex ] = getChainID( data, uniprot )
%GETCHAINID returns the chain Id number of a protein, characterized 
% by its uniprot accession number, in the given pdb matlab struct.
% It also returns the index corresponding to the uniprot number in
% the list of uniprot accession numbers listed in the pdb struct.
%	Arguments:
%		- data : pdb matlab struct
%		- uniprot : char array of a uniprot number
		
	if ischar(uniprot) & strcmp(uniprot, 'protein not in pdb file')
		chainID = 'undefined';
		chainIndex = -1;

	else if ischar(uniprot)
		chainIndex = find(strcmp(uniprot, getUniprot(data)) == 1, 1);
        chainID = data.DBReferences(chainIndex).chainID;

    else
    	chainID = 'protein complex';
    	chainIndex = -1;

    end

end