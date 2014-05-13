function [ chainID, chainIndex ] = getChainID( data, uniprot )

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