function [ pdbStruct ] = getPdbData( pdb )

	% pdb is a char array containing a pdb id
	if ischar(pdb)
		PDBdirectory = '../PDBfiles/';
        if (exist(PDBdirectory, 'file') == 7 && ...
        	exist([PDBdirectory pdb '.pdb'], 'file') == 2)
            pdbStruct = pdbread([PDBdirectory pdb '.pdb']);
        else
        	try
	            % Try to get information from PDB online
	            % and create a local copy of the pdb file for next time
	            pdbStruct = getpdb(pdb);
	            pdbwrite([PDBdirectory pdb '.pdb'], pdbStruct);
	        catch error
	        	rethrow(err);
	        end
        end

    % pdb is already a matlab pdb struct of a pdb file
	else
		pdbStruct = pdb;
	end
end