function [ pdbStruct ] = getPdbData( pdbId )

	% pdb is a char array containing a pdb id
	if ischar(pdbId)
		PDBdirectory = '../PDBfiles/';
        if (exist(PDBdirectory, 'file') == 7 && ...
        	exist([PDBdirectory pdbId '.pdb'], 'file') == 2)
            pdbStruct = pdbread([PDBdirectory pdbId '.pdb']);
        else
%         	try
	            % Try to get information from PDB online
	            % and create a local copy of the pdb file for next time
	            pdbStruct = getpdb(pdbId);
                % To avoid flooding the server
                pause(0.001);
	            pdbwrite([PDBdirectory pdbId '.pdb'], pdbStruct);
% 	        catch err
% 	        	rethrow(err);
% 	        end
        end

    % pdb is already a matlab pdb struct of a pdb file
	else
		pdbStruct = pdb;
	end
end