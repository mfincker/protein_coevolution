function [ pdbStruct ] = getPdbData( pdbId )
%% GETPDBDATA returns the pdb matlab structure corresponding to
% the pdb id entered as input argument.
% The function will try to see if the pdb record file exist locally
% in a PDBfiles directory located 2 nodes up from the working directory.
% If the directory exists, as well as the pdb record, the pdb matlab structure
% will be created from the local record. If the folder exists, but not the record,
% the function will download the record online, create the structure and keep the
% file in the directory. If the directory does not exist, it will just create the
% structure from online PDB.

	% pdb is a char array containing a pdb id
	if ischar(pdbId)
		try
			% Look for a PDBfiles folder 2 nodes up from the working directory
			PDBdirectory = '../../PDBfiles/';
	        if (exist(PDBdirectory, 'file') == 7 && ...
	        	exist([PDBdirectory pdbId '.pdb'], 'file') == 2)
	            pdbStruct = pdbread([PDBdirectory pdbId '.pdb']);
	        else
		            % Try to get information from PDB online
		            % and create a local copy of the pdb file for next time
		            pdbStruct = getpdb(pdbId);
	                % To avoid flooding the server
	                pause(0.001);
		            pdbwrite([PDBdirectory pdbId '.pdb'], pdbStruct);
	        end
    	catch err
			% Try to get information from PDB online
            % and create a local copy of the pdb file for next time
            pdbStruct = getpdb(pdbId);
            % To avoid flooding the server
            pause(0.001);
            % No keeping of the pdb file because the PDBfiles folder was not found.
         end

    % pdb is already a matlab pdb struct of a pdb file
	else
		pdbStruct = pdb;
	end
end