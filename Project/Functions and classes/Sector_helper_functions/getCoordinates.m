function [ coords ] = getCoordinates ( data, ...
					 	residueInd, chainID)
%GETCOORDINATES returns a 3 x n matrix of the x,y and z coordinates
% of the centroids of each residue in a sector.
% The function takes as arguments:
%	- the pdb matlab structure of the protein of which the sector is part of
%	- the residue indexes
%	- the chainID corresponding to the chain on wich the residues are located

	if (sum(strcmp(chainID,'protein complex') + strcmp(chainID,'undefined'))>0)
		coords = 'undefined';

	else 
		% Get atom numbers that correspond to each residue in the sector
		atomChainIndex = find(strcmp({data.Model.Atom.chainID}, chainID) == 1);
		atomResidue = [data.Model.Atom(atomChainIndex).resSeq];
		atomIndex = {};
		for i = 1:numel(residueInd)
		    atomIndex{i} = find(atomResidue == residueInd(i))';

		end
		% For each residue, get X, Y and Z coordinates of all atoms and
		% calculate the centroid
		centroids = [];
		for i = 1:numel(residueInd)
            if isempty(atomIndex{i})
                atomCoord = [-1 ; -1 ; -1];
            else
                atomCoord = [data.Model.Atom(atomIndex{i}).X ; ...
                             data.Model.Atom(atomIndex{i}).Y ; ...
                             data.Model.Atom(atomIndex{i}).Z];
            end
            centroids = [centroids mean(atomCoord,2)];
		end
		coords = centroids;
	end
end