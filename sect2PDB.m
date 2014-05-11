function [ sector_PDB ] = sect2PDB( sector_structure )
% this will overlay the sector centroids onto the PDBid and write an output
% pdb

% get coordinates from sector structure
sect1cord = sector_structure.Coordinates
pdbid = sector_structure.Pdb

sector1overlay = pdbwrite(sect1cord, pdbid)

end

