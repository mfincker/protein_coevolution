function [sectorcordxyz, pdbSectorCoords] = sect2PDB(sector_structure)
% this will overlay the sector centroids onto the PDBid and write an output
% pdb





% call proper PDB from sector_structure
% pdbtomap = sector_structure.Pdb;
% pdbMapped = getpdb(pdbtomap);
sectorcordxyz = sector_structure.Coordinates


% extract coodinates
sectcord_x = sectorcordxyz(1, :);
sectcord_y = sectorcordxyz(2, :);
sectcord_z = sectorcordxyz(3, :);


% field1 = 'X'; value1 = {sectcord_x};
% field2 = 'Y'; value2 = {sectcord_y};
% field3 = 'Z'; value3 = {sectcord_z};

pdbSectorCoords.X = {sectcord_x}

pdbSectorCoords.Y = {sectcord_y}

pdbSectorCoords.Z = {sectcord_z}

% % get coordinates in xyz
% for i = length(sector_structure)
%     sectcord_X = sectorcordxyz(1, :);
%     sectcord.Y = sectorcordxyz(2, :);
%     sectcord.Z = sectorcordxyz(3, :);
% end

    
    % for i = 1:numel(sector_structure)
%     [sectcordxyz] = sector_pca(sector_structure{1,i})
% end


% % get coordinates from sector structure
%  for i = 1:length(sector_structure)
%      sectcord.X = sector_structure{i}.Coordinates
%      sectcord.Y = sec
%      sectcord.Z =
%      pdbid = sector_structure{i}.Pdb
%  end


end

