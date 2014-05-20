function [ pdbCoords, mergedPDB ] = sectorOverlay( sectorStructure )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% get coordinates
sectorcordxyz = sectorStructure.Coordinates;

%extract xyz
sect_x = sectorcordxyz(1,:)
sect_y = sectorcordxyz(2,:)
sect_z = sectorcordxyz(3,:)

% make structure
field1 = 'X'; value1 = {sect_x};
field2 = 'Y'; value2 = {sect_y};
field3 = 'Z'; value3 = {sect_z};

pdbCoords = struct(field1, value1, field2, value2, field3, value3)

% create PDB from coordinates
mat2pdb(pdbCoords);
sectorPDB = importdata('mat2PDB.pdb');

% overlay with pdb
pdbtomap = sectorStructure.Pdb;
pdbMapped = getpdb(pdbtomap);





end

