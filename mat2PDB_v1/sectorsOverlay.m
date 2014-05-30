function sectorsOverlay(sectorDB)
% Protein database per protein
%   Detailed explanation goes here

 field1 = 'X'; 
 field2 = 'Y'; 
 field3 = 'Z';
 field4 = 'betaFactor';
 field5 = 'outfile';
 
 value1 = [];
 value2 = [];
 value3 = [];
 value4 = [];
 value5 = [sectorDB{1}.Pdb,'_sectors.pdb'];
 for i = 1:length(sectorDB)
    sectorStructure = sectorDB{i};

    % get coordinates
    sectorcordxyz = sectorStructure.Coordinates;

    %extract xyz
    sect_x = sectorcordxyz(1,:);
    sect_y = sectorcordxyz(2,:);
    sect_z = sectorcordxyz(3,:);

    % make structure
    value1 = [value1 sect_x];
    value2 = [value2 sect_y];
    value3 = [value3 sect_z];
    value4 = [value4 i*ones(1,length(sect_x))];
    


    end
    
pdbCoords = struct(field1, value1, field2, value2, field3, value3, field4, value4, field5, value5);
% create PDB from coordinates
mat2pdb(pdbCoords);







end

