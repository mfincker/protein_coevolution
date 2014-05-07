classdef Sector
    %SECTOR class: each instance of this class corresponds to a 
    %sector identified through SCA.
    %   A sector is a group of residues from a given protein that
    %	coevolved together. 
   	%	Property list:
   	%	Method list:
    
    properties
        % Id corresponds to the sector id in the sector database
        % It is by default set to -1 when a sector is instanciated
        Id = -1;

    	% Protein properties
    	Pdb
    	Uniprot
    	ProteinLength
    	OrganismID
    	% Kingdom  % Manually curate it ? Script to get from ncbi
    	% Location % How to get it ?
    	EC

    	% Sector sequence properties
    	Length
    	residueIndexes
    	Sequence
        Coordinates
    	% SecondaryStructure  % Create function to extract this

    end
    
    methods
        % Constructor for the sector class
        % residueresidueIndexes is a row vector !!!
        function sector = Sector(pdbId, residueresidueIndexes)
            % Construction from the online pdb database
            % If no argument is given, construct an empty object
            if nargin > 0
                residueIndexes = sort(residueresidueIndexes);
                sector.Pdb = pdbId;
                % In case your internet connection or PDB is down
                try
                    % Try to get information form
                    data = getpdb(pdbId);

                    % Protein information

                    % Extract uniprot accession number
                    sector.Uniprot = data.DBReferences.dbAccession;
                    % Extract protein length
                    sector.ProteinLength = data.DBReferences.seqEnd - ...
                                        data.DBReferences.seqBegin +1;
                    % Extract organism taxonomical id
                    src = data.Source;
                    % Reformat the character array to be searchable
                    src = src';
                    src = reshape(src,1,numel(src));
                    taxonomyId = regexp(src, 'ORGANISM_TAXID:\s(\d*);', 'tokens');
                    sector.OrganismID = str2num(taxonomyId{1}{1});
                    % Extract protein EC number:
                    src = data.Compound;
                    % Reformat the character array to be searchable
                    src = src';
                    src = reshape(src,1,numel(src));
                    ecNum = regexp(src, 'EC:\s(.*);', 'tokens');
                    sector.EC = ecNum{1}{1};

                    % Sector sequence information

                    sector.Length = numel(residueIndexes);
                    sector.residueIndexes = residueIndexes;
                    % Assuming the residue numbers are the same as the numbers for
                    % the atom list in the PDB file. (Start at 27 for the G6PD)
                    sector.Sequence = data.Sequence.Sequence(residueIndexes - ...
                                        data.DBReferences.seqBegin + 1);

                    % Get atom numbers that correspond to each residue in the sector
                    atomResidue = [data.Model.Atom.resSeq];
                    atomIndex = {};
                    for i = 1:numel(residueIndexes)
                        atomIndex{i} = find(atomResidue == residueIndexes(i))';
                    end
                    sector.Coordinates = atomIndex;




                    % For each residue
                        % Extraxt x, y and z coordinates
                        % Calculate centroid
                        % Store it in a row vector at the same index as the corresponding
                            % residue

                catch err
                    disp(err);
                end
            end
        end

    end
    
end
