classdef Sector
    %SECTOR class: each instance of this class corresponds to a 
    %sector identified through SCA. Does not support protein complexes yet !!!!
    %   A sector is a group of residues from a given protein that
    %	coevolved together. 
    %   Ex : sector = Sector('2BH9', [166 170 171 172 173 174 176 183 193 ...
    %                                 198 199 200 201 202 203 204 205 206 ...
    %                                 208 210 215 216 218 225 237 239 242 ...
    %                                 243 245 ])
    %       creates a sector from the pdb file 2BH9 with residues number 166, 170 ...
    %       243, 245.
    %       Once created, you can access the fields like in a structure.
    %       The fields are not protected so be careful not to modify them.
    
    properties
    	% Protein properties
    	Pdb
        ChainID
    	Uniprot
    	ProteinLength
    	OrganismID
        Membrane
    	EC
        Group = '';
        Subgroup = '';
        MemberProteins = {};
        RelatedProteins = {};


    	% Sector sequence properties
    	Length
    	ResidueIndexes
    	Sequence
        Coordinates
    	% SecondaryStructure  % Create function to extract this

    end
    
    methods
        % Constructor for the sector class
        % residueIndexes is a row vector !!!
        function sector = Sector(pdb, residueInd, uniprotNum)
            % Construction from the online pdb database
            % If no argument is given, construct an empty object
            if nargin > 0
                residueInd = sort(residueInd);

                % Create the matlab pdb structure corresponding
                % to the pdb argument:
                if ischar(pdb)
                    data = getPdbData(pdb);
                else
                    data = pdb;
                end


                sector.Pdb = data.Header.idCode;
                    % disp(['PDB id : ' sector.Pdb]);

                uniprot = getUniprot(data);
                if (nargin == 3 & sum(strcmp(uniprotNum, uniprot)) >= 1)
                    sector.Uniprot = uniprotNum;
                elseif nargin == 3 
                    sector.Uniprot = 'protein not in pdb file';
                elseif size(uniprot,2) == 1 
                    sector.Uniprot = uniprot{1};
                else
                    sector.Uniprot = uniprot;
                end
                    % disp(['Uniprot : ' sector.Uniprot]);

                [sector.ChainID, chainIndex] = getChainID(data, sector.Uniprot);
                    % disp(['Chain ID : ' sector.ChainID]);
                    % disp(['Chain index : ' num2str(chainIndex)]);

                sector.ProteinLength = getProteinLength(data, chainIndex);
                    % disp(['Protein length : ' num2str(sector.ProteinLength)]);

                molID = getMolID(data, sector.ChainID);
                    % disp(['Mol ID : ' num2str(molID) ])
                sector.EC = getEC(data, molID);
                    % disp(['EC : ' sector.EC]);
                sector.OrganismID = getOrganismID(data, molID);
                    % disp(['Organism ID : ' num2str(sector.OrganismID)]);

                sector.Membrane = getMembrane(data);
                    % disp(['Membrane : ' num2str(sector.Membrane)]);

                sector.Length = numel(residueInd);
                    % disp(['Length : ' num2str(sector.Length)]);
                sector.ResidueIndexes = residueInd;
                    % disp(['Residue indexes : ' num2str(sector.ResidueIndexes)]);
                sector.Sequence = getSequence(data, sector.ResidueIndexes, ...
                                    chainIndex, sector.ChainID);
                    % disp(['Sequence : ' sector.Sequence]);
                sector.Coordinates = getCoordinates(data, sector.ResidueIndexes, sector.ChainID);
                    % disp('Coords : ');
                    % disp(sector.Coordinates);
            end
        end

    end
    
end
