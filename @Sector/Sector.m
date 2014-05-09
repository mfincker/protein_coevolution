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
        % Id corresponds to the sector id in the sector database
        % It is by default set to -1 when a sector is instanciated
        Id = -1;

    	% Protein properties
    	Pdb
        ChainID
    	Uniprot
    	ProteinLength
    	OrganismID
    	% Kingdom  % Manually curate it ? Script to get from ncbi
        Membrane
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
        % residueIndexes is a row vector !!!
        function sector = Sector(pdbId, residueInd, uniprotNum)
            % Construction from the online pdb database
            % If no argument is given, construct an empty object
            if nargin > 0
                residueInd = sort(residueInd);
                sector.Pdb = pdbId;
                % In case your internet connection or PDB is down
                try
                    % Try to get information form
                    data = getpdb(pdbId);

                    % Protein information

                    % Extract uniprot accession number(s) of the PDB file
                    uniprot = (regexp([data.DBReferences.database],'UNP   ')-1)/6 +1;
                    uniprotAccession = {};
                    for i = 1:size(uniprot,2)
                        uniprotAccession{i} = data.DBReferences(i).dbAccession(1:6);
                    end

                    % monomer
                    if size(uniprotAccession,2) == 1 
                        sector.Uniprot = uniprotAccession{1};
                        sector.ChainID = 'A';

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

                        % Extract membrane location
                        src = data.Keywords;
                        % Reformat the character array to be searchable
                        src = src';
                        src = reshape(src,1,numel(src));
                        membrane = [strfind(src, 'MEMBRANE') strfind(src, 'membrane')];
                        sector.Membrane = ( size(membrane,2) > 0); 

                        % Extract protein EC number:
                        src = data.Compound;
                        % Reformat the character array to be searchable
                        src = src';
                        src = reshape(src,1,numel(src));
                        ecNum = regexp(src, 'EC:\s(.*);', 'tokens');
                        sector.EC = ecNum{1}{1};

                        % Sector sequence information

                        sector.Length = numel(residueInd);
                        sector.residueIndexes = residueInd;
                        % Assuming the residue numbers are the same as the numbers for
                        % the atom list in the PDB file. (Start at 27 for the G6PD)
                        sector.Sequence = data.Sequence.Sequence(residueInd - ...
                                            data.DBReferences.seqBegin + 1);

                        % Get atom numbers that correspond to each residue in the sector
                        atomResidue = [data.Model.Atom.resSeq];
                        atomIndex = {};
                        for i = 1:numel(residueInd)
                            atomIndex{i} = find(atomResidue == residueInd(i))';

                        end
                        % For each residue, get X, Y and Z coordinates of all atoms and
                        % calculate the centroid
                        centroids = [];
                        for i = 1:numel(residueInd)
                            atomCoord = [data.Model.Atom(atomIndex{i}).X ; ...
                                         data.Model.Atom(atomIndex{i}).Y ; ...
                                         data.Model.Atom(atomIndex{i}).Z];
                            centroids = [centroids mean(atomCoord,2)];
                        end
                        sector.Coordinates = centroids;


                    % protein complex
                    else 
                        % Check if a uniprot number as been entered to identify the
                        % protein of interest in the complex
                        if (nargin < 3)
                            % Extract protein length
                            sector.ProteinLength = -1;
                            sector.OrganismID = -1;
                            sector.Uniprot = uniprotAccession;

                            % Extract membrane location
                            src = data.Keywords;
                            % Reformat the character array to be searchable
                            src = src';
                            src = reshape(src,1,numel(src));
                            membrane = [strfind(src, 'MEMBRANE') strfind(src, 'membrane')];
                            sector.Membrane = ( size(membrane,2) > 0); 
                            sector.EC = 'protein complex';

                            % Sector sequence information

                            sector.Length = numel(residueInd);
                            sector.residueIndexes = residueInd;
                            % Assuming the residue numbers are the same as the numbers for
                            % the atom list in the PDB file. (Start at 27 for the G6PD)
                            sector.Sequence = 'protein complex: need a uniprot number';
                            sector.Coordinates = 'protein complex: need a uniprot number'; 

                        else
                            % Extract chainId:
                            % In a multi-mer with multiple occurences of the same protein,
                            % the first occurence of the protein is chosen as chainID.
                            chainIndex = find(strcmp(uniprotNum,uniprotAccession) == 1, 1);
                            sector.ChainID = data.DBReferences(chainIndex).chainID;

                            % Uniprot accession number
                            sector.Uniprot = uniprotNum;

                            % Extract protein length
                            sector.ProteinLength = data.DBReferences(chainIndex).seqEnd - ...
                                            data.DBReferences(chainIndex).seqBegin +1;

                            % Extract EC number if it exists
                                % Find mol_id corresponding to the chainID:
                            cmpd = data.Compound;
                            cmpd = cmpd';
                            cmpd = reshape(cmpd,1,numel(cmpd));
                            cmpdFlip = fliplr(cmpd);
                            molID = regexp(cmpdFlip,['(' sector.ChainID '|' sector.ChainID ...
                                    '[A-Za-a0-9 ,]*?)\s:NIAHC.*?;(?<mol_id>\d*)\s:DI_LOM'],'names');
                            % should only return a single mol_id, but might bug if same chainID in two 
                            % molecules ...
                            molID = fliplr(molID.mol_id);
                                % EC number if it exists:
                            cmpdSplit = regexp(cmpd,'MOL_ID: ','split');
                            if strncmp(cmpdSplit{str2double(molID) + 1}, [molID ';'], ...
                                        size(molID, 2) + 1) == 1
                                ec = regexp(cmpdSplit{str2double(molID) + 1}, ...
                                    'EC:\s(?<ec>\d+\.\d+\.\d+\.\d+);', 'names');
                                if numel(ec) == 1
                                    sector.EC = ec.ec;
                                else
                                    sector.EC = 'undefined';
                                end
                            end

                            % Extract taxonomical data
                                % Match the mol_id and the taxonomical data:
                            src = data.Source;
                            % Reformat the character array to be searchable
                            src = src';
                            src = reshape(src,1,numel(src));


                            sector.OrganismID = -1;

                            % Extract membrane location
                            src = data.Keywords;
                            % Reformat the character array to be searchable
                            src = src';
                            src = reshape(src,1,numel(src));
                            membrane = [strfind(src, 'MEMBRANE') strfind(src, 'membrane')];
                            sector.Membrane = ( size(membrane,2) > 0); 

                            % Sector sequence information
                            sector.Length = numel(residueInd);
                            sector.residueIndexes = residueInd;
                            % Assuming the residue numbers are the same as the numbers for
                            % the atom list in the PDB file. (Start at 27 for the G6PD)
                            sector.Sequence = 'protein complex: need chain id';
                            sector.Coordinates = 'protein complex: need chain id';
                        end
                    end


                catch err
                    rethrow(err);
                    
                end
            end
        end

    end
    
end
