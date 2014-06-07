%% Script making the necessary variables tracking transmembrane
% residue indexation to make the transmembrane SDB
% 06/01/14 - Maeva
%% Paths need to be modified to run again.

%% Create a transmembrane feature dictionnary:
% Uniprot ID --> Array of residue in the transmembrane region

fFeatures = fopen('transmembraneFeatures.txt', 'r');

transDico = containers.Map();

line = fgetl(fFeatures);
while (strcmp(line,'') == 0 && ischar(line))
    tokens = regexp(line, '(?<ac>\w*)\t(?<start>\d*)\t(?<stop>\d*)','names');
    ac = tokens.ac;
    start = str2num(tokens.start);
    stop = str2num(tokens.stop);
    if isKey(transDico, tokens.ac)
        transDico(tokens.ac) = [transDico(tokens.ac) [start:stop]];
    else
        transDico(tokens.ac) = [start:stop];
    end
    
    line = fgetl(fFeatures);
end

save('transDico.mat', 'transDico');


%% Create dictionnary containing pdbID as key and a matrix containing
% the index of the first residue of the pdb sequence in the pdb record 
% and on uniprot. This will be useful to map uniprot indexation onto the
% pdb indexation.
pdbDico = containers.Map();
pdbIds = getPdb(DB);
for i = 1:numel(pdbIds)
    try
        ac = pdbIds{i};
        disp(ac);
        % Very cool website keepinck track of the indexation !!!
        url = ['http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=pdb&id=' lower(ac) '&all=yes' ];
        s = urlread(url);
        s = s(1,5*49+1:6*49);
        tokens = regexp(s, '(?<pdb>.+)[ ]+(?<chain>\w+)[ ]+(?<pdbSeq>\d+)[ ]+(?<aa>\w+)[ ]+(?<num>\d+?)[ ]+(?<ac>\w*)[ ]+(?<aaLetter>\w+)[ ]+(?<acSeq>\d+)','names');
        num = str2num(tokens.num);
        acSeq = str2num(tokens.acSeq);
        if ~isKey(pdbDico, ac)
            pdbDico(ac) = [num acSeq];
        end
    catch err
        if ~isKey(pdbDico, ac)
            pdbDico(ac) = [1 1];
        end
    end
    
end

save('pdbMapNumberDico.mat', 'pdbDico');

%% Modify the transDico to match the numbering of the sectors
transmembraneSectorCell = {};
pdbDico = importdata('pdbMapNumberDico.mat');
transDico = importdata('transDico.mat');
acNotWorking = {};
pdbNotWorking = {};

for i =1:numel(DB)
    pdb = DB{i}.Pdb;
    uniprot = DB{i}.Uniprot;
    try
        diffNumbering = pdbDico(pdb);
        diffNumbering = diffNumbering(2) - diffNumbering(1);
        try 
            transmembraneResidue = transDico(uniprot) - diffNumbering;
            sectorResidue = DB{i}.ResidueIndexes;
            for res = 1:DB{i}.Length
                residue = sectorResidue(1,res);
                if sum(find(transmembraneResidue == residue)) > 0
                    sectorResidue(1,res) = 1;
                else
                    sectorResidue(1,res) = 0;
                end
            end
            transmembraneSectorCell = [transmembraneSectorCell  {sectorResidue}];
        catch err
            acNotWorking = [ acNotWorking ; {uniprot}];
            transmembraneSectorCell = [transmembraneSectorCell  {zeros(size(sectorResidue))}];
        end
    catch err
        pdbNotWorking = [pdbNotWorking ; {pdb}];
        transmembraneSectorCell = [transmembraneSectorCell  {zeros(size(sectorResidue))}];
    end
end

save('transmembraneSectorCell.mat', 'transmembraneSectorCell');


    


% THE END <--- This is the end baby. The end is right here :) Love you!