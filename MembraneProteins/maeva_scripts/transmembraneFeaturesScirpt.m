

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


%% Create a pdb <--> uniprot seq number map
pdbDico = containers.Map();
% pdbAcDico = containers.Map();
% 
% files = dir('pdbsws_res_*');
% 
% for i = 1:numel(files)
%     disp (files(i).name)
%     fMap = fopen(files(i).name, 'r');
% 
%     pdbDico = containers.Map();
%     pdbAcDico = containers.Map();
% 
%     line = fgetl(fMap);
%     while (strcmp(line,'') == 0 && ischar(line))
%         try
%             tokens = regexp(line, '(?<pdb>.+)[ ]+(?<chain>\w+)[ ]+(?<pdbSeq>\d+)[ ]+(?<aa>\w+)[ ]+(?<num>\d+?)[ ]+(?<ac>\w*)[ ]+(?<aaLetter>\w+)[ ]+(?<acSeq>\d+)','names');
%             ac = tokens.ac;
%             pdb = tokens.pdb;
%             pdbSeq = str2num(tokens.pdbSeq);
%             acSeq = str2num(tokens.acSeq);
% 
%             if ~isKey(pdbDico, upper(pdb))
%                 disp(pdb);
%                 pdbDico(upper(pdb)) = [pdbSeq acSeq];
%             end
%         catch err
%         end
% 
%         line = fgetl(fMap);
%     end
% end

pdbIds = getPdb(DB);
for i = 1:numel(pdbIds)
    try
        ac = pdbIds{i};
        disp(ac);
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


%%
% Percentage of the aa in sectors in the transmembrane region:
aaInTrans = 0;
aaInSector = 0;
percentInTransPerSector = zeros(1, numel(DB));
for i = 1:numel(transmembraneSectorCell)
    aaInTrans = sum(transmembraneSectorCell{i}) + aaInTrans;
    aaInSector = aaInSector + numel(transmembraneSectorCell{i});
    percentInTransPerSector(i) = sum(transmembraneSectorCell{i}) / numel(transmembraneSectorCell{i});
end

percentInTrans = aaInTrans / aaInSector;

%%
% Amino acid frequency in the transmembrane part of the sectors
seq = '';
for i = 1:numel(transmembraneSectorCell)
    sequence = DB{i}.Sequence;
    seq = [seq sequence(transmembraneSectorCell{i} == 1)];
end

aaFreqTransmembrane = aacount(seq);
aaFreqTransmembrane = struct2cell(aaFreqTransmembrane);
aaFreqTransmembrane = cell2mat(aaFreqTransmembrane);
aaFreqTransmembrane = aaFreqTransmembrane / sum(aaFreqTransmembrane);
%%
figure('Name', 'Percentage of amino acid in transmembrane region per sector');
hist(percentInTransPerSector);
ylabel('Number of sectors');
xlabel('Percentage of transmembrane amino acids');

%%
% Let's isolate the sectors which have more than 50 % of their aa in the
% transmembrane region.
    


% THE END <--- This is the end baby. The end is right here :) Love you!