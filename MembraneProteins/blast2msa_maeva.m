function [msa,seedseq,residue_numbers] = blast2msa(seqacc,varargin) 
% BLAST2MSA     This function returns a multiple sequence alignment given a
%               sequence or accession number.  A second argument sets the
%               numebr of sequences to return (default 500). A third
%               argument (1/0) determines whether the fasta and aligned
%               fasta files are deleted after the msa is generated. 

%%
disp('Starting blast..');
acc = seqacc;
keepfile = 0;
tokeep = 500;
if nargin == 1
    % do nothing.
elseif nargin==2
    tokeep = varargin{1};
elseif nargin==3    
    tokeep = varargin{1};
    keepfile = varargin{2};
else
    error('MATLAB:TooManyInputs','Too many input arguments.');
end

%% downloading original sequence
disp('Getting original sequence');
seedseq = getgenpept(acc);
seedgi = seedseq.GI;
seedseq = seedseq.Sequence;
residue_numbers = 1:length(seedseq);

disp(['Blasting accession ', acc]);
rid = blastncbi(seedseq,'blastp','Alignments',tokeep);

%% download blast
clear blastresults
disp('Attempting to download blast results.');
returnblast = 0;
while ~returnblast
    try
        blastresults = getblast_hack(rid,'Alignments',tokeep);
        returnblast = 1;
    catch
        for i = 1:52
            fprintf('\b');
        end
        fprintf('\b.');
    end
end
disp('Success!');
%% get GI values from BLAST
ids = '';
for l = 1:length(blastresults.Hits)
    splits = find(blastresults.Hits(l).Name=='|');
    gi = blastresults.Hits(l).Name(4:(splits(2)-1));
    ids = [ids,',',gi];
end
ids = ids(2:end);
%% download their full FASTA sequences
disp('Downloading sequences of blast hits');
divides = find(ids==',');
towriteout = '';
filenameblast = [acc,'-',num2str(tokeep),'-rawblast.fasta'];

fid = fopen(filenameblast,'w');
stepby = 1:300:(length(divides)-1);
if stepby(end)~=(length(divides)-1)
    stepby(end+1) = (length(divides)-1);
end
for l = 1:(length(stepby)-1)
    curids = ids((divides(stepby(l))+1):(divides(stepby(l+1))-1));
    url = ['http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&id=',curids];
    toadd = urlread(url);
    nofasta = strfind(toadd,'Supplied id parameter is empty.');
    if ~isempty(nofasta)
        toadd = toadd(1:(nofasta(1)-1));
    end
    fprintf(fid,toadd);
end
fclose(fid);
%% align and such
disp('Aligning sequences.');
filenameout = [acc,'-',num2str(tokeep),'-aligned.fasta'];
[curdir] = fileparts(which(mfilename));
findspaces = find(curdir == ' ');
for l = length(findspaces):-1:1
    curdir = [curdir(1:(findspaces(l)-1)),'\',curdir(findspaces(l):end)];
end
[~,~] = system([curdir,'/clustalo -i ', filenameblast, ' -o ', filenameout]);
%% get msa
pause(1);
disp('Generating MSA');
rawdata = fastaread(filenameout);
numseq = length(rawdata);
seqlength = length(rawdata(1).Sequence);

%% Convert MSA to numbers
msa = char(zeros(numseq,seqlength));
for curseq = 1:numseq
    % this line converts characters to their number representation
    msa(curseq,:) = rawdata(curseq).Sequence;
end
% To work with Brad function
% msa = aa2int(msa);

%% clean up if necessary
disp('Cleaning up!');
if ~keepfile
    delete(filenameout);
    delete(filenameblast);
end
%% Done
disp('Done!');