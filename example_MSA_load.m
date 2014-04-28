%% Importing MSAs into MATLAB
rawdata = fastaread('1UD2-aligned.fasta');
numseq = length(rawdata);
seqlength = length(rawdata(1).Sequence);


%% Convert them to numbers
msa = zeros(numseq,seqlength);
for curseq = 1:numseq
    % this line converts characters to their number representation
    msa(curseq,:) = double(rawdata(curseq).Sequence);
end
subplot(1,2,1); imagesc(msa)
xlabel('MSA Position');
ylabel('Sequence Number');
%% Prune low-information
toremove = zeros(1,seqlength);
for curpos = 1:seqlength
    % remove gap characters
    if mode(msa(:,curpos))== double('-')
        toremove(curpos) = 1;
    end
    % remove completely conserved sequences
    if nnz(msa(:,curpos) == msa(1,curpos)) == seqlength
        toremove(curpos) = 1;
    end
end
msaprune = msa;
msaprune(:,find(toremove)) = [];
subplot(1,2,2); imagesc(msaprune);
xlabel('Trimmed Position');
ylabel('Sequence');
%% Calculate entropy of each position
% (its like a calculation of possible states)
seqlength2 = size(msaprune,2);
posentropy = zeros(1,seqlength2);

for curpos = 1:seqlength2
    curdata = msaprune(:,curpos);
    states = unique(curdata);
    numstates = length(states);
    probstates = zeros(1,numstates);
    for l = 1:length(probstates)
        probstates(l) = nnz(curdata==states(l))/numseq;
    end
    posentropy(curpos) = -sum(probstates.*log(probstates));
end
figure;
hist(posentropy,20); xlabel('Entropy of position'); ylabel('Frequency');
%% calculate joint entropy (see wikipedia)
% like above, but for two positions)
jentropy = zeros(seqlength2,seqlength2);
for curpos1 = 1:seqlength2
    for curpos2 = curpos1:seqlength2
        curdata = msaprune(:,[curpos1,curpos2]);
        states = unique(curdata,'rows');
        numstates = size(states,1);
        probstates = zeros(1,numstates);
        for l = 1:length(probstates)
            probstates(l) = nnz((curdata(:,1)==states(l,1)).*(curdata(:,2)==states(l,2)))/numseq;
        end
        probstates = probstates.*log(probstates);
        probstates(isnan(probstates)) = 0;
        jentropy(curpos1,curpos2) = -sum(probstates);
        jentropy(curpos2,curpos1) = jentropy(curpos1,curpos2);
    end
end
figure;
xlabel('Position'); ylabel('Position');
imagesc(jentropy)