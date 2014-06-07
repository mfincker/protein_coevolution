% Record the number of sequences used in the MSA pruned by Brad's function
% and save the result for further pruning of the SDB

% List all MSA to process
d =dir('./*-500-aligned.fasta');

% Create a structure that will keep track
% of the name of the MSA and the number of sequences
% in the MSA
sizeMSA = struct([]);
smallMSA = {};
errMSA = struct([]);

% For each file:
for i = 1:numel(d)
    rawdata = fastaread(d(i).name);
    numseq = length(rawdata);
    seqlength = length(rawdata(1).Sequence);

    %Convert MSA to numbers
    msa = char(zeros(numseq,seqlength));
    for curseq = 1:numseq
        % this line converts characters to their number representation
        msa(curseq,:) = rawdata(curseq).Sequence;
    end
    
    % Get the reference sequence to which distance will be calculated
    acc = regexp(d(i).name, '(?<id>.+)-500-aligned.fasta','names');
    acc = acc.id;
    seedseq = getgenpept(acc);
    seedseq = seedseq.Sequence;
    residue_numbers = 1:length(seedseq);
    try
        % Calculate the distance from every sequence in the MSA to the reference seq
        msaKept = msa_sort_maeva(acc, msa, seedseq);
        sizeMSA = [sizeMSA, struct('PDB', acc, 'MSAsize', size(msaKept, 1))];
        disp([acc ' : ' num2str(size(msaKept, 1))]);
        if size(msaKept,1) < 50
            smallMSA = [smallMSA {acc}];
        end
    catch err
        errMSA = [errMSA struct('PDB', acc, 'error', {err})];
        disp([acc ' : ' 'error']);
    end
    
end

% Save results
save('new_uMembraneSizeMSA_Leah.mat', 'sizeMSA');
save('new_uMembraneSmallMSA_Leah.mat', 'smallMSA');
save('new_uMembraneMSAerrors_Leah.mat', 'errMSA');
    