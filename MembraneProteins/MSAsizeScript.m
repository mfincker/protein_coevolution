% Record the number of sequences used in the MSA pruned by Brad function.

d =dir('./*-500-aligned.fasta');
sizeMSA = struct([]);
smallMSA = {};
errMSA = struct([]);

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
    
    acc = regexp(d(i).name, '(?<id>.+)-500-aligned.fasta','names');
    acc = acc.id;
    seedseq = getgenpept(acc);
    seedseq = seedseq.Sequence;
    residue_numbers = 1:length(seedseq);
    try
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

save('uMembraneSizeMSA', 'sizeMSA');
save('uMembraneSmallMSA', 'smallMSA');
save('uMembraneMSAerrors', 'errMSA');
    