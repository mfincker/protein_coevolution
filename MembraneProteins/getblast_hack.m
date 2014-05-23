function blastdata = getblast(rid,varargin)
%GETBLAST retrieves BLAST report files from NCBI.
%
%   GETBLAST parses NCBI BLAST reports, including BLASTN, BLASTP, BLASTX,
%   TBLASTN, TBLASTX and PSI-BLAST. BLAST reports offer a fast and
%   powerful comparative analysis of interesting protein and nucleotide
%   sequences against known structures in existing on-line databases.
%
%   BLASTDATA = GETBLAST(RID) reads a BLAST report RID, returning the data
%   in the report as a structure, DATA. The NCBI Request ID (RID) must be a
%   recently generated report - NCBI purges reports after 24 hours.
%
%   BLASTDATA = GETBLAST(..., 'DESCRIPTIONS', NUM_DESC) includes the
%   specified number of descriptions in the report. Acceptable values are
%   1-500, and the default is 100.
%
%   BLASTDATA = GETBLAST(..., 'ALIGNMENTS', NUM_ALGNMNTS) includes the
%   specified number of alignments in the report. Acceptable values are
%   1-500, and the default is 50.
%
%   BLASTDATA = GETBLAST(..., 'TOFILE', FILENAME) saves the data returned
%   from the NCBI BLAST report in a file, FILENAME. The default format for
%   the file is text, but you can specify HTML with the FILEFORMAT argument.
%
%   BLASTDATA = GETBLAST(..., 'FILEFORMAT', FORMAT) returns the report in
%   the specified format, FORMAT.  Acceptable values are 'TEXT' and 'HTML'.
%
%   BLASTDATA = GETBLAST(..., 'WAITTIME', TIME) pauses MATLAB up to the
%   specified TIME (in minutes) to wait for the report from NCBI. If the
%   report is still not available at that time, it will return an error
%   message. The default behavior is not to wait for the report.
%
%   Examples:
%
%       % Run a BLAST search with an NCBI accession number.
%       RID = blastncbi('AAA59174','blastp','expect',1e-10)
%
%       % Then pass the RID to GETBLAST to parse the report, load it into
%       % a MATLAB structure, and save a copy as a text file.
%       report = getblast(RID,'TOFILE','Report.txt')
%
%   For more information about reading and interpreting BLAST reports, see
%   http://www.ncbi.nlm.nih.gov/Education/BLASTinfo/Blast_output.html.
%
%   See also BLASTFORMAT, BLASTLOCAL, BLASTNCBI, BLASTREAD, BLASTREADLOCAL.

%   Copyright 2004-2008 The MathWorks, Inc.


%   Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A.
%   Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman
%   (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein
%   database search programs",  Nucleic Acids Res. 25:3389-3402.

% verify java is available
if (~usejava('jvm'))
    error(message('bioinfo:getblast:NoJava'))
end

% Set default parameters to retrieve the BLAST report
site = 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&NEW_FORMATTER=on&RID=';
alignview = '&ALIGNMENT_VIEW=Pairwise';
desc = '&DESCRIPTIONS=100';
align = '&ALIGNMENTS=50';
showgi = '&NCBI_GI=on';
format = '&FORMAT_TYPE=text';
tofile = false;
filename = '';
fileformat = 'text';
wait_allow = 0; 

% Check optional input arguments
if nargin > 1
    if rem(nargin,2) == 0
        error(message('bioinfo:getblast:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'descriptions','alignments','tofile','fileformat',...
        'waittilready','waittime'};
    for j=1:2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname, okargs,numel(pname))); 
        % changed waittilready to waituntilready. Allow old code to work.
        if isequal(k(:),[5;6])
            k = 6;
        end
        if isempty(k)
            error(message('bioinfo:getblast:UnknownParameterName', pname));
        elseif length(k)>1
            error(message('bioinfo:getblast:AmbiguousParameterName', pname));
        else
            switch(k)
                case 1  % descriptions
                    if (isnumeric(pval) && (ismember(pval,1:500)))
                        desc = ['&DESCRIPTIONS=' int2str(pval)];
                    else
                        error(message('bioinfo:getblast:BadDescParam'));
                    end
                case 2  % alignments
                    if (isnumeric(pval) && (ismember(pval,1:5000)))
                        align = ['&ALIGNMENTS=' int2str(pval)];
                    else
                        error(message('bioinfo:getblast:BadAlignParam'));
                    end
                case 3  % tofile
                    if ischar(pval)
                        tofile = true;
                        filename = pval;
                    end
                case 4 % fileformat
                    okformats = {'text','html'};
                    val = strmatch(lower(pval),okformats); 
                    if length(val) == 1
                        fileformat = okformats{val};
                    else
                        if isempty(val)
                            error(message('bioinfo:getblast:BadFormat'));
                        else
                            error(message('bioinfo:getblast:AmbiguousFormat', pval));
                        end
                    end
                case {5,6}  % waittilready
                    if ischar(pval)
                        pval = str2double(pval);
                    end
                    if (isnumeric(pval) && (pval > 0) && isreal(pval))
                        wait_allow = 60*pval;
                    else
                        error(message('bioinfo:getblast:wait'));
                    end
            end
        end
    end
end

% Check for a valid RID type and format
if iscellstr(rid)
    rid=char(rid);
end
%if regexpi(rid,'\d+-\d+-\d+\.BLASTQ\d')
%if regexpi(rid, '\<(\w{11})\>') % change in RID format April 16 2007
if regexpi(rid, '\<(\w+)\>') % change in RID format April 16 2008
 
    htmlurl = [site rid alignview desc align showgi];
    texturl = [htmlurl format];
else
    error(message('bioinfo:getblast:InvalidRID', rid));
end

% polls the NCBI BLAST SERVER for results
searchurl = [site rid '&FORMAT_OBJECT=SearchInfo'];
searchInfo = urlread(searchurl);

start = regexp(searchInfo, 'QBlastInfoBegin');
stop = regexp(searchInfo, 'QBlastInfoEnd');
if isempty(start) || isempty(stop)
    error(message('bioinfo:getblast:BlastInfoError', rid, searchurl));
end

status = regexpi(searchInfo(start(1):stop(1)), '(?<=Status=)\w+', 'match', 'once');

if strcmp(status, 'WAITING')
    disp 'Blast results are not available yet. Please wait ...';
    wait_timer = 0;
    while wait_timer <= wait_allow && strcmp(status, 'WAITING')
        pause(10);
        wait_timer = wait_timer + 10;
        searchInfo = urlread(searchurl);
        start = regexp(searchInfo, 'QBlastInfoBegin');
        stop = regexp(searchInfo, 'QBlastInfoEnd');
        status = regexpi(searchInfo(start(1):stop(1)), '(?<=Status=)\w+', 'match', 'once');
    end
end

if strcmp(status, 'WAITING')% examine status
    error(message('bioinfo:getblast:BlastSearchWaiting', rid));
elseif strcmp(status, 'FAILED')
    error(message('bioinfo:getblast:BlastSearchError', rid));   
elseif strcmp(status, 'UNKNOWN')
    error(message('bioinfo:getblast:BlastSearchExpired', rid));
elseif strcmp(status, 'READY')
    blasttext = urlread(texturl);
    blasttext = strrep(blasttext,char(0),'a'); % Correct invalid character in returned text
    %RID = regexpi(blasttext, '(?<=RID\s*:\s*)\w{11}','match','once'); % new RID format, 16 April 2007
    if isempty(regexpi(searchInfo(start(2):stop(2)), 'ThereAreHits=yes'))
        warning(message('bioinfo:getblast:BlastSearchNull', rid));
        blastdata = [];
        return;
    else % hits found
        if tofile == true %  Write out file?
            if strcmpi(fileformat,'html')
                try
                    [fpath,fstatus] = urlwrite(htmlurl, filename); %#ok
                catch theException
                    if strcmpi(theException.identifier,'MATLAB:urlwrite:InvalidOutputLocation')
                        error(message('bioinfo:getblast:CouldNotOpenHTML', filename));
                    else
                        msgId = 'bioinfo:getblast:HTMLWriteError';
                        newException = MException(msgId,getString(message(msgId)));
                        throw(addCause(newException,theException))
                    end
                end
            else
                fid = fopen(filename,'w');
                if fid == (-1)
                    error(message('bioinfo:getblast:CouldNotOpenForWriting', filename));
                else
                    fwrite(fid, blasttext,'char');
                    fclose(fid);
                end
            end
       end
    end
else
    error(message('bioinfo:getblast:BlastUnknownStatus', rid));
end
% Parse the text file output
blastdata = blastread(blasttext);
