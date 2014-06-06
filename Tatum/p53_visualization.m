
%% Sectors to visualize. The first four are PDBs with DNA and the last four
%are without. 3Q05 is pdb for tetramer p53 with DNA, and 3Q01 is tetramer
%without DNA.
sector1wDNA = Sector('3Q05',clusters{1,1});
sector2wDNA = Sector('3Q05',clusters{1,2});
sector3wDNA = Sector('3Q05',clusters{1,3});
sector4wDNA = Sector('3Q05',clusters{1,4});

sector1wo = Sector('3Q01',clusters{1,1});
sector2wo = Sector('3Q01',clusters{1,2});
sector3wo = Sector('3Q01',clusters{1,3});
sector4wo = Sector('3Q01',clusters{1,4});

%% 
%In order to visualize sector 5 we had to trim first and last residue in
%order to make a visualizable sector. The first one is without DNA and the
%second one is with DNA
% sector5_untrimmed = clusters{1,5};
% sector5_trimmed = sector5_untrimmed(1:end-1);
% sector5wDNA = Sector('3Q05',sector5_trimmed);
% 
% sector5wo = ector('3Q01',sector5_trimmed);

%% Generate Visualizable pdb files
% p53_sector_stucture_wDNA is one struct and p53_sector_stucture_wo is
% another

p53_sector_stucture_wDNA = cell(1,4);
p53_sector_stucture_wo = cell(1,4);

p53_sector_stucture_wDNA{1,1} = sector1wDNA;
p53_sector_stucture_wDNA{1,2} = sector2wDNA;
p53_sector_stucture_wDNA{1,3} = sector3wDNA;
p53_sector_stucture_wDNA{1,4} = sector4wDNA;
% p53_sector_stucturewDNA{1,5} = sector5wDNA;
%Now without
% p53_sector_stucture_wo{1,1} = sector1wo;
% p53_sector_stucture_wo{1,2} = sector2wo;
% p53_sector_stucture_wo{1,3} = sector3wo;
% p53_sector_stucture_wo{1,4} = sector4wo;
% % p53_sector_stucturewo{1,5} = sector5wo;


sectorsOverlay(p53_sector_stucture_wDNA);
% 
% sectorsOverlay(p53_sector_stucture_wo);


%Sector 1 visualized as blue, Sector 2 as green, sector 3 as yellow and
%sector 4 as red.