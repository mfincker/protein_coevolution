Folder: TP53_code_database
Last Modified Date: Jun 04, 2014

Description: This folder contains the raw data downloaded from the IARC TP53 Database, the BLAST’ed MSA and the code for the analyses.

Files:

clusters.mat: the matlab workspace containing the clusters obtained from the aligned non-human MSA.
is_in_sector.m: a function that accepts a residue number and an array of clusters, and returns the index of the cluster the residue is in.
p53-1000-aligned.fasta: fast file with MSA of non-human P53 sequences
somatic.xlsx: the Excel version of the raw data downloaded from IARC TP53 Database
somaticMutationDataIARC TP53 Database, R17.txt: the raw somatic mutation data downloaded from IARC TP53 Database on May 30, 2014
TP53_process_data: the bulk of the code for analyses is in this .m file, including the visualization and statistics parts