Folder: TP53_results
Last Modified Date: Jun 04, 2014

Description: This folder contains the completed figures and matlab workspace file containing the statistical data for the analyses.

Files:

Figure 1 enrichment_p53: the mutation enrichments of each sector.

Figure 2 comp_mut_count: the count of sequences with mutations on residues in each sector

Figure 3 clust_mutation_res133: visualization by images of the non-human residues in sector 3 in attempt to identify compensatory mutations within the sector

processed_data.mat: the saved matlab workspace after running the analyses in TP53_process_data of the TP53_code_database folder. This file does not include the last analysis that visualize figure 3.

stats.mat: the saved matlab workspace after running the statistics on enrichment of figure 1 and counts of mutated sequences in figure 2.

compensatory mutation list.mat: the list of compensatory mutations found using our strategy is stored in comp_mut, with the first four columns and last four containing information of the residue index, wild type AA, mutant AA (using MATLAB’s aa2int), and the cluster index. Column 5 is the sequences with the compensatory mutation.

Supplemental Tables TP53.docx: this word file contains the three supplemental tables of statistic tests done on the TP53 data.

Comp Mutation Table TP53.docx: this word file contains the table of compensatory mutation candidates generated from the database.