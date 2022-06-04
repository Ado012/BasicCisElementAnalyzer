
PURPOSE: 

To examine the relationship between cis elements and transcriptional output on a larger scale the CRM structure of a list of upregulated and downregulated wus genes was analyzed computationally.

Processing begins with the main function Arabmotifsearch(). A target file containing the genes to be analyzed is input along with the desired names of the output files. A GFF file containing annotation information for Arabidopsis genes is also read in. The annoation file is broken into several pieces including chromosome, gene names, and start and end positions. A filestream is started for the output files. 

Targets are extracted from the target list one by one. For each target the function targetScanner() is run associating the target with a chromosome and sequence. The function MotifPrelimScanner() scans the sequence defined from 3000 before the gene start to 3000 after the gene end for TAAT/ATTA cis elements. The cis elements list is fed into ClusterScanner() to detect cis element clusters/CRMs which are defined as strings of at least 4 cis elements which are within 50 bp of the last. 

Each CRM is then scanned for complex cis elements in Complexcorescanner() which is defined as a string of cis elements which each cis elements is 4bp or less from the last one. Motifscorer() elements calculates phasing score which is defined by how well consecutive cis elements adhere to a 10.5x bp spacing relationship. Then complexcore score; which is the length of complex cis elements summed, the phasing per base; where the phasing is divided by the number of bp of the cluster is calculated.

These results are arranged by gene and by cluster into columns and printed out by Motifwriter().  


USEAGE: 

Necessary files: 
target file : ex upregulatedwusgenes_cyclo.csv
Complete set of TAIR10 Chromosome files: ex Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa (from http://ftp.gramene.org/CURRENT_RELEASE/fasta/arabidopsis_thaliana/dna/). Place in Rdatafiles directory. Alternative chromosome fastas can be substituted by modifying the script.
script: BasicCisElementAnalyzer.R. 

Install R and ensure target file is in working directory. Run BasicCisElementAnalyzer.R. 

sample command: arabMotifSearch(targetfile = "Rdatafiles/upregulatedwusgenes.csv", resultsOutput="Rdatafiles/MotifOutputUP.txt", resultsBed="Rdatafiles/MotifOutputBEDUP.bed", resultsSorted="Rdatafiles/MotifOutputSortedUP.txt")

targetfile: list of genes to examine
resultsOutput: desired name of main output file
resultsBed: desired name of output file in BED format.
resultsSorted: desired name sorted file (currently nonfunctional)

chainstart: start of cluster
chainend: end of cluster	
phasescore: phasing score (how well cores adhere to a 10.5x bp spacing 	
phasePerBase: phasing score divided by cluster bp
corNum: number of cis elements
complexCores: number of complex cis elements (multiple cis elements 4bp or less apart) 	
complexCoreScore: length of complex cis elements summed across a cluster. 


Loading BED files
BED files can be loaded through a genome browser such as IGV from the Broad Institute. Simply download and run the browser. Load the TAIR10 Arabidopsis genome and then load the BED file as a custom track. 




