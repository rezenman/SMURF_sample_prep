# SMURF_sample_prep
This repository contatins all scripts and auxilliary files to Prepare 16S samples obtained through the swift pipeline for smurf analysis.

packages needed: 
 - cutadapt - here it is installed on a conda environment called cutadaptenv which is activated first - version 2.6
 - R - version used here is 3.5.1 with the following packages need to be installed:
             dada2_1.10.1, ShortRead_1.40.0, Biostrings_2.50.2, stringr_1.4.0, phylotools_0.2.2, seqRFLP_1.0.1 
             
external scripts needed: 
- Config file called "aux_01_config_file" - which list paths to all files needed
- R script called "aux_02_dada2_all_script.R"    
- Script called  "aux_03_fasta_to_fastq.pl"
- Primer files "primers_16S_V1-9_anchored_for.fasta" "primers_16S_V1-9_anchored_rev.fasta"
                         
Inputs needed:
- Path to the folder where read 1 and read 2 in fastq.gz format 
- Path to the config file

Most relevant output: two fastq.gz files - one for read1 and one for read2
log output - a file called "log.txt" contains all logs

to run use: 
```
./Prepare_samples_for_smurf.sh "path_to_folder_where_R1_and_R2" "path_to_config_file"
```
