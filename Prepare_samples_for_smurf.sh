#!/bin/bash

#purpose: prepare 16S samples obtained through the swift pipeline for smurf analysis
#packages needed: cutadapt - here it is installed on a conda environment called cutadaptenv which is activated first - version 2.6
#                 R - version used here is 3.5.1 with the following packages need to be installed:
#                     dada2_1.10.1, ShortRead_1.40.0, Biostrings_2.50.2, stringr_1.4.0, phylotools_0.2.2, seqRFLP_1.0.1 
                        
#external scripts needed: the config file called "aux_01_config_file" - which list paths to the following files:
#                         the R script called "aux_02_dada2_all_script.R"    
#                         the script called  "aux_03_fasta_to_fastq.pl"
#                         the primer files "primers_16S_V1-9_anchored_for.fasta" "primers_16S_V1-9_anchored_rev.fasta"
#Input: path to the folder where read 1 and read 2 in fastq.gz format 
#input: path to the config file
#most relevant output: two fastq.gz files - one for read1 and one for read2
#log output - a file called "log.txt" contains all logs

#to run use: ./Prepare_samples_for_smurf.sh "path_to_folder_where_R1_and_R2" "path_to_config_file"

#loading all parameters from the config file
source $2

#loading relevant packages
module load R/3.5.1 #abovementioned relevant packages that need to be installed beforehand
source activate cutadaptenv #use cutadapt version 2.6

folder=$1
cd $folder

read1=$(ls | grep _R1.*fastq.gz$)
read2=$(ls | grep _R2.*fastq.gz$)

##1 - demultiplexing files to 16S variable regions
echo -e "\n-----DEMULTIPLEXING_READS_TO_16S_REGIONS-----\n" >> log.txt

cutadapt \
  -e 0.15 --no-indels \
  -g file:$fwd_primer \
  -G file:$rev_primer\
  -o {name1}-{name2}.1.fastq.gz -p {name1}-{name2}.2.fastq.gz \
  $read1 $read2 &>>log.txt 

comp=$(grep -c -E "=== Summary ===" log.txt)
if (( $comp == 1))
then 
  echo "Successfully demultiplexed files to 16S variable regions" >> log.txt
else 
  echo "ERROR - Something went wrong with demultiplexing files to 16S variable regions, please check log files" >> log.txt && exit
fi

 ##2 - creating a list of only relevant primer combinations for later use and report of how many reads per region 
echo -e "\n-----CREATING_A_LIST_OF RELEVEANT_FILES_TO_USE_FOLLOWING_REGION_DEMULTIPLEXING-----\n" >> log.txt
ls  | grep -E "^V[1_f|3_f|4_f|4_r|5_r].*-V[4_r|5_r|1_f|3_f|4_f]|^V[6_f|7_f|8_r|9_r].*-V[6_f|7_f|8_r|9_r]" | grep -E -v "_r-.*_r|_f.*_f" > valid_combs
if [ $? -eq 0 ] ; then echo "Finished creating valid comb file" >> log.txt ; else echo "Err - couldn't create valid combinations file" >>log.txt && exit  ; fi

##3 - creating report of how many reads per region 
echo -e "\n-----READS_PER_REGION-----\n" >> log.txt
echo "sample_name" "read_count" "total_reads" "percent_of_total" > read_report

total=$(cat log.txt  | grep "Total read pairs processed" | cut -d':' -f2 | sed -e 's/^[ \t]*//' | sed -e 's/,//g')   
for file in $(cat valid_combs| grep .1.fastq.gz | awk '{print $1}')
    do
        read_count=$(zcat $file | grep -c ^@)
        percent=$(echo "scale=5; $read_count / $total" | bc)
        echo $file $read_count $total $percent >> read_report
    done 
    
total_reads=$(awk '{ sum += $2 } END { print sum }' read_report)
total_percent=$(awk '{ sum += $4 } END { print sum }' read_report)

echo "Total" $total_reads $total $total_percent >> read_report

cat read_report >> log.txt

##4 - Moving files to folders called FWD and REV for later use with the dada package
echo -e "\n-----MOVING_AND_ARRANGING_FILES_FOR_DADA2-----\n" >> log.txt

mkdir FWD REV
cat valid_combs | grep -E "_f.*_r.1.fastq.gz|_r.*_f.2.fastq.gz" | awk '{system("mv "$1" FWD/"$1)}'
cat valid_combs | grep -E "_f.*_r.2.fastq.gz|_r.*_f.1.fastq.gz"  | awk '{system("mv "$1" REV/"$1)}'   
ls | grep -E "^[V|unknow].*fastq.gz$" | awk '{system("rm "$1)}'

echo "Finished arranging fastqs successfully" >> log.txt

##5 - concatenate files with the same primer combination but reversed reads - this step might be relevant only to the swift pipeline
echo -e "\n-----CONCATENATING_FILES_WITH_SAME_PRIMER_COMBINATIONS_BUT_REVERSED-----\n" >> log.txt
    cd FWD
    cat V1_f-V4_r.1.fastq.gz V4_r-V1_f.2.fastq.gz > Cat_V1_f-V4_r.1.fastq.gz &
    cat V1_f-V5_r.1.fastq.gz V5_r-V1_f.2.fastq.gz > Cat_V1_f-V5_r.1.fastq.gz &
    cat V3_f-V4_r.1.fastq.gz V4_r-V3_f.2.fastq.gz > Cat_V3_f-V4_r.1.fastq.gz &
    cat V3_f-V5_r.1.fastq.gz V5_r-V3_f.2.fastq.gz > Cat_V3_f-V5_r.1.fastq.gz &
    cat V4_f-V4_r.1.fastq.gz V4_r-V4_f.2.fastq.gz > Cat_V4_f-V4_r.1.fastq.gz &
    cat V4_f-V5_r.1.fastq.gz V5_r-V4_f.2.fastq.gz > Cat_V4_f-V5_r.1.fastq.gz &
    cat V6_f-V8_r.1.fastq.gz V8_r-V6_f.2.fastq.gz > Cat_V6_f-V8_r.1.fastq.gz &
    cat V6_f-V9_r.1.fastq.gz V9_r-V6_f.2.fastq.gz > Cat_V6_f-V9_r.1.fastq.gz &
    cat V7_f-V8_r.1.fastq.gz V8_r-V7_f.2.fastq.gz > Cat_V7_f-V8_r.1.fastq.gz &
    cat V7_f-V9_r.1.fastq.gz V9_r-V7_f.2.fastq.gz > Cat_V7_f-V9_r.1.fastq.gz &
    wait
    cd ..

    cd REV
    cat V1_f-V4_r.2.fastq.gz V4_r-V1_f.1.fastq.gz > Cat_V1_f-V4_r.2.fastq.gz &
    cat V1_f-V5_r.2.fastq.gz V5_r-V1_f.1.fastq.gz > Cat_V1_f-V5_r.2.fastq.gz &
    cat V3_f-V4_r.2.fastq.gz V4_r-V3_f.1.fastq.gz > Cat_V3_f-V4_r.2.fastq.gz &
    cat V3_f-V5_r.2.fastq.gz V5_r-V3_f.1.fastq.gz > Cat_V3_f-V5_r.2.fastq.gz &
    cat V4_f-V4_r.2.fastq.gz V4_r-V4_f.1.fastq.gz > Cat_V4_f-V4_r.2.fastq.gz &
    cat V4_f-V5_r.2.fastq.gz V5_r-V4_f.1.fastq.gz > Cat_V4_f-V5_r.2.fastq.gz &
    cat V6_f-V8_r.2.fastq.gz V8_r-V6_f.1.fastq.gz > Cat_V6_f-V8_r.2.fastq.gz &
    cat V6_f-V9_r.2.fastq.gz V9_r-V6_f.1.fastq.gz > Cat_V6_f-V9_r.2.fastq.gz &
    cat V7_f-V8_r.2.fastq.gz V8_r-V7_f.1.fastq.gz > Cat_V7_f-V8_r.2.fastq.gz &
    cat V7_f-V9_r.2.fastq.gz V9_r-V7_f.1.fastq.gz > Cat_V7_f-V9_r.2.fastq.gz &
    wait
    cd ..

echo "Finished concatenating files successfully" >> log.txt

##6- Running dada2 script
echo -e "\n-----RINNING_DADA2_SCRIPT-----\n" >> log.txt

Rscript $dada_script &>> log.txt

comp=$(grep -c -E "FINISHED ANALYZING" log.txt)
if (( $comp == 1))
then 
  echo "Successfully finished dada2" >> log.txt
else 
  echo "ERROR - Something went wrong with dada2, please check log files" >> log.txt && exit
fi


##7 - modifying files to match the fromat required by smurf

echo -e "\n-----ARRANGING_FILES_FOR_SMURF-----\n" >> log.txt

#replacing all spaces by new lines to create a fasta file
echo "--replacing all spaces to new lines to create fasta file" >> log.txt
for i in $(ls SampCat*.fasta); do sed 's/ /\n/g' $i > "$i"_mod ; done
  
echo "--adding back primers for smurf analysis" >> log.txt

for fa in $(ls *_for.fasta_mod | awk '{print $1}')
  do
    fwd=$(echo $fa | grep -E -o V[0-9]_f)
    seq_for=$(grep -A1 $fwd $fwd_primer | grep -v $fwd | sed 's/\^//g')
    seq_for=$(echo $seq_for | sed 's/M/A/g')
    seq_for=$(echo $seq_for | sed 's/Y/C/g')
    sed '2~2 s/^/'$seq_for'/' $fa > $fa"_primers_added"&  
  done

for fa in $(ls *_rev.fasta_mod | awk '{print $1}')
  do
    rev=$(echo $fa | grep -E -o V[0-9]_r)
    seq_rev=$(grep -A1 $rev $rev_primer | grep -v $rev | sed 's/\^//g')
    seq_rev=$(echo $seq_rev | sed 's/M/A/g')
    seq_rev=$(echo $seq_rev | sed 's/Y/C/g')
    sed '2~2 s/^/'$seq_rev'/' $fa > $fa"_primers_added"& 
  done

echo "--finished adding back primers" >> log.txt
## run the following loop to turn files from fasta to fastq with q-score of 40 for all bases
echo "reformating fasta files as pseudo fastq files" >> log.txt
for i in $(ls *primers_added); do perl $fasta_to_fastq $i > "$i".fastq ; done    
echo "--finished reformatting to fastq files" >> log.txt

##move all files with added primers to a new folder
echo "--moving all files to designated foldes" >> log.txt

mkdir fastq_with_primers
mkdir fasta_files
find . -maxdepth 1 -name "*primers_added.fastq" -exec mv {} "./fastq_with_primers" \; 
find . -maxdepth 1 -name "*.fasta*" -exec mv {} "./fasta_files" \; 

##8 - concatenate all regions to one fastq file with sample name, compress it to create final read1 and read2 files
echo -e "\n-----CONCATENATING_ALL_REGIONS_TO_CREATE_FINAL_FILES-----\n" >> log.txt
cd fastq_with_primers
samp=$(pwd | rev | cut -d'/' -f2 | rev)
name_1=$samp"_L001_R1_001.fastq"
name_2=$samp"_L001_R2_001.fastq"
     
cat *_for*fastq > $name_1 
cat *_rev*fastq > $name_2 

gzip -c  $name_1 > $samp"_R1_swift.gz" 
gzip -c  $name_2 > $samp"_R2_swift.gz"
rm $name_1 $name_2 
cd ..
echo "--finished creating R1 and R2 concatenated files" >> log.txt

echo -e "\n-----FINSIHED_PREPARING_SAMPLE_FOR_SMURF_ANALYSIS-----\n" >> log.txt
echo -e "--final sample names are\n" $samp"_R1_swift.gz\n" $samp"_R2_swift.gz"  >> log.txt