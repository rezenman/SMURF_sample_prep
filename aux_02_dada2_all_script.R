print("^>^>^>^>^>^>^>^>^> First of all, session info: ")
print(sessionInfo())
print("^>^>^>^>^>^>^>^>^>")

getwd()
library(stringr)
library(ShortRead)
library(Biostrings)
library(phylotools)
library(seqRFLP)
library(dada2)

wd = getwd()

path <- wd

# File parsing
pathF <- file.path(wd, "FWD/") 
pathR <- file.path(wd, "REV/") 
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...

fastqFs <- sort(list.files(pathF, pattern="Cat.*fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="Cat.*fastq.gz"))


if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

print("Starting filtering process, please wait :)")
out = filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
                    rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
                    truncLen=c(120,120), maxEE=c(2,2), truncQ=2, maxN=0,
                    compress=TRUE, verbose=TRUE)
print("Finished filtering process, thank you for waiting :)")

out_df = as.data.frame(out)
colnames(out_df)
out_df$percent = round(out_df$reads.out / out_df$reads.in, digits = 2)
print(out_df)
print("FINISHED FILTERING")

##################################################################

filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)


sample.names <- sapply(strsplit(basename(filtFs), "\\.[1|2]\\.fastq\\.gz"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "\\.[1|2]\\.fastq\\.gz"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz


if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names


set.seed(100)

# Learn forward error rates3
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)


# Sample inference and merger of paired-end reads
fwd_list = vector("list", length(sample.names))
rev_list = vector("list", length(sample.names))
mergers = vector("list", length(sample.names))
names(fwd_list) <- sample.names
names(rev_list) <- sample.names
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR, justConcatenate = T)
    fwd_list[[sam]] = ddF
    rev_list[[sam]] = ddR
    mergers[[sam]] = merger
}

# Construct sequence table and remove chimeras
save.image(file = "infer_variants.Rdata")

print("FINISHED INFERING VARIANTS")
#############################################################################

keep_index = NULL
for(i in 1:length(mergers)){
  df = mergers[[i]]
  if(dim(df)[1] > 0){keep_index[i] = i}
}
ind_to_keep = as.numeric(na.omit(keep_index))
mergers = mergers[ind_to_keep]



for(sample in 1:length(mergers)){
  ##subsetting all data frames to contain only relevant coulumns - sequence and abundance
  print(paste("analyzing sample:", names(mergers[sample])))
  m = mergers[[sample]][,c(1,2)]
  m$cluster_num = 1:length(m$abundance)  #adding cluster number to each unique combination of for and rev
  m = m[,c(2,3,1)] #reordering the data frame
  ##replicating each merged sequence times its abundance
  multi_m = m[rep(1:nrow(m), m$abundance),]
  print(paste("replicating sequneces went well?", nrow(multi_m) == sum(m$abundance)))
  
  ##splitting merged sequence to Forward and Reverse and reveres complementing the reverse
  multi_m$forward = sapply(strsplit(as.character(multi_m$sequence), "NNNNNNNNNN"), `[`, 1)
  multi_m$reverse = sapply(strsplit(as.character(multi_m$sequence), "NNNNNNNNNN"), `[`, 2)
  multi_m$reverse = as.character(reverseComplement(DNAStringSet(multi_m$reverse)))
  
  ##adding Ids to each read
  multi_m$for_name = paste0(">", 1:length(multi_m$forward),"_Samp", names(mergers[sample]), "_1")
  multi_m$rev_name = paste0(">", 1:length(multi_m$forward),"_Samp", names(mergers[sample]), "_2")
  
  ##creating seperate data frame for reverse and forward
  for_df = multi_m[,c(6,4)]
  rev_df = multi_m[,c(7,5)]
  
  ##writing data frames as fasta sequences
  file_name_for = paste0("Samp", names(mergers[sample]), "_for.fasta")
  file_name_rev = paste0("Samp", names(mergers[sample]), "_rev.fasta")
  write.table(for_df, file = file_name_for, quote = F, row.names = F, col.names = F)
  write.table(rev_df, file = file_name_rev, quote = F, row.names = F, col.names = F)
  
}

print("FINISHED ANALYZING")