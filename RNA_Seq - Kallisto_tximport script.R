#Heidi Town 26th May

#Aim combine samples to get gene expression from transcript level, importing data from Kallisto - seperating treated and untreated samples. 
#Adapted from Philippa Borrill Script 02_tximport_summarise_counts_tpm_per_gene 
getwd()

library(tximportData)
library(readr)

#Conditions used
study_list <-c("Treated","Control")
study_list

#Read in Tx2gene table (transcript to gene table)
#Tx2gene table made using: 
#url <-"https://raw.githubusercontent.com/Borrill-Lab/WheatFlagLeafSenescence/master/data/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt"
#destfile <- file.path(mydir,"transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt")
#destfile
#download.file(url,destfile)

tx2gene <- read.table("transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt")
head(tx2gene)
getwd()
#header=T means first row contains column names 
tx2gene <- read.table("C:/Users/town/OneDrive - Norwich Bioscience Institutes/RNA_Seq analysis_ treated and untreated seperated at tximport/transcripts_to_genes_RefSeqv1.0_annot_v1.1.txt",header=T)

getwd()

#locations of Kallisto files for treated samples
mydir <- getwd()
samples <- read.csv("Treated_sample_index.csv", header=TRUE)
samples
files <-file.path(mydir,samples$SampleName,"abundance.tsv")
files
names(files)<-samples$SampleName
head(files)
all(file.exists(files))

library(tximport)
txi <- tximport(files, type="kallisto", tx2gene=tx2gene)
names(txi)

#saving counts per summarised gene - treated
setwd("C:/Users/town/OneDrive - Norwich Bioscience Institutes/RNA_Seq analysis_ treated and untreated seperated at tximport/expression per gene")
head(txi$counts)
colnames(txi$counts)
write.table(txi$counts, file="Treated_timecourse_count.tsv", sep ="\t")

#saving tpm summarised per gene 
head(txi$abundance)
colnames(txi$abundance)
write.table(txi$abundance, file="Treated_timecourse_tpm.tsv", sep = "\t")

#lengths per summarised gene
head(txi$length)

#calculate avg. gene length across samples 
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
write.csv(gene_lengths,file="Treated_timecourse_gene_lengths.csv")



#locations of Kallisto files for control samples
setwd("C:/Users/town/OneDrive - Norwich Bioscience Institutes/RNA_Seq analysis_ treated and untreated seperated at tximport")
getwd()
mydir <- getwd()
samples <- read.csv("Control_sample_index.csv", header=TRUE)
samples
files <-file.path(mydir,samples$SampleName,"abundance.tsv")
files
names(files)<-samples$SampleName
head(files)
files
all(file.exists(files))

library(tximport)
txi<- tximport(files, type="kallisto", tx2gene=tx2gene)
names(txi)

#saving counts per summarised gene - control
setwd("C:/Users/town/OneDrive - Norwich Bioscience Institutes/RNA_Seq analysis_ treated and untreated seperated at tximport/expression per gene")
head(txi$counts)
colnames(txi$counts)
write.table(txi$counts, file="control_timecourse_count.tsv", sep ="\t")

#saving tpm summarised per gene 
head(txi$abundance)
colnames(txi$abundance)
write.table(txi$abundance, file="control_timecourse_tpm.tsv", sep = "\t")


#lengths per summarised gene
head(txi$length)

#calculate avg. gene length across samples 
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
write.csv(gene_lengths,file="control_timecourse_gene_lengths.csv")
