## Running DEseq2 Treated vs Control 
## Aim to carry out differential expression analysis comparing samples treated with CHA, Croisor100  and control samples - and plotting a PCA plot. 
## Run after using tximport to merge counts for TPM from Kallisto data. 
## Script adapted from Philippa Borrill Script '03_WT_vs_RNAi_DESeq2_goseq.R'

install.packages("BiocManager")
library(BiocManager)

#Check version
BiocManager::version()

#Install DESeq2 packages
BiocManager::install("DESeq2")

library("DESeq2")
library("dplyr")
library("tidyr")
library(ggplot2)

#set data directory for expression per gene from tximport

data_dir <- "U:/RNA seq/RNA_Seq analysis_ treated and untreated seperated at tximport/expression per gene/"
setwd("U:/RNA seq/RNA_Seq analysis_ treated and untreated seperated at tximport/DESeq2")
#Check working directory
#getwd()

#read in data and select genes with >0.5 transcript per million (>0.5 high confidence genes)
Control.tpm.data <- read.csv(file=paste0(data_dir,"control_timecourse_tpm.tsv"), sep="\t")
head(Control.tpm.data)
Treated.tpm.data <- read.csv(file=paste0(data_dir,"Treated_timecourse_tpm.tsv"), sep="\t")
head(Treated.tpm.data)                            

#selecting genes >0.5 tpm control (>0.5 filtered to select high confidence genes)
control.tpm.data <- Control.tpm.data[,1:15]
head(control.tpm.data)

#Average per timepoint
control.tpm.data$T0 <- (control.tpm.data[,1] + control.tpm.data[,2] + control.tpm.data[,3])
control.tpm.data$T1 <- (control.tpm.data[,4] + control.tpm.data[,5] + control.tpm.data[,6])
control.tpm.data$T2 <- (control.tpm.data[,7] + control.tpm.data[,8] + control.tpm.data[,9])
control.tpm.data$T3 <- (control.tpm.data[,10] + control.tpm.data[,11] + control.tpm.data[,12])
control.tpm.data$T4 <- (control.tpm.data[,13] + control.tpm.data[,14] + control.tpm.data [,15])

# rename columns to keep Average per timepoint
colnames(control.tpm.data)[16:20]
tmp_data_avg_con <- control.tpm.data[,16:20]
head(tmp_data_avg_con)
tmp_data_avg_con

tmp_data_avg_con$maxtpm <- apply(tmp_data_avg_con[1:5],1,max)
head(tmp_data_avg_con)

#Selecting genes >0.5 for treated: 
treated.tpm.data <- Treated.tpm.data[,1:15]
head(treated.tpm.data)

#creating average per timepoint for treated
treated.tpm.data$T0 <- (treated.tpm.data[,1] + treated.tpm.data[,2] + treated.tpm.data[,3])
treated.tpm.data$T1 <- (treated.tpm.data[,4] + treated.tpm.data[,5] + treated.tpm.data[,6])
treated.tpm.data$T2 <- (treated.tpm.data[,7] + treated.tpm.data[,8] + treated.tpm.data[,9])
treated.tpm.data$T3 <- (treated.tpm.data[,10] + treated.tpm.data[,11] + treated.tpm.data[,12])
treated.tpm.data$T4 <- (treated.tpm.data[,13] + treated.tpm.data[,14] + treated.tpm.data[,15])

#rename columns to keep average time point 
colnames(treated.tpm.data)[16:20]
tpm_data_avg_treated <- treated.tpm.data[,16:20]
head(tpm_data_avg_treated)

#Making a new column with the max 
tpm_data_avg_treated$maxtpm <- apply(tpm_data_avg_treated[1:5],1,max)
head(tpm_data_avg_treated)

#select only genes >0.5 max tpm for control and treated 
control_genes_0.5tpm <- tmp_data_avg_con[tmp_data_avg_con$maxtpm > 0.5,]
dim(tmp_data_avg_con)
dim(control_genes_0.5tpm)

Treated_genes_0.5tpm <-tpm_data_avg_treated[tpm_data_avg_treated$maxtpm > 0.5,]
dim(tpm_data_avg_treated)
dim(Treated_genes_0.5tpm)

#merging >0.5tpm of control and treated samples 

merged_genes_0.5tpm <- merge(control_genes_0.5tpm, Treated_genes_0.5tpm, by = 0, all = T)
head(merged_genes_0.5tpm)
dim(merged_genes_0.5tpm)
colnames(merged_genes_0.5tpm)[1] <- "gene"
head(merged_genes_0.5tpm)

#remove low confidence genes using grepl function to pattern match
merged_genes_0.5tpm <- merged_genes_0.5tpm[!grepl("LC",merged_genes_0.5tpm$gene), ]
head(merged_genes_0.5tpm)
dim(merged_genes_0.5tpm)

#read in count data
getwd()
Control.count.data <- read.csv(file=paste0(data_dir,"control_timecourse_count.tsv"), sep = "\t")
head(Control.count.data)
Treated.count.data <- read.csv(file=paste0(data_dir,"Treated_timecourse_count.tsv"), sep = "\t")
head(Treated.count.data)

#filtering counts to match filtered 0.5 tpm merged files
FLB.control.count.data <- Control.count.data[,1:15]
head(FLB.control.count.data)
dim(FLB.control.count.data)
FLB.control.count.data.0.5.tpm <- FLB.control.count.data[rownames(FLB.control.count.data) %in% merged_genes_0.5tpm$gene,]
head(FLB.control.count.data.0.5.tpm)
dim(FLB.control.count.data.0.5.tpm)

FLB.treated.count.data <- Treated.count.data[,1:15]
head(FLB.treated.count.data)
dim (FLB.treated.count.data)
FLB.treated.count.data.0.5tpm <- FLB.treated.count.data[rownames(FLB.treated.count.data) %in% merged_genes_0.5tpm$gene,]
dim(FLB.treated.count.data.0.5tpm)

head(FLB.treated.count.data.0.5tpm)
head(FLB.control.count.data.0.5.tpm)

###DESeq2

counts_for_DESeq2 <- merge(FLB.control.count.data.0.5.tpm, FLB.treated.count.data.0.5tpm, by = 0)
head(counts_for_DESeq2)
rownames(counts_for_DESeq2) <- counts_for_DESeq2[,1]
counts_for_DESeq2 <- counts_for_DESeq2[,-1]
head(counts_for_DESeq2)
dim(counts_for_DESeq2)

library(DESeq2)
colnames(counts_for_DESeq2)
timepoints <- c("T0_WT" , "T1_WT", "T2_WT" , "T3_WT", "T4_WT", "T0_CHA", "T1_CHA", "T2_CHA", "T3_CHA", "T4_CHA")
CondVector <- rep(timepoints, each=3)
CondVector

#column and row names 

SampleTable <- data.frame(row.names=colnames(counts), condition=as.factor(CondVector))
SampleTable

#Rounding counts to specified no of decimal places - default=0
counts <- round(counts_for_DESeq2)
head(counts)

dds<- DESeqDataSetFromMatrix(countData = counts, colData = SampleTable, design=~condition)
dds

#T0_WT treated as reference for relevel
dds$condition <- relevel(dds$condition, "T0_WT")

dim(dds)

#save copy of dds 
dds_copy <- dds

dds <- DESeq(dds)

bckCDS_1 <- dds

#save the results of the DEseq (BckCDS_1)

# timepoints and tissue samples 

timepoints

times_list <- c("0","1","2","3","4")
times_list

#get data for each timepoint for WT Control and CHA Treated

library(DESeq2)

for (i in times_list){
  bck_res <- results(bckCDS_1, contrast = c("condition",paste0("T",i,"_CHA"),paste0("T",i,"_WT")))

#Sorting results on padj (Adjusted p-value)
ordered_res <- bck_res[order(bck_res$padj),]
head(ordered_res)
tail(ordered_res)
# removing incomplete data in ordered file (Padj sorted)

ordered_res_na.rm <- na.omit(ordered_res)
ordered_res_na.rm 
head(ordered_res_na.rm)
tail(ordered_res_na.rm)


#output ordered res to csv 
write.csv(ordered_res_na.rm[ordered_res_na.rm$padj<0.05,],file=paste0(i,"DESeq_CHA_vs_WT_results.csv"))
}

##### Making PCA plot of DESeq data

library(ggplot2)
library(DESeq2)

#Normalisation with rlog - look at clustering 
rld <- rlog(dds)
plotPCA(rld)

#Normalisation with VST - variance stabilizing transformation - gives better plot.
rld_2 <-vst(dds)
plotPCA(rld_2)


#trying to change shape and labels
library(tidyr)
PCA_data <- plotPCA(rld_2, intgroup=c("condition"), returnData=TRUE)
PCA_data <- separate(PCA_data, group, c("time","Treatment"),"_")
head(PCA_data)
PCA_data$time <- factor(PCA_data$time, levels=c("T0","T1","T2","T3","T4"))

#different labelling styles: 

pdf(file ="deseq2_pca_VST_ggplot_no_names.pdf", height = 8, width = 10)
PCA_no_names <- ggplot(PCA_data,aes(PC1,PC2, color=time, group=Treatment))+geom_point(aes(shape=Treatment),size=3)+ scale_color_brewer(palette="Set2")
PCA_no_names

PCA_sample_names <- ggplot(PCA_data,aes(PC1,PC2, color=time, group=Treatment))+geom_point(aes(shape=Treatment),size=3)+ geom_text(aes(label=condition),hjust=-0.3, vjust=0.0, size=2)+ scale_color_brewer(palette="Set2")
PCA_sample_names

PCA_treatment <- ggplot(PCA_data,aes(PC1,PC2, color=time, group=Treatment))+geom_point(aes(shape=Treatment),size=2.5)+ geom_text(aes(label=Treatment),hjust=-0.4, vjust=-0.4, size=2)+ scale_color_brewer(palette="Set2")
PCA_treatment

