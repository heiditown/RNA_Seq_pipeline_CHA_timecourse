# Graphical display of Threshold table - showing upreg/downreg DE genes at different strigency p-adj values. 

library(readr)
library(readxl)
library(ggplot2)

getwd()

threshold_data<- read.csv("threshold_data.xlsx.csv")
threshold_data

threshold_data2 <- read.csv("threshold_data2.csv")
threshold_data2

barchart1 <- ggplot(threshold_data2, aes(fill=Threshold, y=Number.of.genes, x=Time.Point))+ geom_bar(position="dodge", stat="identity")

Timepoints <-(threshold_data2$Time.Point)
No.of.genes <- (threshold_data2$Number.of.genes)
threshold_data2$factor_Threshold <- as.factor(threshold_data2$Threshold)
threshold_data2$Up.Down.regulated_factor <-as.factor(threshold_data2$Up.Down.regulated)
Regulation_effect <- (threshold_data2$Up.Down.regulated_factor)

Regulation_effect
p_value_threshold <-(threshold_data2$factor_Threshold)

data <- data.frame(Timepoints,No.of.genes,Regulation_effect,p_value_threshold)

barchart1 <- ggplot(data, aes(fill=p_value_threshold, y=No.of.genes, x=Timepoints))+ geom_bar(position="dodge", stat="identity")


All_genes_DE <- barchart1 + labs(x="Time points", y="Number of genes DE")
AG_DE_2 <- All_genes_DE + labs(fill="p-Value Thresholds")

threshold_data3 <- read.csv("threshold data 3.csv")

D3 <- ggplot(threshold_data3, aes(fill=Thresholds, y=Number.of.genes, x=Time.Point))+ geom_bar(position="dodge", stat="identity")
d3
D3 + labs(x="Time points", y="Number of genes")


one <- ggplot(threshold_data2)+ geom_bar(aes(x=Time.Point, y=Number.of.genes, fill=Threshold), position="dodge", stat='identity')
one
two <- ggplot(threshold_data2)+
  geom_bar(aes(x=Time.Point, y=Number.of.genes, fill=Up.Down.regulated_factor), position = "stack", stat="identity")
gridExtra::grid.arrange(one,two, nrow=2)

split <-ggplot(threshold_data2)+
  geom_bar(aes(x=Time.Point, y=Number.of.genes, fill=Up.Down.regulated_factor), position ="stack", stat="identity")+
  facet_wrap(~Threshold)
split + labs(x="Timepoints", y="Number of genes", fill="DE Regulation Effect")
