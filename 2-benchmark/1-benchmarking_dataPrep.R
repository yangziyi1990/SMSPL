#################################################
# Code version V1.0
# 1-benchmarking_dataPrep.R
# Date: December 26, 2019
# updated: December 26, 2019
# Author: Zi-Yi Yang
#
#################################################
WhereAmI <- "D:/path/1-data"

## load libraries
library(amritr)
library(mixOmics)
library(tidyverse)

# import SNF datasets
## COAD

COAD.group0 <- read.delim(paste0(WhereAmI, "colon/colon_Survival.txt"), row.names = 1)
table(COAD.group0$Death)
COAD.group <- rep(NA, nrow(COAD.group0))
COAD.group[COAD.group0$Survival <= median(COAD.group0$Survival)] <- "low"
COAD.group[COAD.group0$Survival > median(COAD.group0$Survival)] <- "high"
table(COAD.group)
names(COAD.group) <- unlist(lapply(strsplit(rownames(COAD.group0), "-"), function(i) i[3]))

# mRNA
COAD.mrna0 <- read.delim(paste0(WhereAmI, "colon/colon_Gene_Expression.txt"))
colnames(COAD.mrna0) <- c(colnames(COAD.mrna0)[-1], NA)
COAD.mrna <- as.matrix(na.omit(COAD.mrna0[, -ncol(COAD.mrna0)]))
colnames(COAD.mrna) <- unlist(lapply(strsplit(colnames(COAD.mrna), "\\."), function(i) i[3]))
rownames(COAD.mrna) <- paste("mrna", rownames(COAD.mrna), sep = "_")
# miRNA
COAD.mirna0 <- read.delim(paste0(WhereAmI, "colon/colon_Mirna_Expression.txt"))
colnames(COAD.mirna0) <- c(colnames(COAD.mirna0)[-1], NA)
COAD.mirna <- as.matrix(na.omit(COAD.mirna0[, -ncol(COAD.mirna0)]))
colnames(COAD.mirna) <- unlist(lapply(strsplit(colnames(COAD.mirna), "\\."), function(i) i[3]))
rownames(COAD.mirna) <- paste("mirna", rownames(COAD.mirna), sep = "_")
# Methylation
COAD.methylation <- read.delim(paste0(WhereAmI, "colon/colon_Methy_Expression.txt"), row.names = NULL)
colnames(COAD.methylation) <- c(colnames(COAD.methylation)[-2], NA)
rownames(COAD.methylation) <- paste("methylation", COAD.methylation$row.names, 1:nrow(COAD.methylation), sep = "_")
COAD.methylation <- as.matrix(COAD.methylation[, -c(1, ncol(COAD.methylation))])
colnames(COAD.methylation) <- unlist(lapply(strsplit(colnames(COAD.methylation), "\\."), function(i) i[3]))
## check ordering
all(names(COAD.group) == colnames(COAD.mrna))
all(names(COAD.group) == colnames(COAD.mirna))
all(names(COAD.group) == colnames(COAD.methylation))

## check to see if there are no NAs
sum(is.na(COAD.mrna)); sum(is.na(COAD.mirna)); sum(is.na(COAD.methylation));
dim(COAD.mrna)  #  17814    92
dim(COAD.mirna) # 312  92
dim(COAD.methylation)  # 23088    92

### dataDistribution-COAD
pdf(paste0(WhereAmI, "/dataDistribution-COAD.pdf"),width = 8, height = 3)
par(mfrow = c(1, 3))
hist(COAD.mrna,col=rgb(0.3,0.6,1,1));hist(COAD.mirna,col=rgb(0.3,0.6,1,1));hist(COAD.methylation,col=rgb(0.3,0.6,1,1))
dev.off()

## KRCCC

KRCCC.group0 <- read.delim(paste0(WhereAmI, "kidney/kidney_Survival.txt"), row.names = 1)
KRCCC.group <- rep(NA, nrow(KRCCC.group0))
names(KRCCC.group) <- unlist(lapply(strsplit(rownames(KRCCC.group0), "-"), function(i) i[3]))
KRCCC.group[KRCCC.group0$Survival <= median(KRCCC.group0$Survival)] <- "low"
KRCCC.group[KRCCC.group0$Survival > median(KRCCC.group0$Survival)] <- "high"
table(KRCCC.group)

## mRNA
KRCCC.mrna0 <- read.delim(paste0(WhereAmI, "kidney/kidney_Gene_Expression.txt"), row.names = 1) %>% 
  mutate(Gene = unlist(lapply(strsplit(rownames(.), "LLL"), function(i) i[1]))) %>% group_by(Gene) %>% 
  dplyr::summarise_all(funs(mean))
KRCCC.mrna <- as.matrix(KRCCC.mrna0[, -1])
rownames(KRCCC.mrna) <- as.character(KRCCC.mrna0$Gene)
colnames(KRCCC.mrna) <- c(colnames(KRCCC.mrna)[-1], NA)
KRCCC.mrna <- as.matrix(KRCCC.mrna[, colSums(is.na(KRCCC.mrna)) == 0])
colnames(KRCCC.mrna) <- unlist(lapply(strsplit(colnames(KRCCC.mrna), "\\."), function(i) i[3]))
rownames(KRCCC.mrna) <- paste("mrna", rownames(KRCCC.mrna), sep = "_")
# miRNA
KRCCC.mirna0 <- read.delim(paste0(WhereAmI, "kidney/kidney_Mirna_Expression.txt"), row.names = 1)
KRCCC.mirna <- as.matrix(KRCCC.mirna0[, colSums(is.na(KRCCC.mirna0)) == 0])
colnames(KRCCC.mirna) <- unlist(lapply(strsplit(colnames(KRCCC.mirna0)[-1], "\\."), function(i) i[3]))
rownames(KRCCC.mirna) <- paste("mirna", rownames(KRCCC.mirna), sep = "_")
## methylations
KRCCC.methylation <- read.delim(paste0(WhereAmI, "kidney/kidney_Methy_Expression.txt"), row.names = 1)
colnames(KRCCC.methylation) <- unlist(lapply(strsplit(colnames(KRCCC.methylation)[-1], "\\."), function(i) i[3]))
KRCCC.methylation <- as.matrix(KRCCC.methylation[, -ncol(KRCCC.methylation)])
rownames(KRCCC.methylation) <- sapply(strsplit(rownames(KRCCC.methylation), " "), function(i) paste(c("methylation", rev(i)), collapse="_"))

dim(KRCCC.mrna); dim(KRCCC.mirna); dim(KRCCC.methylation); 
table(KRCCC.group); length(KRCCC.group)
all(names(KRCCC.group) == colnames(KRCCC.mrna))
all(names(KRCCC.group) == colnames(KRCCC.mirna))
all(names(KRCCC.group) == colnames(KRCCC.methylation))

## check to see if there are no NAs
sum(is.na(KRCCC.mrna)); sum(is.na(KRCCC.mirna)); sum(is.na(KRCCC.methylation));
dim(KRCCC.mrna) # 17665   122
dim(KRCCC.mirna) # 329 122
dim(KRCCC.methylation) # 24960   122


### dataDistribution-KRCCC
pdf(paste0(WhereAmI, "/dataDistribution-KRCCC.pdf"),width = 8, height = 3)
par(mfrow = c(1, 3))
hist(KRCCC.mrna,col=rgb(0.3,0.6,1,1)); hist(KRCCC.mirna,col=rgb(0.3,0.6,1,1)); hist(KRCCC.methylation,col=rgb(0.3,0.6,1,1))
dev.off()


## GBM

GBM.group0 <- read.delim(paste0(WhereAmI, "GBM/GLIO_Survival.txt"))
GBM.group <- rep(NA, nrow(GBM.group0))
names(GBM.group) <- gsub("-", ".", as.character(GBM.group0$PatientID))
names(GBM.group) <- unlist(lapply(strsplit(as.character(GBM.group0$PatientID), "-"), function(i) i[3]))
GBM.group[GBM.group0$Survival <= median(GBM.group0$Survival)] <- "low"
GBM.group[GBM.group0$Survival > median(GBM.group0$Survival)] <- "high"
table(GBM.group)

## remove doubles
which(names(GBM.group) %in% names(table(names(GBM.group)))[table(names(GBM.group)) > 1])
names(GBM.group)[which(names(GBM.group) %in% names(table(names(GBM.group)))[table(names(GBM.group)) > 1])]
-c(52,80)

## mRNA
GBM.mrna0 <- read.delim(paste0(WhereAmI, "GBM/GLIO_Gene_Expression.txt"), row.names = 1)
colnames(GBM.mrna0) <- c(colnames(GBM.mrna0)[-1], NA)
GBM.mrna <- as.matrix(na.omit(GBM.mrna0[, -ncol(GBM.mrna0)]))
colnames(GBM.mrna) <- unlist(lapply(strsplit(colnames(GBM.mrna), "\\."), function(i) i[3]))
rownames(GBM.mrna) <- paste("mrna", rownames(GBM.mrna), sep = "_")
# miRNA
GBM.mirna0 <- read.delim(paste0(WhereAmI, "GBM/GLIO_Mirna_Expression.txt"))
colnames(GBM.mirna0) <- c(colnames(GBM.mirna0)[-1], NA)
GBM.mirna <- as.matrix(na.omit(GBM.mirna0[, -ncol(GBM.mirna0)]))
colnames(GBM.mirna) <- unlist(lapply(strsplit(colnames(GBM.mirna), "\\."), function(i) i[3]))
rownames(GBM.mirna) <- paste("mirna", rownames(GBM.mirna), sep = "_")
# Methylation
GBM.methylation <- read.delim(paste0(WhereAmI, "GBM/GLIO_Methy_Expression.txt"), row.names = 1)
colnames(GBM.methylation) <- c(colnames(GBM.methylation)[-1], NA)
GBM.methylation <- as.matrix(GBM.methylation[,-ncol(GBM.methylation)])
colnames(GBM.methylation) <- unlist(lapply(strsplit(colnames(GBM.methylation), "\\."), function(i) i[3]))
rownames(GBM.methylation) <- paste("methylation", rownames(GBM.methylation), sep="_")


## remove duplicates samples
GBM.group <- GBM.group[-c(52,80)]
GBM.mrna <- GBM.mrna[, -c(52,80)]
GBM.mirna <- GBM.mirna[, -c(52,80)]
GBM.methylation <- GBM.methylation[, -c(52,80)]

all(names(GBM.group) == colnames(GBM.mrna))
all(names(GBM.group) == colnames(GBM.mirna))
all(names(GBM.group) == colnames(GBM.methylation))

## check to see if there are no NAs
sum(is.na(GBM.mrna)); sum(is.na(GBM.mirna)); sum(is.na(GBM.methylation));
dim(GBM.mrna) # 12042   213
dim(GBM.mirna) # 534 213
dim(GBM.methylation) # 1305  213

### dataDistribution-GBM

pdf(paste0(WhereAmI, "/dataDistribution-GBM.pdf"),width = 8, height = 3)
par(mfrow = c(1, 3))
hist(GBM.mrna,col=rgb(0.3,0.6,1,1)); hist(GBM.mirna,col=rgb(0.3,0.6,1,1)); hist(GBM.methylation,col=rgb(0.3,0.6,1,1))
dev.off()

# LSCC

LSCC.group0 <- read.delim(paste0(WhereAmI, "Lung/LUNG_Survival.txt"), row.names = 1)
LSCC.group <- rep(NA, nrow(LSCC.group0))
LSCC.group[LSCC.group0$Survival <= median(LSCC.group0$Survival)] <- "low"
LSCC.group[LSCC.group0$Survival > median(LSCC.group0$Survival)] <- "high"
table(LSCC.group)
names(LSCC.group) <- unlist(lapply(strsplit(rownames(LSCC.group0), "-"), function(i) i[3]))

## mRNA
LSCC.mrna0 <- read.delim(paste0(WhereAmI, "Lung/LUNG_Gene_Expression.txt"), row.names = 1)
LSCC.mrna <- as.matrix(LSCC.mrna0[, colSums(is.na(LSCC.mrna0)) == 0])
colnames(LSCC.mrna) <- unlist(lapply(strsplit(colnames(LSCC.mrna0)[-1], "\\."), function(i) i[3]))
rownames(LSCC.mrna) <- paste("mrna", rownames(LSCC.mrna), sep = "_")
# miRNA
LSCC.mirna0 <- read.delim(paste0(WhereAmI, "Lung/LUNG_Mirna_Expression.txt"), row.names = 1)
LSCC.mirna <- as.matrix(LSCC.mirna0[, colSums(is.na(LSCC.mirna0)) == 0])
colnames(LSCC.mirna) <- unlist(lapply(strsplit(colnames(LSCC.mirna0)[-1], "\\."), function(i) i[3]))
rownames(LSCC.mirna) <- paste("mirna", rownames(LSCC.mirna), sep = "_")
## methylations
LSCC.methylation0 <- read.delim(paste0(WhereAmI, "Lung/LUNG_Methy_Expression.txt"), row.names = 1)
LSCC.methylation <- as.matrix(LSCC.methylation0[, colSums(is.na(LSCC.methylation0)) == 0])
colnames(LSCC.methylation) <- unlist(lapply(strsplit(colnames(LSCC.methylation0)[-1], "\\."), function(i) i[3]))

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
rownames(LSCC.methylation) <- paste("methylation", unlist(lapply(ann[rownames(LSCC.methylation), "UCSC_RefGene_Name"], function(i){
  paste(unique(unlist(strsplit(as.character(i), ";"))), collapse = ";")
})), rownames(LSCC.methylation), sep = "_")
nrow(LSCC.methylation)

dim(LSCC.mrna); dim(LSCC.mirna); dim(LSCC.methylation); 
table(LSCC.group); length(LSCC.group)
all(names(LSCC.group) == colnames(LSCC.mrna))
all(names(LSCC.group) == colnames(LSCC.mirna))
all(names(LSCC.group) == colnames(LSCC.methylation))

## check to see if there are no NAs
sum(is.na(LSCC.mrna)); sum(is.na(LSCC.mirna)); sum(is.na(LSCC.methylation));
dim(LSCC.mrna)  # 12042   106
dim(LSCC.mirna) # 352 106
dim(LSCC.methylation) # 23074   106

### dataDistribution-LSCC
pdf(paste0(WhereAmI, "/dataDistribution-LSCC.pdf"),width = 8, height = 3)
par(mfrow = c(1, 3))
hist(LSCC.mrna,col=rgb(0.3,0.6,1,1)); hist(LSCC.mirna,col=rgb(0.3,0.6,1,1)); hist(LSCC.methylation,col=rgb(0.3,0.6,1,1))
dev.off()



## Compare survival times

survivalDat <- data.frame(survival_times=c(COAD.group0$Survival,KRCCC.group0$Survival,GBM.group0$Survival[-c(52,80)],LSCC.group0$Survival),
                          group=factor(rep(c("COAD","KRCCC","GBM","LSCC"), c(length(COAD.group), length(KRCCC.group), length(GBM.group), length(LSCC.group))), levels=c("KRCCC","LSCC","GBM","COAD")),
                          survival=c(COAD.group, KRCCC.group, GBM.group, LSCC.group))

survivalDat_summary <- survivalDat %>% 
  group_by(group, survival) %>% 
  summarise(survival_times = 4000, n = paste("n", n(), sep="="))
survivalDat_summary$survival_times[seq(2,nrow(survivalDat_summary),2)] <- -100

pdf(paste0(WhereAmI, "label_generation.pdf"), width = 6, height = 4)
ggplot(survivalDat, aes(x=group, y=survival_times, color=survival, fill=survival)) +
  geom_boxplot() +
  theme_bw() +
  customTheme(sizeStripFont = 15, xAngle = 0, hjust = 0.5, vjust = 0.5, 
              xSize = 10, ySize = 10, xAxisSize = 10, yAxisSize = 10) +
  ylab("Survival times (days)") +
  xlab("Cancer datasets") +
  geom_text(data = survivalDat_summary, aes(x=group, y=survival_times, label = n), position=position_dodge(width=0.7)) +
  scale_fill_brewer(palette = "Pastel1")
dev.off()


## save datasets to rdata file

snf_data = list(COAD = list(mrna = t(COAD.mrna), mirna = t(COAD.mirna), methylation = t(COAD.methylation)),
                KRCCC = list(mrna = t(KRCCC.mrna), mirna = t(KRCCC.mirna), methylation = t(KRCCC.methylation)),
                GBM = list(mrna = t(GBM.mrna), mirna = t(GBM.mirna), methylation = t(GBM.methylation)),
                LSCC = list(mrna = t(LSCC.mrna), mirna = t(LSCC.mirna), methylation = t(LSCC.methylation)))
snf_group = list(COAD = COAD.group,
                 KRCCC = KRCCC.group,
                 GBM = GBM.group,
                 LSCC = LSCC.group)

save(snf_data=snf_data, snf_group=snf_group,
     file = paste0(WhereAmI, "SNFdatasets.RDATA"))

