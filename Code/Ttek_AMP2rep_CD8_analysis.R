library(ggplot2)
library(tidyr)
library(dplyr)
library(Seurat)
library(scRepertoire)
library(harmony)
library(RColorBrewer)
library(EnhancedVolcano)
library(tibble)
library(patchwork)
library(cowplot)
library(ggpubr)
library(circlize)
library(slingshot)
library(viridis)


# Initial loading of TCR/RNA/ADT data -----------------------------------------------------

#loading concatenated QC metrics from 10x output
phase2_TCR_QC <- read.csv(file = "AMP_Phase_2_TCRmetricsQC.csv")

#loading TCR files
P_300_0150_TCR <- read.csv("data/Sample_300_0150_PBL_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
S_300_0150_TCR <- read.csv("data/Sample_300_0150_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

P_300_0171_1_TCR <- read.csv("data/Sample_300_0171_PBL_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
P_300_0171_2_TCR <- read.csv("data/Sample_300_0171_PBL_BT_2/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
S_300_0171_TCR <- read.csv("data/Sample_300_0171_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

P_300_0173_TCR <- read.csv("data/Sample_300_0173_PBL_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
S_300_0173_TCR <- read.csv("data/Sample_300_0173_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

P_300_0174_TCR <- read.csv("data/Sample_300_0174_PBL_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
S_300_0174_TCR <- read.csv("data/Sample_300_0174_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

P_300_0392_TCR <- read.csv("data/Sample_300_0392_PBL_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
S_300_0392_TCR <- read.csv("data/Sample_300_0392_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

P_300_0410_TCR <- read.csv("data/Sample_300_0410_PBL_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
S_300_0410_TCR <- read.csv("data/Sample_300_0410_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

P_300_0414_TCR <- read.csv("data/Sample_300_0414_PBL_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
S_300_0414_TCR <- read.csv("data/Sample_300_0414_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

S_300_0416_TCR <- read.csv("data/Sample_300_0416_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

S_300_1883_TCR <- read.csv("data/Sample_300_1883_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

P_300_1930_TCR <- read.csv("data/Sample_300_1930_PBL_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
S_300_1930_TCR <- read.csv("data/Sample_300_1930_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

P_301_0174_TCR <- read.csv("data/Sample_301_0174_PBL_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
S_301_0174_TCR <- read.csv("data/Sample_301_0174_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

P_301_0270_TCR <- read.csv("data/Sample_301_0270_PBL_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
S_301_0270_TCR <- read.csv("data/Sample_301_0270_Syn_BT/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)

TCR_list <- list(P_300_0150_TCR, S_300_0150_TCR, P_300_0171_1_TCR, P_300_0171_2_TCR, S_300_0171_TCR,
                 P_300_0173_TCR, S_300_0173_TCR, P_300_0174_TCR, S_300_0174_TCR, P_300_0392_TCR, S_300_0392_TCR,
                 P_300_0410_TCR, S_300_0410_TCR, P_300_0414_TCR, S_300_0414_TCR, S_300_0416_TCR, S_300_1883_TCR,
                 P_300_1930_TCR, S_300_1930_TCR, P_301_0174_TCR, S_301_0174_TCR, P_301_0270_TCR, S_301_0270_TCR)

AMP_2_RA_TCR <- combineTCR(TCR_list, samples = phase2_TCR_QC$Sample, ID = phase2_TCR_QC$Patient, cells = "T-AB")



#loading RNA and ADT

#PBL_300_0150
PBL_300_0150_GEX <- Read10X_h5("data/Sample_300_0150_PBL_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_300_0150_genedata <- PBL_300_0150_GEX$`Gene Expression`
colnames(PBL_300_0150_genedata) <- paste0("PBL_300_0150_", colnames(PBL_300_0150_genedata))
PBL_300_0150_protdata <- PBL_300_0150_GEX$`Antibody Capture`
colnames(PBL_300_0150_protdata) <- paste0("PBL_300_0150_", colnames(PBL_300_0150_protdata))

PBL_300_0150 <- CreateSeuratObject(PBL_300_0150_genedata)
PBL_300_0150[["ADT"]] <- CreateAssayObject(PBL_300_0150_protdata)
PBL_300_0150$Patient <- "300_0150"
PBL_300_0150$Sample <- "PBL"

PBL_300_0150 <- combineExpression(AMP_2_RA_TCR$PBL_300_0150, PBL_300_0150, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_300_0150$cloneProp <- PBL_300_0150$Frequency/sum(table(PBL_300_0150$CTaa))

#SYN_300_0150
SYN_300_0150_GEX <- Read10X_h5("data/Sample_300_0150_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_300_0150_genedata <- SYN_300_0150_GEX$`Gene Expression`
colnames(SYN_300_0150_genedata) <- paste0("SYN_300_0150_", colnames(SYN_300_0150_genedata))
SYN_300_0150_protdata <- SYN_300_0150_GEX$`Antibody Capture`
colnames(SYN_300_0150_protdata) <- paste0("SYN_300_0150_", colnames(SYN_300_0150_protdata))

SYN_300_0150 <- CreateSeuratObject(SYN_300_0150_genedata)
SYN_300_0150[["ADT"]] <- CreateAssayObject(SYN_300_0150_protdata)
SYN_300_0150$Patient <- "300_0150"
SYN_300_0150$Sample <- "SYN"

SYN_300_0150 <- combineExpression(AMP_2_RA_TCR$SYN_300_0150, SYN_300_0150, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_300_0150$cloneProp <- SYN_300_0150$Frequency/sum(table(SYN_300_0150$CTaa))

#PBL_300_0171_1
PBL_300_0171_GEX <- Read10X_h5("data/Sample_300_0171_PBL_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_300_0171_genedata <- PBL_300_0171_GEX$`Gene Expression`
colnames(PBL_300_0171_genedata) <- paste0("PBL_300_0171_", colnames(PBL_300_0171_genedata))
PBL_300_0171_protdata <- PBL_300_0171_GEX$`Antibody Capture`
colnames(PBL_300_0171_protdata) <- paste0("PBL_300_0171_", colnames(PBL_300_0171_protdata))

PBL_300_0171 <- CreateSeuratObject(PBL_300_0171_genedata)
PBL_300_0171[["ADT"]] <- CreateAssayObject(PBL_300_0171_protdata)
PBL_300_0171$Patient <- "300_0171"
PBL_300_0171$Sample <- "PBL"

PBL_300_0171 <- combineExpression(AMP_2_RA_TCR$PBL_300_0171, PBL_300_0171, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_300_0171$cloneProp <- PBL_300_0171$Frequency/sum(table(PBL_300_0171$CTaa))

#PBL_300_0171_2
PBL_300_0171_2_GEX <- Read10X_h5("data/Sample_300_0171_PBL_BT_2/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_300_0171_2_genedata <- PBL_300_0171_2_GEX$`Gene Expression`
colnames(PBL_300_0171_2_genedata) <- paste0("PBL_300_0171_2_", colnames(PBL_300_0171_2_genedata))
PBL_300_0171_2_protdata <- PBL_300_0171_2_GEX$`Antibody Capture`
colnames(PBL_300_0171_2_protdata) <- paste0("PBL_300_0171_2_", colnames(PBL_300_0171_2_protdata))

PBL_300_0171_2 <- CreateSeuratObject(PBL_300_0171_2_genedata)
PBL_300_0171_2[["ADT"]] <- CreateAssayObject(PBL_300_0171_2_protdata)
PBL_300_0171_2$Patient <- "300_0171"
PBL_300_0171_2$Sample <- "PBL"

PBL_300_0171_2 <- combineExpression(AMP_2_RA_TCR$PBL_300_0171, PBL_300_0171_2, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_300_0171_2$cloneProp <- PBL_300_0171_2$Frequency/sum(table(PBL_300_0171_2$CTaa))

#SYN_300_0171
SYN_300_0171_GEX <- Read10X_h5("data/Sample_300_0171_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_300_0171_genedata <- SYN_300_0171_GEX$`Gene Expression`
colnames(SYN_300_0171_genedata) <- paste0("SYN_300_0171_", colnames(SYN_300_0171_genedata))
SYN_300_0171_protdata <- SYN_300_0171_GEX$`Antibody Capture`
colnames(SYN_300_0171_protdata) <- paste0("SYN_300_0171_", colnames(SYN_300_0171_protdata))

SYN_300_0171 <- CreateSeuratObject(SYN_300_0171_genedata)
SYN_300_0171[["ADT"]] <- CreateAssayObject(SYN_300_0171_protdata)
SYN_300_0171$Patient <- "300_0171"
SYN_300_0171$Sample <- "SYN"

SYN_300_0171 <- combineExpression(AMP_2_RA_TCR$SYN_300_0171, SYN_300_0171, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_300_0171$cloneProp <- SYN_300_0171$Frequency/sum(table(SYN_300_0171$CTaa))

#PBL_300_0173
PBL_300_0173_GEX <- Read10X_h5("data/Sample_300_0173_PBL_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_300_0173_genedata <- PBL_300_0173_GEX$`Gene Expression`
colnames(PBL_300_0173_genedata) <- paste0("PBL_300_0173_", colnames(PBL_300_0173_genedata))
PBL_300_0173_protdata <- PBL_300_0173_GEX$`Antibody Capture`
colnames(PBL_300_0173_protdata) <- paste0("PBL_300_0173_", colnames(PBL_300_0173_protdata))

PBL_300_0173 <- CreateSeuratObject(PBL_300_0173_genedata)
PBL_300_0173[["ADT"]] <- CreateAssayObject(PBL_300_0173_protdata)
PBL_300_0173$Patient <- "300_0173"
PBL_300_0173$Sample <- "PBL"

PBL_300_0173 <- combineExpression(AMP_2_RA_TCR$PBL_300_0173, PBL_300_0173, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_300_0173$cloneProp <- PBL_300_0173$Frequency/sum(table(PBL_300_0173$CTaa))

#SYN_300_0173
SYN_300_0173_GEX <- Read10X_h5("data/Sample_300_0173_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_300_0173_genedata <- SYN_300_0173_GEX$`Gene Expression`
colnames(SYN_300_0173_genedata) <- paste0("SYN_300_0173_", colnames(SYN_300_0173_genedata))
SYN_300_0173_protdata <- SYN_300_0173_GEX$`Antibody Capture`
colnames(SYN_300_0173_protdata) <- paste0("SYN_300_0173_", colnames(SYN_300_0173_protdata))

SYN_300_0173 <- CreateSeuratObject(SYN_300_0173_genedata)
SYN_300_0173[["ADT"]] <- CreateAssayObject(SYN_300_0173_protdata)
SYN_300_0173$Patient <- "300_0173"
SYN_300_0173$Sample <- "SYN"

SYN_300_0173 <- combineExpression(AMP_2_RA_TCR$SYN_300_0173, SYN_300_0173, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_300_0173$cloneProp <- SYN_300_0173$Frequency/sum(table(SYN_300_0173$CTaa))

#PBL_300_0174
PBL_300_0174_GEX <- Read10X_h5("data/Sample_300_0174_PBL_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_300_0174_genedata <- PBL_300_0174_GEX$`Gene Expression`
colnames(PBL_300_0174_genedata) <- paste0("PBL_300_0174_", colnames(PBL_300_0174_genedata))
PBL_300_0174_protdata <- PBL_300_0174_GEX$`Antibody Capture`
colnames(PBL_300_0174_protdata) <- paste0("PBL_300_0174_", colnames(PBL_300_0174_protdata))

PBL_300_0174 <- CreateSeuratObject(PBL_300_0174_genedata)
PBL_300_0174[["ADT"]] <- CreateAssayObject(PBL_300_0174_protdata)
PBL_300_0174$Patient <- "300_0174"
PBL_300_0174$Sample <- "PBL"

PBL_300_0174 <- combineExpression(AMP_2_RA_TCR$PBL_300_0174, PBL_300_0174, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_300_0174$cloneProp <- PBL_300_0174$Frequency/sum(table(PBL_300_0174$CTaa))

#SYN_300_0174
SYN_300_0174_GEX <- Read10X_h5("data/Sample_300_0174_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_300_0174_genedata <- SYN_300_0174_GEX$`Gene Expression`
colnames(SYN_300_0174_genedata) <- paste0("SYN_300_0174_", colnames(SYN_300_0174_genedata))
SYN_300_0174_protdata <- SYN_300_0174_GEX$`Antibody Capture`
colnames(SYN_300_0174_protdata) <- paste0("SYN_300_0174_", colnames(SYN_300_0174_protdata))

SYN_300_0174 <- CreateSeuratObject(SYN_300_0174_genedata)
SYN_300_0174[["ADT"]] <- CreateAssayObject(SYN_300_0174_protdata)
SYN_300_0174$Patient <- "300_0174"
SYN_300_0174$Sample <- "SYN"

SYN_300_0174 <- combineExpression(AMP_2_RA_TCR$SYN_300_0174, SYN_300_0174, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_300_0174$cloneProp <- SYN_300_0174$Frequency/sum(table(SYN_300_0174$CTaa))

#PBL_300_0392
PBL_300_0392_GEX <- Read10X_h5("data/Sample_300_0392_PBL_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_300_0392_genedata <- PBL_300_0392_GEX$`Gene Expression`
colnames(PBL_300_0392_genedata) <- paste0("PBL_300_0392_", colnames(PBL_300_0392_genedata))
PBL_300_0392_protdata <- PBL_300_0392_GEX$`Antibody Capture`
colnames(PBL_300_0392_protdata) <- paste0("PBL_300_0392_", colnames(PBL_300_0392_protdata))

PBL_300_0392 <- CreateSeuratObject(PBL_300_0392_genedata)
PBL_300_0392[["ADT"]] <- CreateAssayObject(PBL_300_0392_protdata)
PBL_300_0392$Patient <- "300_0392"
PBL_300_0392$Sample <- "PBL"

PBL_300_0392 <- combineExpression(AMP_2_RA_TCR$PBL_300_0392, PBL_300_0392, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_300_0392$cloneProp <- PBL_300_0392$Frequency/sum(table(PBL_300_0392$CTaa))

#SYN_300_0392
SYN_300_0392_GEX <- Read10X_h5("data/Sample_300_0392_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_300_0392_genedata <- SYN_300_0392_GEX$`Gene Expression`
colnames(SYN_300_0392_genedata) <- paste0("SYN_300_0392_", colnames(SYN_300_0392_genedata))
SYN_300_0392_protdata <- SYN_300_0392_GEX$`Antibody Capture`
colnames(SYN_300_0392_protdata) <- paste0("SYN_300_0392_", colnames(SYN_300_0392_protdata))

SYN_300_0392 <- CreateSeuratObject(SYN_300_0392_genedata)
SYN_300_0392[["ADT"]] <- CreateAssayObject(SYN_300_0392_protdata)
SYN_300_0392$Patient <- "300_0392"
SYN_300_0392$Sample <- "SYN"

SYN_300_0392 <- combineExpression(AMP_2_RA_TCR$SYN_300_0392, SYN_300_0392, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_300_0392$cloneProp <- SYN_300_0392$Frequency/sum(table(SYN_300_0392$CTaa))

#PBL_300_0410
PBL_300_0410_GEX <- Read10X_h5("data/Sample_300_0410_PBL_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_300_0410_genedata <- PBL_300_0410_GEX$`Gene Expression`
colnames(PBL_300_0410_genedata) <- paste0("PBL_300_0410_", colnames(PBL_300_0410_genedata))
PBL_300_0410_protdata <- PBL_300_0410_GEX$`Antibody Capture`
colnames(PBL_300_0410_protdata) <- paste0("PBL_300_0410_", colnames(PBL_300_0410_protdata))

PBL_300_0410 <- CreateSeuratObject(PBL_300_0410_genedata)
PBL_300_0410[["ADT"]] <- CreateAssayObject(PBL_300_0410_protdata)
PBL_300_0410$Patient <- "300_0410"
PBL_300_0410$Sample <- "PBL"

PBL_300_0410 <- combineExpression(AMP_2_RA_TCR$PBL_300_0410, PBL_300_0410, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_300_0410$cloneProp <- PBL_300_0410$Frequency/sum(table(PBL_300_0410$CTaa))

#SYN_300_0410
SYN_300_0410_GEX <- Read10X_h5("data/Sample_300_0410_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_300_0410_genedata <- SYN_300_0410_GEX$`Gene Expression`
colnames(SYN_300_0410_genedata) <- paste0("SYN_300_0410_", colnames(SYN_300_0410_genedata))
SYN_300_0410_protdata <- SYN_300_0410_GEX$`Antibody Capture`
colnames(SYN_300_0410_protdata) <- paste0("SYN_300_0410_", colnames(SYN_300_0410_protdata))

SYN_300_0410 <- CreateSeuratObject(SYN_300_0410_genedata)
SYN_300_0410[["ADT"]] <- CreateAssayObject(SYN_300_0410_protdata)
SYN_300_0410$Patient <- "300_0410"
SYN_300_0410$Sample <- "SYN"

SYN_300_0410 <- combineExpression(AMP_2_RA_TCR$SYN_300_0410, SYN_300_0410, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_300_0410$cloneProp <- SYN_300_0410$Frequency/sum(table(SYN_300_0410$CTaa))

#PBL_300_0414
PBL_300_0414_GEX <- Read10X_h5("data/Sample_300_0414_PBL_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_300_0414_genedata <- PBL_300_0414_GEX$`Gene Expression`
colnames(PBL_300_0414_genedata) <- paste0("PBL_300_0414_", colnames(PBL_300_0414_genedata))
PBL_300_0414_protdata <- PBL_300_0414_GEX$`Antibody Capture`
colnames(PBL_300_0414_protdata) <- paste0("PBL_300_0414_", colnames(PBL_300_0414_protdata))

PBL_300_0414 <- CreateSeuratObject(PBL_300_0414_genedata)
PBL_300_0414[["ADT"]] <- CreateAssayObject(PBL_300_0414_protdata)
PBL_300_0414$Patient <- "300_0414"
PBL_300_0414$Sample <- "PBL"

PBL_300_0414 <- combineExpression(AMP_2_RA_TCR$PBL_300_0414, PBL_300_0414, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_300_0414$cloneProp <- PBL_300_0414$Frequency/sum(table(PBL_300_0414$CTaa))

#SYN_300_0414
SYN_300_0414_GEX <- Read10X_h5("data/Sample_300_0414_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_300_0414_genedata <- SYN_300_0414_GEX$`Gene Expression`
colnames(SYN_300_0414_genedata) <- paste0("SYN_300_0414_", colnames(SYN_300_0414_genedata))
SYN_300_0414_protdata <- SYN_300_0414_GEX$`Antibody Capture`
colnames(SYN_300_0414_protdata) <- paste0("SYN_300_0414_", colnames(SYN_300_0414_protdata))

SYN_300_0414 <- CreateSeuratObject(SYN_300_0414_genedata)
SYN_300_0414[["ADT"]] <- CreateAssayObject(SYN_300_0414_protdata)
SYN_300_0414$Patient <- "300_0414"
SYN_300_0414$Sample <- "SYN"

SYN_300_0414 <- combineExpression(AMP_2_RA_TCR$SYN_300_0414, SYN_300_0414, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_300_0414$cloneProp <- SYN_300_0414$Frequency/sum(table(SYN_300_0414$CTaa))

#SYN_300_0416
SYN_300_0416_GEX <- Read10X_h5("data/Sample_300_0416_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_300_0416_genedata <- SYN_300_0416_GEX$`Gene Expression`
colnames(SYN_300_0416_genedata) <- paste0("SYN_300_0416_", colnames(SYN_300_0416_genedata))
SYN_300_0416_protdata <- SYN_300_0416_GEX$`Antibody Capture`
colnames(SYN_300_0416_protdata) <- paste0("SYN_300_0416_", colnames(SYN_300_0416_protdata))

SYN_300_0416 <- CreateSeuratObject(SYN_300_0416_genedata)
SYN_300_0416[["ADT"]] <- CreateAssayObject(SYN_300_0416_protdata)
SYN_300_0416$Patient <- "300_0416"
SYN_300_0416$Sample <- "SYN"

SYN_300_0416 <- combineExpression(AMP_2_RA_TCR$SYN_300_0416, SYN_300_0416, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_300_0416$cloneProp <- SYN_300_0416$Frequency/sum(table(SYN_300_0416$CTaa))

#SYN_300_1883
SYN_300_1883_GEX <- Read10X_h5("data/Sample_300_1883_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_300_1883_genedata <- SYN_300_1883_GEX$`Gene Expression`
colnames(SYN_300_1883_genedata) <- paste0("SYN_300_1883_", colnames(SYN_300_1883_genedata))
SYN_300_1883_protdata <- SYN_300_1883_GEX$`Antibody Capture`
colnames(SYN_300_1883_protdata) <- paste0("SYN_300_1883_", colnames(SYN_300_1883_protdata))

SYN_300_1883 <- CreateSeuratObject(SYN_300_1883_genedata)
SYN_300_1883[["ADT"]] <- CreateAssayObject(SYN_300_1883_protdata)
SYN_300_1883$Patient <- "300_1883"
SYN_300_1883$Sample <- "SYN"
SYN_300_1883 <- combineExpression(AMP_2_RA_TCR$SYN_300_1883, SYN_300_1883, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_300_1883$cloneProp <- SYN_300_1883$Frequency/sum(table(SYN_300_1883$CTaa))

#PBL_300_1930
PBL_300_1930_GEX <- Read10X_h5("data/Sample_300_1930_PBL_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_300_1930_genedata <- PBL_300_1930_GEX$`Gene Expression`
colnames(PBL_300_1930_genedata) <- paste0("PBL_300_1930_", colnames(PBL_300_1930_genedata))
PBL_300_1930_protdata <- PBL_300_1930_GEX$`Antibody Capture`
colnames(PBL_300_1930_protdata) <- paste0("PBL_300_1930_", colnames(PBL_300_1930_protdata))

PBL_300_1930 <- CreateSeuratObject(PBL_300_1930_genedata)
PBL_300_1930[["ADT"]] <- CreateAssayObject(PBL_300_1930_protdata)
PBL_300_1930$Patient <- "300_1930"
PBL_300_1930$Sample <- "PBL"

PBL_300_1930 <- combineExpression(AMP_2_RA_TCR$PBL_300_1930, PBL_300_1930, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_300_1930$cloneProp <- PBL_300_1930$Frequency/sum(table(PBL_300_1930$CTaa))

#SYN_300_1930
SYN_300_1930_GEX <- Read10X_h5("data/Sample_300_1930_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_300_1930_genedata <- SYN_300_1930_GEX$`Gene Expression`
colnames(SYN_300_1930_genedata) <- paste0("SYN_300_1930_", colnames(SYN_300_1930_genedata))
SYN_300_1930_protdata <- SYN_300_1930_GEX$`Antibody Capture`
colnames(SYN_300_1930_protdata) <- paste0("SYN_300_1930_", colnames(SYN_300_1930_protdata))

SYN_300_1930 <- CreateSeuratObject(SYN_300_1930_genedata)
SYN_300_1930[["ADT"]] <- CreateAssayObject(SYN_300_1930_protdata)
SYN_300_1930$Patient <- "300_1930"
SYN_300_1930$Sample <- "SYN"

SYN_300_1930 <- combineExpression(AMP_2_RA_TCR$SYN_300_1930, SYN_300_1930, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_300_1930$cloneProp <- SYN_300_1930$Frequency/sum(table(SYN_300_1930$CTaa))

#PBL_301_0174
PBL_301_0174_GEX <- Read10X_h5("data/Sample_301_0174_PBL_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_301_0174_genedata <- PBL_301_0174_GEX$`Gene Expression`
colnames(PBL_301_0174_genedata) <- paste0("PBL_301_0174_", colnames(PBL_301_0174_genedata))
PBL_301_0174_protdata <- PBL_301_0174_GEX$`Antibody Capture`
colnames(PBL_301_0174_protdata) <- paste0("PBL_301_0174_", colnames(PBL_301_0174_protdata))

PBL_301_0174 <- CreateSeuratObject(PBL_301_0174_genedata)
PBL_301_0174[["ADT"]] <- CreateAssayObject(PBL_301_0174_protdata)
PBL_301_0174$Patient <- "301_0174"
PBL_301_0174$Sample <- "PBL"

PBL_301_0174 <- combineExpression(AMP_2_RA_TCR$PBL_301_0174, PBL_301_0174, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_301_0174$cloneProp <- PBL_301_0174$Frequency/sum(table(PBL_301_0174$CTaa))

#SYN_301_0174
SYN_301_0174_GEX <- Read10X_h5("data/Sample_301_0174_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_301_0174_genedata <- SYN_301_0174_GEX$`Gene Expression`
colnames(SYN_301_0174_genedata) <- paste0("SYN_301_0174_", colnames(SYN_301_0174_genedata))
SYN_301_0174_protdata <- SYN_301_0174_GEX$`Antibody Capture`
colnames(SYN_301_0174_protdata) <- paste0("SYN_301_0174_", colnames(SYN_301_0174_protdata))

SYN_301_0174 <- CreateSeuratObject(SYN_301_0174_genedata)
SYN_301_0174[["ADT"]] <- CreateAssayObject(SYN_301_0174_protdata)
SYN_301_0174$Patient <- "301_0174"
SYN_301_0174$Sample <- "SYN"

SYN_301_0174 <- combineExpression(AMP_2_RA_TCR$SYN_301_0174, SYN_301_0174, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_301_0174$cloneProp <- SYN_301_0174$Frequency/sum(table(SYN_301_0174$CTaa))

#PBL_301_0270
PBL_301_0270_GEX <- Read10X_h5("data/Sample_301_0270_PBL_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
PBL_301_0270_genedata <- PBL_301_0270_GEX$`Gene Expression`
colnames(PBL_301_0270_genedata) <- paste0("PBL_301_0270_", colnames(PBL_301_0270_genedata))
PBL_301_0270_protdata <- PBL_301_0270_GEX$`Antibody Capture`
colnames(PBL_301_0270_protdata) <- paste0("PBL_301_0270_", colnames(PBL_301_0270_protdata))

PBL_301_0270 <- CreateSeuratObject(PBL_301_0270_genedata)
PBL_301_0270[["ADT"]] <- CreateAssayObject(PBL_301_0270_protdata)
PBL_301_0270$Patient <- "301_0270"
PBL_301_0270$Sample <- "PBL"

PBL_301_0270 <- combineExpression(AMP_2_RA_TCR$PBL_301_0270, PBL_301_0270, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
PBL_301_0270$cloneProp <- PBL_301_0270$Frequency/sum(table(PBL_301_0270$CTaa))

#SYN_301_0270
SYN_301_0270_GEX <- Read10X_h5("data/Sample_301_0270_SYN_BT/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
SYN_301_0270_genedata <- SYN_301_0270_GEX$`Gene Expression`
colnames(SYN_301_0270_genedata) <- paste0("SYN_301_0270_", colnames(SYN_301_0270_genedata))
SYN_301_0270_protdata <- SYN_301_0270_GEX$`Antibody Capture`
colnames(SYN_301_0270_protdata) <- paste0("SYN_301_0270_", colnames(SYN_301_0270_protdata))

SYN_301_0270 <- CreateSeuratObject(SYN_301_0270_genedata)
SYN_301_0270[["ADT"]] <- CreateAssayObject(SYN_301_0270_protdata)
SYN_301_0270$Patient <- "301_0270"
SYN_301_0270$Sample <- "SYN"

SYN_301_0270 <- combineExpression(AMP_2_RA_TCR$SYN_301_0270, SYN_301_0270, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
SYN_301_0270$cloneProp <- SYN_301_0270$Frequency/sum(table(SYN_301_0270$CTaa))

AMP_2_RA_Seurat <- merge(PBL_300_0150, y = c(SYN_300_0150, PBL_300_0171, PBL_300_0171_2, SYN_300_0171, 
                                             PBL_300_0173, SYN_300_0173, PBL_300_0174, SYN_300_0174, PBL_300_0392, SYN_300_0392,
                                             PBL_300_0410, SYN_300_0410, PBL_300_0414, SYN_300_0414, SYN_300_0416, SYN_300_1883,
                                             PBL_300_1930, SYN_300_1930, PBL_301_0174, SYN_301_0174, PBL_301_0270, SYN_301_0270),project = "AMP2_RNA_TCR")

AMP_2_RA_Seurat$SamplePatient <- paste0(AMP_2_RA_Seurat$Patient, "_", AMP_2_RA_Seurat$Sample)

Idents(AMP_2_RA_Seurat) <- AMP_2_RA_Seurat$SamplePatient
AMP_2_RA_Seurat$SamplePatient <- as.character(AMP_2_RA_Seurat$SamplePatient)
AMP_2_RA_Seurat$SamplePatient <- factor(AMP_2_RA_Seurat$SamplePatient, levels=unique(AMP_2_RA_Seurat$SamplePatient))

AMP_2_RA_patient_metadata <- as.data.frame(read.csv("metadata/AMP_2_RA_repertoire_patient_metadata.csv"))
rownames(AMP_2_RA_patient_metadata) <- AMP_2_RA_patient_metadata[,1]
AMP_2_RA_patient_metadata[,1] <- NULL
AMP_2_RA_Seurat <- AddMetaData(AMP_2_RA_Seurat, AMP_2_RA_patient_metadata)

AMP_2_RA_Seurat$cloneType <- as.character(AMP_2_RA_Seurat$cloneType)
AMP_2_RA_Seurat$cloneType <- replace_na(AMP_2_RA_Seurat$cloneType, "No TCR")
Idents(AMP_2_RA_Seurat) <- AMP_2_RA_Seurat$cloneType

slot(AMP_2_RA_Seurat, "meta.data")$cloneType <- factor(slot(AMP_2_RA_Seurat, "meta.data")$cloneType, 
                                                       levels = c("Hyperexpanded (100 < X <= 500)", 
                                                                  "Large (20 < X <= 100)", 
                                                                  "Medium (5 < X <= 20)", 
                                                                  "Small (1 < X <= 5)", 
                                                                  "Single (0 < X <= 1)",
                                                                  "No TCR"))

current.clone.names <- c("Hyperexpanded (100 < X <= 500)", "Large (20 < X <= 100)", "Medium (5 < X <= 20)", "Small (1 < X <= 5)", "Single (0 < X <= 1)", "No TCR")
new.clone.names <- c("Hyperexpanded (100 < X)", "Large (20 < X <= 100)", "Medium (5 < X <= 20)", "Small (1 < X <= 5)", "Single (0 < X <= 1)", "No TCR information")
AMP_2_RA_Seurat@active.ident <- plyr::mapvalues(x = AMP_2_RA_Seurat@active.ident, from = current.clone.names, to = new.clone.names)
AMP_2_RA_Seurat$cloneType <- AMP_2_RA_Seurat@active.ident

slot(AMP_2_RA_Seurat, "meta.data")$cloneType <- factor(slot(AMP_2_RA_Seurat, "meta.data")$cloneType, 
                                                       levels = c("Hyperexpanded (100 < X)", 
                                                                  "Large (20 < X <= 100)", 
                                                                  "Medium (5 < X <= 20)", 
                                                                  "Small (1 < X <= 5)", 
                                                                  "Single (0 < X <= 1)", 
                                                                  "No TCR information"))

AMP_2_RA_Seurat$clonal <-ifelse(AMP_2_RA_Seurat$Frequency > 1, "Yes", "No")

AMP_2_RA_Seurat$clonal <-ifelse(AMP_2_RA_Seurat$cloneType ==  "Single (0 < X <= 1)" | AMP_2_RA_Seurat$cloneType == "No TCR information", "No", "Yes")

rm(PBL_300_0150, SYN_300_0150, PBL_300_0171, PBL_300_0171_2, SYN_300_0171, 
   PBL_300_0173, SYN_300_0173, PBL_300_0174, SYN_300_0174, PBL_300_0392, SYN_300_0392,
   PBL_300_0410, SYN_300_0410, PBL_300_0414, SYN_300_0414, SYN_300_0416, SYN_300_1883,
   PBL_300_1930, SYN_300_1930, PBL_301_0174, SYN_301_0174, PBL_301_0270, SYN_301_0270,
   PBL_300_0150_GEX, SYN_300_0150_GEX, PBL_300_0171_GEX, PBL_300_0171_2_GEX, SYN_300_0171_GEX, 
   PBL_300_0173_GEX, SYN_300_0173_GEX, PBL_300_0174_GEX, SYN_300_0174_GEX, PBL_300_0392_GEX, SYN_300_0392_GEX,
   PBL_300_0410_GEX, SYN_300_0410_GEX, PBL_300_0414_GEX, SYN_300_0414_GEX, SYN_300_0416_GEX, SYN_300_1883_GEX,
   PBL_300_1930_GEX, SYN_300_1930_GEX, PBL_301_0174_GEX, SYN_301_0174_GEX, PBL_301_0270_GEX, SYN_301_0270_GEX,
   PBL_300_0150_genedata, SYN_300_0150_genedata, PBL_300_0171_genedata, PBL_300_0171_2_genedata, SYN_300_0171_genedata, 
   PBL_300_0173_genedata, SYN_300_0173_genedata, PBL_300_0174_genedata, SYN_300_0174_genedata, PBL_300_0392_genedata, SYN_300_0392_genedata,
   PBL_300_0410_genedata, SYN_300_0410_genedata, PBL_300_0414_genedata, SYN_300_0414_genedata, SYN_300_0416_genedata, SYN_300_1883_genedata,
   PBL_300_1930_genedata, SYN_300_1930_genedata, PBL_301_0174_genedata, SYN_301_0174_genedata, PBL_301_0270_genedata, SYN_301_0270_genedata,
   PBL_300_0150_protdata, SYN_300_0150_protdata, PBL_300_0171_protdata, PBL_300_0171_2_protdata, SYN_300_0171_protdata, 
   PBL_300_0173_protdata, SYN_300_0173_protdata, PBL_300_0174_protdata, SYN_300_0174_protdata, PBL_300_0392_protdata, SYN_300_0392_protdata,
   PBL_300_0410_protdata, SYN_300_0410_protdata, PBL_300_0414_protdata, SYN_300_0414_protdata, SYN_300_0416_protdata, SYN_300_1883_protdata,
   PBL_300_1930_protdata, SYN_300_1930_protdata, PBL_301_0174_protdata, SYN_301_0174_protdata, PBL_301_0270_protdata, SYN_301_0270_protdata)

Idents(AMP_2_RA_Seurat) <- AMP_2_RA_Seurat$SamplePatient

AMP_2_RA_Seurat[["percent.mito"]] <- PercentageFeatureSet(AMP_2_RA_Seurat, pattern = "^MT-")

VlnPlot(AMP_2_RA_Seurat, features = c("nFeature_RNA"), group.by = "SamplePatient") + 
  theme(legend.position = "none") + 
  ggsave("analysis/figs/features_bySample.eps", width = 10, height = 6)

VlnPlot(AMP_2_RA_Seurat, features = c("percent.mito"), group.by = "SamplePatient") + 
  theme(legend.position = "none") + 
  ggsave("analysis/figs/mito_bySample.eps", width = 10, height = 6)

VlnPlot(AMP_2_RA_Seurat, features = c("percent.mito"), group.by = "Sample", pt.size = 0) + 
  theme(legend.position = "none") + 
  ggsave("analysis/figs/mito_bySampleType.eps", width = 4, height = 4)

VlnPlot(AMP_2_RA_Seurat, features = c("nCount_RNA"), group.by = "SamplePatient") + 
  theme(legend.position = "none") + 
  ggsave("analysis/figs/counts_bySample.eps", width = 10, height = 6)

cells_sample <- table(AMP_2_RA_Seurat$SamplePatient)

# Initial QC and Broad Clustering -----------------------------------------------------

AMP_2_RA_Seurat_sub <- subset(AMP_2_RA_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mito < 15)

rm(AMP_2_RA_Seurat)

cells_sample <- cbind(cells_sample, table(AMP_2_RA_Seurat_sub$SamplePatient))
colnames(cells_sample) <- c("pre-QC", "post-QC")
write.csv(cells_sample, file = "QC_nums.csv")

AMP_2_RA_Seurat_sub <- NormalizeData(AMP_2_RA_Seurat_sub, normalization.method = "LogNormalize", scale.factor = 10000)
AMP_2_RA_Seurat_sub <- FindVariableFeatures(AMP_2_RA_Seurat_sub, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AMP_2_RA_Seurat_sub), 10)

all.genes <- rownames(AMP_2_RA_Seurat_sub)
AMP_2_RA_Seurat_sub <- ScaleData(AMP_2_RA_Seurat_sub)
AMP_2_RA_Seurat_sub <- RunPCA(AMP_2_RA_Seurat_sub, features = VariableFeatures(object = AMP_2_RA_Seurat_sub))

dims_use <- 1:30
AMP_2_RA_Seurat_sub <- RunHarmony(AMP_2_RA_Seurat_sub, group.by.vars=c("SamplePatient"))
AMP_2_RA_Seurat_sub <- RunUMAP(object=AMP_2_RA_Seurat_sub, reduction="harmony", dims=dims_use, verbose=FALSE)
AMP_2_RA_Seurat_sub <- FindNeighbors(object=AMP_2_RA_Seurat_sub, reduction="harmony", dims=dims_use, verbose=FALSE)
AMP_2_RA_Seurat_sub <- FindClusters(object=AMP_2_RA_Seurat_sub, resolution=.2, verbose=FALSE)

DefaultAssay(AMP_2_RA_Seurat_sub) <- "ADT"
AMP_2_RA_Seurat_sub <- NormalizeData(AMP_2_RA_Seurat_sub, normalization.method = "CLR", margin = 2)
DefaultAssay(AMP_2_RA_Seurat_sub) <- "RNA"

DimPlot(AMP_2_RA_Seurat_sub) + 
  scale_color_brewer(palette = "Paired") 
ggsave("analysis/figs/UMAP_QCsubset_all_byCluster.eps", width = 5, height = 4)

FeaturePlot(AMP_2_RA_Seurat_sub, features = c("adt_CD4-TotalA", "adt_CD8a-TotalA", "CD4", "CD8A")) 
ggsave("analysis/figs/CD8_CD4_proteinVsRNA.png", width = 6, height = 5)

FeaturePlot(AMP_2_RA_Seurat_sub, features = c("adt_CD4-TotalA", "adt_CD8a-TotalA", "adt_CD27-TotalA", "adt_CD20-TotalA", "adt_IgD-TotalA", "adt_CD11c-TotalA"), min.cutoff = "q05", ncol = 3) 
ggsave("analysis/figs/FeaturePlot_CITEseq_panel.png", width = 9, height = 5)

DotPlot(AMP_2_RA_Seurat_sub, features = c("CD3E", "CD8A", "CD4", "IL7R", "NKG7", "FCER2", "MS4A1", "IGHD", 
                                          "CD79A", "IGHG1", "TNFRSF4", "PDCD1", "CTLA4", "FOXP3", "S100A8", 
                                          "LYZ", "CD68", "STMN1", "MKI67", "SPARC")) + RotatedAxis() + coord_flip() 
ggsave("analysis/figs/dotplot_broadCluster.eps", width = 5.7, height = 5)

# Subclustering of all T cells -----------------------------------------------------
#subsetting out T cell populations
AMP_2_RA_Seurat_sub_tcells <- subset(AMP_2_RA_Seurat_sub, idents = c(0,1,4,5,7,10))

AMP_2_RA_Seurat_sub_tcells <- NormalizeData(AMP_2_RA_Seurat_sub_tcells, normalization.method = "LogNormalize", scale.factor = 10000)
AMP_2_RA_Seurat_sub_tcells <- FindVariableFeatures(AMP_2_RA_Seurat_sub_tcells, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(AMP_2_RA_Seurat_sub_tcells), 10)

all.genes <- rownames(AMP_2_RA_Seurat_sub_tcells)
AMP_2_RA_Seurat_sub_tcells <- ScaleData(AMP_2_RA_Seurat_sub_tcells)
AMP_2_RA_Seurat_sub_tcells <- RunPCA(AMP_2_RA_Seurat_sub_tcells, features = VariableFeatures(object = AMP_2_RA_Seurat_sub_tcells))

dims_use <- 1:30
AMP_2_RA_Seurat_sub_tcells <- RunHarmony(AMP_2_RA_Seurat_sub_tcells, group.by.vars=c("SamplePatient"))
AMP_2_RA_Seurat_sub_tcells <- RunUMAP(object=AMP_2_RA_Seurat_sub_tcells, reduction="harmony", dims=dims_use, verbose=FALSE)
AMP_2_RA_Seurat_sub_tcells <- FindNeighbors(object=AMP_2_RA_Seurat_sub_tcells, reduction="harmony", dims=dims_use, verbose=FALSE)
AMP_2_RA_Seurat_sub_tcells <- FindClusters(object=AMP_2_RA_Seurat_sub_tcells, resolution=.5, verbose=FALSE)

nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
DimPlot(AMP_2_RA_Seurat_sub_tcells) + 
  scale_color_manual(values = mycolors) 
  ggsave("analysis/figs/UMAP_Tcellsubset_byCluster.eps", width = 5, height = 4)

cluster_markers_tcells <- FindAllMarkers(object=AMP_2_RA_Seurat_sub_tcells, only.pos = FALSE, logfc.threshold = 0.25, min.pct = 0.20)
write.csv(cluster_markers_tcells, "AMP_2_RA_tcell_cluster_markers.csv")

DotPlot(AMP_2_RA_Seurat_sub_tcells, features = c("CD3E", "CD8A", "CD4", "IL7R", "CD69", "TCF7", "SELL", "LEF1", 
                                                 "CCL4", "CCL5", "NKG7", "GZMK", "GZMA", "GZMB", "IFNG", "GNLY", "HLA-DRA", "EOMES", 
                                                 "CXCL13", "LAG3", "TNFRSF4", "PDCD1", "CTLA4", "TOX", "TIGIT", "MAF",
                                                 "FOXP3", "IL2RA", "ZNF683", "KLRD1", "KLRG1", "ZEB2", "KLRK1", "TRDC", "TRGC1",
                                                 "LGALS1", "CXCR3", "MKI67", "HAVCR2", "CD38", "CD79A", "MS4A1")) + RotatedAxis() 
ggsave("analysis/figs/dotplot_TCellCluster.eps", width = 11, height = 4)

#getting final idea of CD8 clusters before subsetting
DotPlot(AMP_2_RA_Seurat_sub_tcells, features = c("CD3E", "CD8A", "CD4", "GZMK", "GZMB", "GZMA", "IFNG", "GNLY")) + RotatedAxis() 
ggsave("analysis/figs/cd8/dotplot_GZMB_GZMK_TCellCluster.eps", width = 5, height = 4)

DotPlot(AMP_2_RA_Seurat_sub_tcells, features = c("CD3E", "CD8A", "CD4", "GZMK", "GZMB", "GZMA", "IFNG", "GNLY"), split.by = "Sample", cols =c("#00AFBB", "#D55E00")) + RotatedAxis()
ggsave("analysis/figs/cd8/dotplot_GZMB_GZMK_TCellCluster_splitSample.eps", width = 5, height = 8)

nb.cols <- 15
mycolors_cd8 <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
mycolors_cd8[c(1:2,4:6,8,11:12,14:15)] <- "lightgrey"
DimPlot(AMP_2_RA_Seurat_sub_tcells) + 
  scale_color_manual(values = mycolors_cd8) 
ggsave("analysis/figs/cd8/UMAP_Tcellsubset_CD8highlight.eps", width = 5, height = 4)


# Subclustering of CD8 cells -----------------------------------------------------

AMP_2_RA_Seurat_sub_CD8 <- subset(AMP_2_RA_Seurat_sub_tcells, idents = c(2,6,8,9,12))
AMP_2_RA_Seurat_sub_CD8 <- subset(AMP_2_RA_Seurat_sub_CD8, subset = `adt_CD8a-TotalA` > 1)

#reclustering
AMP_2_RA_Seurat_sub_CD8 <- FindVariableFeatures(AMP_2_RA_Seurat_sub_CD8, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(AMP_2_RA_Seurat_sub_CD8), 10)

all.genes <- rownames(AMP_2_RA_Seurat_sub_CD8)
AMP_2_RA_Seurat_sub_CD8 <- ScaleData(AMP_2_RA_Seurat_sub_CD8)
AMP_2_RA_Seurat_sub_CD8 <- RunPCA(AMP_2_RA_Seurat_sub_CD8, features = VariableFeatures(object = AMP_2_RA_Seurat_sub_CD8))

dims_use = 1:20
AMP_2_RA_Seurat_sub_CD8 <- RunHarmony(AMP_2_RA_Seurat_sub_CD8, group.by.vars=c("SamplePatient"))
AMP_2_RA_Seurat_sub_CD8 <- RunUMAP(object=AMP_2_RA_Seurat_sub_CD8, reduction="harmony", dims=dims_use, verbose=FALSE)
AMP_2_RA_Seurat_sub_CD8 <- FindNeighbors(object=AMP_2_RA_Seurat_sub_CD8, reduction="harmony", dims=dims_use, verbose=FALSE)
AMP_2_RA_Seurat_sub_CD8 <- FindClusters(object=AMP_2_RA_Seurat_sub_CD8, resolution=.2, verbose=FALSE)

mycolorsCD8recluster <- c("#74C476", "#999999", "#084594", "#ca636d", "#eebf7a", "#7a3669", "#d0c4a6")

#Figure 5A
DimPlot(AMP_2_RA_Seurat_sub_CD8) +
  scale_color_manual(values = mycolorsCD8recluster, labels = c("0 - GZMK+", "1 - Naive", "2 - GZMB+ PRF1+", "3 - Early Memory CD8", "4 - gdT cells", "5 - NK cells/gdT cells", "6 - Proliferating cells")) +
  theme( legend.key.size = unit(1, 'cm'), #change legend key size
         legend.key.height = unit(.4, 'cm'), #change legend key height
         legend.key.width = unit(.4, 'cm'), #change legend key width
         legend.text = element_text(size=9),
         axis.line.x.bottom = element_blank(),
         axis.line.y.left   = element_blank(),
         axis.line.y.right  = element_blank(),
         axis.line.x.top  = element_blank(),
         panel.border = element_rect(color = "black", size = 0.5)) +
  xlab("UMAP1") + ylab("UMAP2")
ggsave("analysis/figs/cd8/UMAP_CD8subset_recluster_byCluster.eps", width = 6, height = 4)

patientcols <- c(colorRampPalette(brewer.pal(8, "Set2"))(13)[1:5], colorRampPalette(brewer.pal(8, "Set1"))(13)[1:6], colorRampPalette(brewer.pal(8, "Set3"))(13)[1:6])

#Figure S11A
DimPlot(AMP_2_RA_Seurat_sub_CD8, group.by = "Patient", shuffle = T) +
  scale_color_manual(values = patientcols) +
  ggtitle(NULL)
ggsave("analysis/figs/cd8/UMAP_CD8subset_recluster_byPatient.eps", width = 6, height = 4)

#Figure S11B
DimPlot(AMP_2_RA_Seurat_sub_CD8, group.by = "Sample", shuffle = T) +
  scale_fill_manual(values = c("#00AFBB", "#D55E00")) +
  ggtitle(NULL)
ggsave("analysis/figs/cd8/UMAP_CD8subset_recluster_bySample.eps", width = 6, height = 4)

DimPlot(AMP_2_RA_Seurat_sub_CD8, split.by = "SamplePatient", shuffle = T, ncol = 6) +
  scale_color_manual(values = mycolorsCD8recluster, labels = c("0 - GZMK+", "1 - Naive", "2 - GZMB+ PRF1+", "3 - Early Memory CD8", "4 - gdT cells", "5 - NK cells/gdT cells", "6 - Proliferating cells")) +
  ggtitle(NULL)
ggsave("analysis/figs/cd8/UMAP_CD8subset_recluster_byCluster_splitSamplePatient.eps", width = 14, height = 8)

#DEGs for clusters
cluster_markers_CD8 <- FindAllMarkers(object=AMP_2_RA_Seurat_sub_CD8, only.pos=FALSE, logfc.threshold=0.25, min.pct = 0.50)
write.csv(cluster_markers_CD8, "AMP_2_RA_CD8_recluster_markers.csv")

AMP_2_RA_Seurat_sub_CD8$sampleTypeCluster <- paste0(AMP_2_RA_Seurat_sub_CD8$seurat_clusters, "-", AMP_2_RA_Seurat_sub_CD8$Sample)
write.csv(table(AMP_2_RA_Seurat_sub_CD8$Patient, AMP_2_RA_Seurat_sub_CD8$sampleTypeCluster), "AMP_2_CD8_cellcount_patient_clusterSampleSplit.csv")

DotPlot(AMP_2_RA_Seurat_sub_CD8, features = c("CD3E", "CD8A", "CD4", "SELL", "TCF7", "LEF1", "CD69", "IL7R", "ENTPD1", "GZMA", "GZMB", "GZMH", "GZMK", "IFNG", "GNLY", "PRF1", "NKG7", "TYROBP", "ZNF683", "ICOS", "TIGIT", "PDCD1", "CTLA4", "TNFRSF9", "TRGC1", "TRDC","MKI67")) + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("analysis/figs/cd8/dotplot_CD8subset_recluster.eps", width = 8, height = 4)

#Figure 5B
DotPlot(AMP_2_RA_Seurat_sub_CD8, features = rev(c("CD3E", "CD8A", "CD4", "SELL", "TCF7", "LEF1", "CD69", "IL7R", "ENTPD1", "GZMA", "GZMB", "GZMH", "GZMK", "IFNG", "GNLY", "PRF1", "NKG7", "TYROBP", "ZNF683", "ICOS", "TIGIT", "PDCD1", "CTLA4", "TNFRSF9", "TRGC1", "TRDC","MKI67"))) + 
  coord_flip() + 
  scale_x_discrete(position = "top") + 
  scale_y_discrete(position = "right") + 
  theme(legend.position="left") 
ggsave("analysis/figs/cd8/dotplot_CD8subset_recluster_rotated.eps", width = 6, height = 7)

#Figure S11C
AMP_2_RA_Seurat_sub_CD8@meta.data %>%
  group_by(SamplePatient, seurat_clusters) %>%
  dplyr::count() %>%
  group_by(SamplePatient) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=SamplePatient,y=Percent, fill=seurat_clusters)) +
  geom_col() +
  scale_fill_manual(values = mycolorsCD8recluster) +
  ggtitle("Cluster Composition per Sample") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
ggsave("analysis/figs/cd8/CD8broadClusterDist_bySample.eps", width = 8, height = 4)

#Figure S11D
AMP_2_RA_Seurat_sub_CD8@meta.data %>%
  group_by(seurat_clusters, Sample) %>%
  dplyr::count() %>%
  group_by(seurat_clusters) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=Percent, fill=Sample)) +
  scale_fill_manual(values = rev(c("#00AFBB", "#D55E00"))) +
  geom_col() +
  ggtitle("Cluster Composition per Cell Source") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  
ggsave("analysis/figs/cd8/CD8ClusterDist_bySampleType.eps", width = 5, height = 4)

AMP_2_RA_Seurat_sub_CD8$CTaa_sampleType <- paste0(AMP_2_RA_Seurat_sub_CD8$CTaa, "-", AMP_2_RA_Seurat_sub_CD8$Sample)

Idents(AMP_2_RA_Seurat_sub_CD8) <- AMP_2_RA_Seurat_sub_CD8$CTaa_sampleType


#Figure 5C

#300-0150
PBL1_PBLplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CASYSGSARQLTF_CASSDYRDRGANQPQHF-PBL", "NA_CASSDYRDRGANQPQHF-PBL"))
P1_0150_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL1_PBLplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL1_SYNplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CASYSGSARQLTF_CASSDYRDRGANQPQHF-SYN", "NA_CASSDYRDRGANQPQHF-SYN"))
P1_0150_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL1_SYNplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL2_PBLplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CILGLMDSNYQLIW_CATSRGSIGVSNTGELFF-PBL"))
P2_0150_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL2_PBLplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL2_SYNplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CILGLMDSNYQLIW_CATSRGSIGVSNTGELFF-SYN"))
P2_0150_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL2_SYNplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL3_PBLplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAAYQTGANNLFF_CAIRTGQDNEQFF-PBL"))
P3_0150_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL3_PBLplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL3_SYNplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAAYQTGANNLFF_CAIRTGQDNEQFF-SYN"))
P3_0150_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL3_SYNplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)


SYN1_PBLplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CLVGSSGGSYIPTF_CASSPRGGGTGDTGELFF-PBL"))
S1_0150_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN1_SYNplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CLVGSSGGSYIPTF_CASSPRGGGTGDTGELFF-SYN"))
S1_0150_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN1_SYNplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN2_PBLplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVEATTSGTYKYIF_CAISDRGLQGSEKLFF-PBL"))
S2_0150_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN2_PBLplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN2_SYNplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVEATTSGTYKYIF_CAISDRGLQGSEKLFF-SYN"))
S2_0150_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN2_SYNplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN3_PBLplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVNRGRDAGGTSYGKLTF_CASSTTSGSYNEQFF-PBL"))
S3_0150_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN3_PBLplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN3_SYNplot_0150 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVNRGRDAGGTSYGKLTF_CASSTTSGSYNEQFF-SYN"))
S3_0150_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN3_SYNplot_0150), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)


#300-1930
PBL1_PBLplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAENGPGGGKIIF_CASSYSTDEQYF-PBL"))
P1_1930_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL1_PBLplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL1_SYNplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAENGPGGGKIIF_CASSYSTDEQYF-SYN"))
P1_1930_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL1_SYNplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL2_PBLplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CARGSQGNLIF_CASSSWASPNTGELFF-PBL"))
P2_1930_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL2_PBLplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL2_SYNplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CARGSQGNLIF_CASSSWASPNTGELFF-SYN"))
P2_1930_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL2_SYNplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL3_PBLplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAYRRGGSEKLVF_CASSQLPRTEAFF-PBL"))
P3_1930_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL3_PBLplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL3_SYNplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAYRRGGSEKLVF_CASSQLPRTEAFF-SYN"))
P3_1930_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL3_SYNplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)


SYN1_PBLplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVRDNTDKLIF_CATSDSRNPRGGPSTYNEQFF-PBL"))
S1_1930_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN1_PBLplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN1_SYNplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVRDNTDKLIF_CATSDSRNPRGGPSTYNEQFF-SYN"))
S1_1930_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN1_SYNplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN2_PBLplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVSARYNTDKLIF_CASSSMGPPSYEQYF-PBL"))
S2_1930_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN2_PBLplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN2_SYNplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVSARYNTDKLIF_CASSSMGPPSYEQYF-SYN"))
S2_1930_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN2_SYNplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN3_PBLplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CARGSQGNLIF_CASSSWASPNTGELFF-PBL"))
S3_1930_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN3_PBLplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN3_SYNplot_1930 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CARGSQGNLIF_CASSSWASPNTGELFF-SYN"))
S3_1930_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN3_SYNplot_1930), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

a <- (P1_0150_P | P1_0150_S) / (P2_0150_P | P2_0150_S) / (P3_0150_P | P3_0150_S) 
b <- (S1_0150_P | S1_0150_S) / (S2_0150_P | S2_0150_S) / (S3_0150_P | S3_0150_S) 
c <- (P1_1930_P | P1_1930_S) / (P2_1930_P | P2_1930_S) / (P3_1930_P | P3_1930_S) 
d <- (S1_1930_P | S1_1930_S) / (S2_1930_P | S2_1930_S) / (S3_1930_P | S3_1930_S) 

a | b | c | d


#301-0174
PBL1_PBLplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVRVYSTLTF_CASRTGTAGNTIYF-PBL"))
P1_0174_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL1_PBLplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL1_SYNplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVRVYSTLTF_CASRTGTAGNTIYF-SYN"))
P1_0174_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL1_SYNplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL2_PBLplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVITFSGGYNKLIF;CAVRPSNAGKSTF_CASSRLAGGPSEQFF-PBL"))
P2_0174_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL2_PBLplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL2_SYNplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAVITFSGGYNKLIF;CAVRPSNAGKSTF_CASSRLAGGPSEQFF-SYN")) #not found in this
P2_0174_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL3_PBLplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAGPKIQGAQKLVF_CSARDGTGGREQYF-PBL"))
P3_0174_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL3_PBLplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

PBL3_SYNplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAGPKIQGAQKLVF_CSARDGTGGREQYF-SYN"))
P3_0174_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(PBL3_SYNplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)


SYN1_PBLplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAAQETSYDKVIF_CASSLTASGNQPQHF-PBL"))
S1_0174_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN1_PBLplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN1_SYNplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAAQETSYDKVIF_CASSLTASGNQPQHF-SYN"))
S1_0174_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN1_SYNplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN2_PBLplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAASETSYDKVIF_CASSPTQAGNQPQHF-PBL"))
S2_0174_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN2_PBLplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN2_SYNplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CAASETSYDKVIF_CASSPTQAGNQPQHF-SYN"))
S2_0174_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN2_SYNplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN3_PBLplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CVVKGREYGNKLVF;CAAINAGKSTF_CASSLGRDPLNEQFF-PBL"))
S3_0174_P <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN3_PBLplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#101010") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)

SYN3_SYNplot_0174 <- WhichCells(AMP_2_RA_Seurat_sub_CD8, idents = c("CVVKGREYGNKLVF;CAAINAGKSTF_CASSLGRDPLNEQFF-SYN"))
S3_0174_S <- DimPlot(AMP_2_RA_Seurat_sub_CD8, pt.size = .5, cells.highlight= list(SYN3_SYNplot_0174), cols= "#dedede", sizes.highlight = 1.5, cols.highlight = "#8b13ef") +
  theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), panel.border = element_rect(color = "black", size = 1)) + ylab(NULL) + xlab(NULL)


a <- (P1_0150_P | P1_0150_S) / (P2_0150_P | P2_0150_S) / (P3_0150_P | P3_0150_S) 
b <- (S1_0150_P | S1_0150_S) / (S2_0150_P | S2_0150_S) / (S3_0150_P | S3_0150_S) 
e <- (P1_0174_P | P1_0174_S) / (P2_0174_P | P2_0174_S) / (P3_0174_P | P3_0174_S) 
f <- (S1_0174_P | S1_0174_S) / (S2_0174_P | S2_0174_S) / (S3_0174_P | S3_0174_S) 

a | b | e | f

#clonal overlap
test <- subset(AMP_2_RA_Seurat_sub_tcells, idents = c(2,6,12))

current.cluster.ids <- c(2,6,12)
new.cluster.ids <- c("GZMK+", "GZMB+", "Cycling")
test@active.ident <- plyr::mapvalues(x = test@active.ident, from = current.cluster.ids, to = new.cluster.ids)
test$clusters_names <- test@active.ident

#Figure S11E
test$clusters_names_tissue <- paste0(test$clusters_names, "_", test$Sample)
GZM_clonotypes_bySample <- expression2List(test, group = "clusters_names_tissue")
clonalOverlap(GZM_clonotypes_bySample, cloneCall = "aa") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_x_discrete(limits = c("GZMK+_SYN","GZMK+_PBL", "GZMB+_SYN", "GZMB+_PBL", "Cycling_SYN", "Cycling_PBL")) +
  labs(fill = "Overlap") +
  theme(legend.position = c(0.10, 0.35)) +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 5))
ggsave("analysis/figs/cd8/CD8_GZM_overlap_Tissue.eps", width = 5, height = 3)

#Figure S11G
AMP_2_RA_Seurat_sub_CD8@meta.data %>%
  group_by(sampleTypeCluster, cloneType) %>%
  count() %>%
  group_by(sampleTypeCluster) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=sampleTypeCluster,y=Percent, fill=cloneType)) +
  geom_col() +
  scale_fill_brewer(palette = "Reds", direction = -1, name = "Clone Type") +
  ggtitle(NULL) +
  xlab("Cluster") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(.4, 'cm'), #change legend key height
        legend.key.width = unit(.4, 'cm'), #change legend key width
        legend.text = element_text(size=9)) 
ggsave("analysis/figs/cd8/CD8_subcluster_Clonality_byCluster_sampleType.eps", width = 5, height = 2.5)


CD8_clonotypes <- CD8_clonotypes[c("GZMK+", "GZMB+", "GZM Mix", "Cycling")]

write.csv(table(AMP_2_RA_Seurat_sub_CD8$CTaa, AMP_2_RA_Seurat_sub_CD8$seurat_clusters), "CD8clusters_clonotypes.csv")

#preparing data for Figure 5D-H calculations
AMP_CD8_meta <- AMP_2_RA_Seurat_sub_CD8@meta.data

AMP_CD8_counts <- AMP_2_RA_Seurat_sub_CD8@assays$RNA@counts
AMP_CD8_counts <- as.data.frame(AMP_CD8_counts)
AMP_CD8_counts <- t(AMP_CD8_counts)
AMP_CD8_counts <- as_tibble(AMP_CD8_counts, rownames = NA)

AMP_CD8_counts_GZM <- cbind(AMP_CD8_counts$GZMB, AMP_CD8_counts$GZMK)
rownames(AMP_CD8_counts_GZM) <- rownames(AMP_CD8_counts)
colnames(AMP_CD8_counts_GZM) <- c("GZMB", "GZMK")

AMP_CD8_meta_GZM <- merge(AMP_CD8_meta, AMP_CD8_counts_GZM, by = 0)

top3CDR3 <- AMP_CD8_meta_GZM %>%
  filter(CTaa != "<NA>") %>%
  add_count(seurat_clusters, name = "cluster_totaln") %>%
  add_count(SamplePatient, name = "samplepatient_totaln") %>%
  add_count(Sample, name = "sampletype_totaln") %>%
  add_count(Patient, name = "patient_totaln") %>%
  dplyr::count(CTaa, SamplePatient, Patient, Sample, seurat_clusters, cluster_totaln, samplepatient_totaln, sampletype_totaln, patient_totaln, sort = TRUE) %>%
  group_by(SamplePatient) %>% 
  dplyr::slice(1:3) %>% 
  pull(CTaa)

top6CDR3 <- AMP_CD8_meta_GZM %>%
  filter(CTaa != "<NA>") %>%
  add_count(seurat_clusters, name = "cluster_totaln") %>%
  add_count(SamplePatient, name = "samplepatient_totaln") %>%
  add_count(Sample, name = "sampletype_totaln") %>%
  add_count(Patient, name = "patient_totaln") %>%
  dplyr::count(CTaa, SamplePatient, Patient, Sample, seurat_clusters, cluster_totaln, samplepatient_totaln, sampletype_totaln, patient_totaln, sort = TRUE) %>%
  group_by(SamplePatient) %>% 
  dplyr::slice(1:6) %>% 
  pull(CTaa)

#top6
AMP_CD8_meta_GZM_filtered_samplePatient <- AMP_CD8_meta_GZM %>%
  filter(CTaa != "<NA>") %>%
  mutate(cluster_samplepatient = paste0(seurat_clusters, "-", SamplePatient)) %>%
  add_count(seurat_clusters, name = "cluster_totaln") %>%
  add_count(SamplePatient, name = "samplepatient_totaln") %>%
  add_count(Sample, name = "sampletype_totaln") %>%
  add_count(Patient, name = "patient_totaln") %>%
  add_count(cluster_samplepatient, name = "cluster_samplepatient_totaln") %>%
  dplyr::count(CTaa, SamplePatient, Patient, Sample, seurat_clusters, cluster_totaln, samplepatient_totaln, sampletype_totaln, patient_totaln, cluster_samplepatient_totaln, sort = TRUE) %>%
  group_by(SamplePatient) %>%
  filter(CTaa %in% top6CDR3) %>%
  filter(seurat_clusters %in% c(0,2)) %>%
  mutate(cluster_freq = (n / cluster_totaln)) %>%
  mutate(samplepatient_freq = (n / samplepatient_totaln)) %>%
  mutate(sampletype_freq = (n / sampletype_totaln)) %>%
  mutate(patient_freq = (n / patient_totaln)) %>%
  mutate(clusterpatient_percent = (n / cluster_samplepatient_totaln)) %>%
  mutate(SampleCluster = paste(Sample, seurat_clusters, sep = "_")) %>%
  select(Patient, CTaa, clusterpatient_percent, SampleCluster) %>%
  spread(SampleCluster, clusterpatient_percent)

#top6
AMP_CD8_meta_GZM_filtered_patient <- AMP_CD8_meta_GZM %>%
  filter(CTaa != "<NA>") %>%
  mutate(cluster_patient = paste0(seurat_clusters, "-", Patient)) %>%
  add_count(seurat_clusters, name = "cluster_totaln") %>%
  add_count(SamplePatient, name = "samplepatient_totaln") %>%
  add_count(Sample, name = "sampletype_totaln") %>%
  add_count(Patient, name = "patient_totaln") %>%
  add_count(cluster_patient, name = "cluster_patient_totaln") %>%
  dplyr::count(CTaa, SamplePatient, Patient, Sample, seurat_clusters, cluster_totaln, samplepatient_totaln, sampletype_totaln, patient_totaln, cluster_patient_totaln, sort = TRUE) %>%
  group_by(SamplePatient) %>%
  filter(CTaa %in% top6CDR3) %>%
  filter(seurat_clusters %in% c(0,2)) %>%
  mutate(cluster_freq = (n / cluster_totaln)) %>%
  mutate(samplepatient_freq = (n / samplepatient_totaln)) %>%
  mutate(sampletype_freq = (n / sampletype_totaln)) %>%
  mutate(patient_freq = (n / patient_totaln)) %>%
  mutate(clusterpatient_percent = (n / cluster_patient_totaln)) %>%
  mutate(SampleCluster = paste(Sample, seurat_clusters, sep = "_")) %>%
  select(Patient, CTaa, clusterpatient_percent, SampleCluster) %>%
  spread(SampleCluster, clusterpatient_percent)

clusterFreqAMP <- read.csv("AMP_2_CD8_GZM_filtered_freqbyClusterbySample_20210823.csv")
clusterFreqAMP <- as.tibble(clusterFreqAMP)

cd8_clonotypes_byCluster <- table(AMP_2_RA_Seurat_sub_CD8$CTaa, AMP_2_RA_Seurat_sub_CD8$seurat_clusters)
cd8clones <- read.csv("AMP_2_CD8_clonotypeClusterCounts_all.csv")
cd8clones <- tibble(cd8clones)

cd8clones_selected <- cd8clones %>%
  filter(CTaa %in% clusterFreqAMP$CTaa)
write.csv(cd8clones_selected, file = "AMP_2_CD8_clonotypeClusterCounts_selected.csv")

cd8clones_all <- cd8clones

cd8clones <- read.csv("AMP_2_CD8_clonotypeClusterCounts_all.csv")
cd8clones_naive <- tibble(cd8clones)
cd8clones_naive <- cd8clones_naive %>%
  select(CTaa, C1) %>%
  add_tally(C1, name = "naive_totaln") %>%
  mutate(naivecluster_freq = (C1 / naive_totaln)) %>%
  select(CTaa, naivecluster_freq) %>%
  filter(CTaa %in% clusterFreqAMP$CTaa)
write.csv(cd8clones_naive, "AMP_2_CD8_naive_cloneFreq_selected.csv")

AMP_2_RA_Seurat_sub_CD8$seurat_clusters_sampleType <- paste0("C", AMP_2_RA_Seurat_sub_CD8$seurat_clusters, "-", AMP_2_RA_Seurat_sub_CD8$Sample)
cd8_clonotypes_byClusterSample <- table(AMP_2_RA_Seurat_sub_CD8$CTaa, AMP_2_RA_Seurat_sub_CD8$seurat_clusters_sampleType)
cd8clones_ClusterSample <- read.csv("AMP_2_TCR_clonotypeClusterSampleTypeCounts_all.csv")
cd8clones_ClusterSample <- tibble(cd8clones_ClusterSample)
cd8clones_ClusterSample_selected <- cd8clones_ClusterSample %>%
  filter(CTaa %in% clusterFreqAMP$CTaa)
write.csv(cd8clones_ClusterSample_selected, file = "AMP_2_TCR_clonotypeClusterSampleTypeCounts_selected.csv")

cd8_clonotypes_bySampleType <- table(AMP_2_RA_Seurat_sub_CD8$CTaa, AMP_2_RA_Seurat_sub_CD8$Sample)
cd8clones_bySampleType <- read.csv("AMP_2_TCR_clonotypeSampleTypeCounts_all.csv")
cd8clones_bySampleType <- tibble(cd8clones_bySampleType)
cd8clones_bySampleType_selected <- cd8clones_bySampleType %>%
  filter(CTaa %in% clusterFreqAMP$CTaa)
write.csv(cd8clones_bySampleType_selected, file = "AMP_2_TCR_clonotypeSampleTypeCounts_selected.csv")

