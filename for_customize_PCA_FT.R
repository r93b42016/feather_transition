##customized PCA
# Load dds data into R

setwd("C:\\Users\\r93b4\\Desktop\\Collaborations\\feather_transition_keratinocyte\\For_plot")
getwd()
library("DESeq2")

# We have used gene names as row names  and records are separated by space
data <- read.table(file.choose(), header=TRUE, row.names="gene_id",sep=",")
head(data)
summary(data)

#log transform (optional)
#data[data < 0.000000001] <- 0.0000000001
#data<-log2(data)

condition <- factor(c("D3","D3","D3","D4","D4","D4","D5","D5","D5","D6","D6","D6","D7","D7","D7","E12","E12","E14","E14","E16","E16"))
condition <- factor(c("E12","E12","E14","E14","E16","E16"))

(coldata <- data.frame(row.names=colnames(data), condition))

#dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeqDataSetFromMatrix(countData=data, colData=coldata, design=~ condition)

dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))

# remove read dds sum rows< ?
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
head(dds)

##to output dds (optional)
countsdds <- counts(dds, normalized = TRUE)
write.csv(countsdds, file="NGanno_e_allexpressed_dds.csv")

#output rld
rld_assay <- assay(rlog(dds, blind=FALSE))
write.table(rld_assay,"NGanno_e_gene_count_matrix_Monly_rld.csv")


# Create a prcomp object. For this we need to transpose the data frame so that sample
# names are in rows and gene names are in columns
input <- read.table(file.choose(), header=TRUE, row.names="gene_id",sep = ",")
head(input)

input <- as.matrix(input)

pca_data=prcomp(t(input))

# Let us calculated the variances covered by components.

pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)

# create a data frame with principal component 1 (PC1), PC2, Conditions and sample names
df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = c("D3_11","D3_12","D3_13","D4_11","D4_12","D4_14","D5_8","D5_9","D5_10","D6_11","D6_12","D6_13","D7_8","D7_9","D7_11"))
df_pca_data

#plot, color by color, shape by region
library(ggplot2)
library(grid)
library(gridExtra)
#df <- as.data.frame(data)
ggplot(df_pca_data, aes(PC1,PC2, color = sample, ntop=1000))+geom_point(size=3,alpha=1) + geom_text(angle = -20, size = 3, vjust = 0.2, hjust = -0.5,  check_overlap = F, aes(label = rownames(pca_data$x)))+
  scale_color_manual(values = c("D3_11"= "red","D3_12"= "red","D3_13" = "red", "D4_11" = "black","D4_12" = "black","D4_14" = "black", "D5_8" = "gray", "D5_9" = "gray","D5_10" = "gray","D6_11" = "blue","D6_12" = "blue","D6_13" = "blue","D7_8" = "orange","D7_9" = "orange","D7_11" = "orange"))+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))+
  xlim(-40, 75)+
  ylim(-25, 25)+
  theme_bw()
theme_classic()

