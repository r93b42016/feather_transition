
# load package
library("DESeq2")

setwd("C:\\Users\\r93b4\\Desktop\\Myprojects\\feather_transition_keratinocyte\\feather transition_YoMing\\To_YM\\DEG")
getwd()

input <- read.table(file.choose(), header=TRUE, row.names="gene_id", sep=",")
head(input)
summary(input)


#remove rows by external list
input <- input[ ! row.names(input) %in% row.names(input_filter), ]

#subset rows by external list
input<-subset(input, row.names(input) %in% row.names(input_filter))
head(input)
summary(input)


# remove read counts sum rows< 10
input = input[ rowSums(input)>10, ]

##contour D3 vs D5
input_D<-input[,c(1,2,3,7,8,9)]
head(input_D)
#contour_comp
condition <- factor(c("D3","D3","D3","D5","D5","D5"))

##contour D3 vs D7
input_D<-input[,c(1,2,3,13,14,15)]
head(input_D)
#contour_comp
condition <- factor(c("D3","D3","D3","D7","D7","D7"))

##contour D6_7 vs D3_4_5
input_D<-input[11:25]
head(input_D)

##down E14_9E vs E7_9D
input_N<-input[1:10]
head(input_N)

#contour_comp
condition <- factor(c("345","345","345","345","345","345","345","345","345","67","67","67","67","67","67"))

#down_comp
condition <- factor(c("79","79","79","79","79","149","149","149","149","79"))



# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(input_D), condition))
#coldata$treatment <- treatment
#coldata$time <- time

#dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeqDataSetFromMatrix(countData=input_D, colData=coldata, design=~ condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds, contrast=c("condition","D3","D7"))
#res <- results(dds, contrast=c("condition","D3","E16"),lfcThreshold= log2(1.2),alpha = 0.1)
summary(res)
#remove any rows with NA
res <- res[complete.cases(res),]  
summary(res)

## Order by adjusted p-value
res <- res[order(res$padj), ]


## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "gene_id"
head(resdata)

## Write results
write.csv(resdata, file="DEGs_D3vsD7.csv")

## Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))

##to output dds (optional)
countsdds <- counts(dds, normalized = TRUE)
write.csv(countsdds, file="FTdowny_DEGs_dds.csv")

#output rld
rld_assay <- assay(rlog(dds, blind=FALSE))
write.table(rld_assay,"FTdowny_DEGs_rld.csv")
