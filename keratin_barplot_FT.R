install.packages(c("ggplot2", "reshape2", "dplyr"))

library(ggplot2)

#input "all_bKRT_all_tissues_TPM_for_plot.csv"
setwd("C:\\Users\\r93b4\\Desktop\\Myprojects\\feather_transition_keratinocyte\\Figures\\Fig_EDC_Bkeratins\\Bar_plot")
setwd("C:\\Users\\r93b4\\Desktop\\Myprojects\\feather_transition_keratinocyte\\Figures\\Fig_Bkrts_exp\\Bar_plot")
getwd()
input_both <- read.table(file.choose(), header=TRUE,sep = ",")
head(input_both)
input_filtered <- filter(input_both, stage != "E7E" | stage != "E7D"|stage != "E9E"|stage != "E9D")
input_both_filtered<-subset(input_both, stage!="E7E" & stage!="E7D" &stage!="E9E" &stage!="E9D")
#level the gene names
input_both$stage <- factor(input_both$stage,levels=c("E12","E14","E16","D3","D4","D5","D6","D7"))
input_both$Chr <- factor(input_both$Chr,levels=c("Chr1","Chr2","Chr6","Chr7","Chr10","Chr25","Chr27"))
#level the gene names
input_both_filtered$stage <- factor(input_both_filtered$stage,levels=c("E12","E14","E16","D3","D4","D5","D6","D7"))
input_both_filtered$Chr <- factor(input_both_filtered$Chr,levels=c("Chr1","Chr2","Chr6","Chr7","Chr10","Chr25","Chr27"))
view(input_both_filtered)
View(input_both_filtered)
View(input_both)
input_both <- read.table(file.choose(), header=TRUE,sep = ",")
head(input_both)
input_both_filtered<-subset(input_both, stage!="E7E" & stage!="E7D" &stage!="E9E" &stage!="E9D")
View(input_both)
View(input_both_filtered)
#level the gene names
input_both_filtered$stage <- factor(input_both_filtered$stage,levels=c("E12","E14","E16","D3","D4","D5","D6","D7"))
input_both_filtered$Chr <- factor(input_both_filtered$Chr,levels=c("Chr1","Chr2","Chr6","Chr7","Chr10","Chr25","Chr27"))
input_both_filtered$chr <- factor(input_both_filtered$chr,levels=c("Chr1","Chr2","Chr6","Chr7","Chr10","Chr25","Chr27"))
head(input_both_filtered)
input_both_filtered<-subset(input_both, stage!="E7E" & stage!="E7D" &stage!="E9E" &stage!="E9D")
head(input_both_filtered)
#level the gene names
input_both_filtered$stage <- factor(input_both_filtered$stage,levels=c("E12","E14","E16","D3","D4","D5","D6","D7"))
input_both_filtered$chr <- factor(input_both_filtered$chr,levels=c("Chr1","Chr2","Chr6","Chr7","Chr10","Chr25","Chr27"))
head(input_both_filtered)
View(input_both_filtered)
input_both_filtered<-subset(input_both, stage!="E7E" & stage!="E7D" &stage!="E9E" &stage!="E9D")
View(input_both)
uniq(input_both_filtered$chr)
unique(input_both_filtered$chr)
View(input_both_filtered)
#level the gene names
input_both_filtered$stage <- factor(input_both_filtered$stage,levels=c("E12F","E14F","E16F","D3","D4","D5","D6","D7"))
input_both_filtered$chr <- factor(input_both_filtered$chr,levels=c("Chr1","Chr2","Chr6","Chr7","Chr10","Chr25","Chr27"))
head(input_both_filtered)
input_both_filtered<-subset(input_both, stage!="E7E" & stage!="E7D" &stage!="E9E" &stage!="E9D")
unique(input_both_filtered$chr)
#level the gene names
input_both_filtered$stage <- factor(input_both_filtered$stage,levels=c("E12F","E14F","E16F","D3","D4","D5","D6","D7"))
input_both_filtered$chr <- factor(input_both_filtered$chr,levels=c("chr1","chr2","chr6","chr7","chr10","chr25","chr27"))
head(input_both_filtered)
barplot<-ggplot(input_both, aes(x = chr, y = mean,fill=type))+
  geom_col()+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))
barplot
barplot<-ggplot(input_both, aes(x = id, y = mean,fill=type))+
  geom_col()+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))
barplot
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 50000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 14),
        legend.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45 , hjust = 1 , color = "red", size = 6.5),
        axis.text.y = element_text(color = "black", size = 10))
# make facet
barplot_theme + facet_grid(stage ~.)
barplot<-ggplot(input_both_filtered, aes(x = id, y = mean,fill=type))+
  geom_col()+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 50000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 14),
        legend.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45 , hjust = 1 , color = "red", size = 6.5),
        axis.text.y = element_text(color = "black", size = 10))
# make facet
barplot_theme + facet_grid(stage ~.)
# make facet
barplot_theme + facet_grid(stage ~chr)
input_both_filtered<-subset(input_both, stage!="E7E" & stage!="E7D" &stage!="E9E" &stage!="E9D")
input_both_filtered = input_both_filtered[ rowSums(input_both_filtered)!=0, ]
head(input_both_filtered)
input_both_filtered<-subset(input_both_filtered, mean!=0)
#level the gene names
input_both_filtered$stage <- factor(input_both_filtered$stage,levels=c("E12F","E14F","E16F","D3","D4","D5","D6","D7"))
input_both_filtered$chr <- factor(input_both_filtered$chr,levels=c("chr1","chr2","chr6","chr7","chr10","chr25","chr27"))
barplot<-ggplot(input_both_filtered, aes(x = id, y = mean,fill=type))+
  geom_col()+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))
barplot
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 50000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 14),
        legend.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45 , hjust = 1 , color = "red", size = 6.5),
        axis.text.y = element_text(color = "black", size = 10))
# make facet
barplot_theme + facet_grid(stage ~chr)
View(input_both_filtered)
input_both_filtered<-subset(input_both, stage!="E7E" & stage!="E7D" &stage!="E9E" &stage!="E9D")
input_both_filtered<-subset(input_both_filtered, mean<3)
#level the gene names
input_both_filtered$stage <- factor(input_both_filtered$stage,levels=c("E12F","E14F","E16F","D3","D4","D5","D6","D7"))
input_both_filtered$chr <- factor(input_both_filtered$chr,levels=c("chr1","chr2","chr6","chr7","chr10","chr25","chr27"))
barplot<-ggplot(input_both_filtered, aes(x = id, y = mean,fill=type))+
  geom_col()+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))
barplot
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 50000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 14),
        legend.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45 , hjust = 1 , color = "red", size = 6.5),
        axis.text.y = element_text(color = "black", size = 10))
# make facet
barplot_theme + facet_grid(stage ~chr)
input_both_filtered<-subset(input_both, stage!="E7E" & stage!="E7D" &stage!="E9E" &stage!="E9D")
#level the gene names
input_both_filtered$stage <- factor(input_both_filtered$stage,levels=c("E12F","E14F","E16F","D3","D4","D5","D6","D7"))
input_both_filtered$chr <- factor(input_both_filtered$chr,levels=c("chr1","chr2","chr6","chr7","chr10","chr25","chr27"))
barplot<-ggplot(input_both_filtered, aes(x = id, y = mean,fill=type))+
  geom_col()+ geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 50000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 14),
        legend.text = element_text(color = "black", size = 12),
        axis.text.x = element_text(angle = 45 , hjust = 1 , color = "red", size = 5),
        axis.text.y = element_text(color = "black", size = 10))
# make facet
barplot_theme + facet_grid(stage ~.)
# make facet
barplot_theme + facet_grid(stage ~.) + facet_grid(~chr, scales = "free", switch = "x")
# make facet
barplot_theme + facet_grid(stage ~chr, scales = "free", switch = "x")
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 10000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 14),
        legend.text = element_text(color = "black", size = 12),
        axis.text.x=element_blank(),
        axis.text.y = element_text(color = "black", size = 10))
# make facet
barplot_theme + facet_grid(stage ~chr, scales = "free", switch = "x")
# make facet
barplot_theme + facet_grid(stage ~chr, switch = "x")
# make facet
barplot_theme + facet_grid(stage ~chr, scale="free",switch = "x")
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 10000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 14),
        legend.text = element_text(color = "black", size = 12),
        axis.text.x=element_blank(),
        axis.text.y = element_text(color = "black", size = 6))
# make facet
barplot_theme + facet_grid(stage ~chr, scale="free",switch = "x")
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 10000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 14),
        legend.text = element_text(color = "black", size = 12),
        axis.text.x=element_blank(),
        axis.text.y = element_text(color = "black", size = 8))
# make facet
barplot_theme + facet_grid(stage ~chr, scale="free",switch = "x")
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 10000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 14),
        legend.text = element_text(color = "black", size = 14),
        axis.text.x=element_blank(),
        axis.text.y = element_text(color = "black", size = 8))
# make facet
barplot_theme + facet_grid(stage ~chr, scale="free",switch = "x")
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 10000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 12),
        legend.text = element_text(color = "black", size = 14),
        axis.text.x=element_blank(),
        axis.text.y = element_text(color = "black", size = 8))
# make facet
barplot_theme + facet_grid(stage ~chr, scale="free",switch = "x") + theme(legend.position = "bottom")
barplot_theme<-barplot+
  coord_cartesian(ylim = c(1, 10000))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.title = element_text(color = "blue", size = 12),
        legend.text = element_text(color = "black", size = 10),
        axis.text.x= element_blank(),
        axis.text.y = element_text(color = "black", size = 8))
# make facet
barplot_theme + facet_grid(stage ~chr, scale="free",switch = "x") + theme(legend.position = "bottom")
