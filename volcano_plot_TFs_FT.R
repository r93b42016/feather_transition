setwd("C:\\Users\\r93b4\\Desktop\\Myprojects\\feather_transition_keratinocyte\\HINT_ATAC_opt")
getwd()

     library(ggplot2)
     library("RColorBrewer")
     library(ggrepel)

#input HINT-ATAC output.txt (differential_statistics.txt)
tfs <- read.table(file.choose(), header=T, sep="\t")
head(tfs)
#pick one below!!
tfs$Expression = ifelse(tfs$P_values < 0.06,
                             ifelse(tfs$TF_Activity > 0 ,'juvenile_up','downy_up'),
                             'None')
#overlap with color TFs
tfs$Label <- ifelse(tfs$Expression != "None",tfs$Motif,'')
#open a big window to let this plot work (for ggrepel)
windows(9,6)
#plot
ggplot(data=tfs, aes(x=TF_Activity, y=-log10(P_values), col=Expression, label=Label)) +
       geom_point(size=2) +
       geom_label_repel(max.overlaps =3) +
       scale_color_manual(values=c("blue","red", "grey"))+
       theme_bw()+
       theme(legend.text = element_text(size = 12),axis.title = element_text(size = 16))
     #open a big window to let this plot work (for ggrepel)
     windows(9,6)
     #plot
     ggplot(data=tfs, aes(x=TF_Activity, y=-log10(P_values), col=Expression, label=Label)) +
       geom_point(size=2) +
       geom_label_repel(max.overlaps =1000) +
       scale_color_manual(values=c("blue","red", "grey"))+
       theme_bw()+
       theme(legend.text = element_text(size = 12),axis.title = element_text(size = 16))
     #plot
     ggplot(data=tfs, aes(x=TF_Activity, y=-log10(P_values), col=Expression, label=Label)) +
       geom_point(size=2) +
       geom_text_repel(max.overlaps =1000) +
       scale_color_manual(values=c("blue","red", "grey"))+
       theme_bw()+
       theme(legend.text = element_text(size = 12),axis.title = element_text(size = 16))
     #open a big window to let this plot work (for ggrepel)
     windows(9,6)
     #plot
     ggplot(data=tfs, aes(x=TF_Activity, y=-log10(P_values), col=Expression, label=Label)) +
       geom_point(size=2) +
       geom_text_repel(max.overlaps =1000) +
       scale_color_manual(values=c("blue","red", "grey"))+
       theme_bw()+
       theme(legend.text = element_text(size = 12),axis.title = element_text(size = 16))
     