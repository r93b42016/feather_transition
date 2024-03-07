library(ggplot2)

setwd("C:\\Users\\r93b4\\Desktop\\Myprojects\\feather_transition_keratinocyte\\Figures\\Violin_plot")
getwd()

input_1 <- read.table(file.choose(), header=TRUE,sep = ",")
head(input_1)

#level the gene names
input_1$Stage <- factor(input_1$Stage,levels=c("E14","P6","P7"))
input_1$Locus <- factor(input_1$Locus,levels=c("Chr25FK","Chr25FL","Chr25Scale","Chr25Claw","Chr27FK","Chr10FK","Chr2FK"))

ggplot() + geom_violin(aes(c(cormat), c(cormat)))

# Set trim argument to FALSE and assign basic color
p <- ggplot(input_1, aes(x=Stage, y=UQTPM, color=Stage)) + 
  geom_violin(trim=FALSE)
p

# present dots
pdots <- p + geom_jitter(shape=16, position=position_jitter(0.2))+ theme_classic()

#set log value
logp <- pdots + scale_y_continuous(trans = 'log10')
logp

### Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# make facet
logp + facet_wrap(~ Locus) 
