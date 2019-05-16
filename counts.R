#http://ggplot2.tidyverse.org/reference/theme.html
#http://bioinfo.umassmed.edu/bootstrappers/bootstrappers-courses/pastCourses/rCourse_2014-09/Session_6/Session_6.html
ggplot(Satb2, aes(x = v1, y = COUNTS))+geom_bar(stat="identity")+geom_errorbar(aes(ymin = COUNTS-SD, ymax = COUNTS+SD,width = 0.2))+geom_line()+ theme(panel.background = element_rect(fill = "white", colour = "grey50"))

ggplot(Satb2, aes(x = v1, y = COUNTS))+geom_bar(stat="identity")+geom_errorbar(aes(ymin = COUNTS-SD, ymax = COUNTS+SD,width = 0.2))+geom_line()+ theme(panel.background = element_rect(fill = "white", colour = "grey50"))






ggplot(Satb2, aes(x = v1, y = COUNTS))+geom_bar(stat="identity")+geom_errorbar(aes(ymin = COUNTS-SD, ymax = COUNTS+SD,width = 0.2))+geom_line()+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white",axis.line = element_line()))

p1 = ggplot(Satb2, aes(x = v1, y = COUNTS))+geom_bar(stat="identity")+geom_errorbar(aes(ymin =COUNTS-SD, ymax = COUNTS+SD,width = 0.2))+geom_line()


p1+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white"))

p1+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white"),axis.line = element_line(size = 0.5, colour = "black"))

p2 = p1+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white"),axis.line = element_line(size = 0.5, colour = "black"))

p2+ labs( x = "STAGES",y ="COUNTS",colour ="Cylinders",shape ="Transmission")
p3 = p2+ labs( title = "Satb2",x = "STAGES",y ="COUNTS",colour ="Cylinders",shape ="Transmission")
########################################################################

idk2 <- read.table("~/zebra_stage_RNA/C.txt","\t", header=T,row.names=1)
data2 = data.matrix(idk2[,1:ncol(idk2)])


my_palette <- colorRampPalette(c("blue","white","green"))(n =299)
heatmap(data2, Colv=NA, Rowv=NA, col=my_palette)

heatmap(data2, Colv=NA, Rowv=NA, col=my_palette)
heatmap.2(data2, Rowv = NA, Colv = NA, trace = "none", col = my_palette2,
          key.xlab = "cpm" ,key.ylab = NULL,main = NULL, xlab = "STAGES", ylab = "COUNTS")

my_palette2 <- colorRampPalette(c("white","red", "green"))(n = 299)
heatmap.2(data2, Rowv = NA, Colv = NA, trace = "none", col = my_palette2)


Breaks1=seq(0,400,50)
heatmap.2(data2, Rowv = NA, Colv = NA, trace = "none", col = my_palette2, breaks = Breaks1)

heatmap.2(data2, Rowv = FALSE, 
          Colv = FALSE, col = my_palette2 ,
          density.info="none",trace="none")
##############################################################################

heatmap(data2, Colv=NA, Rowv=NA, col=my_palette, key= TRUE)
my_palette <- colorRampPalette(c("blue","white","green"))
heatmap(data2, Colv=NA, Rowv=NA, col=my_palette,key=TRUE)
heatmap.2(data2, Colv=NA, Rowv=NA,key=T,keysize=1)
heatmap(data2, Colv=NA, Rowv=NA, col=my_palette, ColSideColors, RowSideColors)

#GG
library(ggplot2)
library(reshape2)
data2.m = melt(data2)

data2.m <- ddply(data2.m, .(variable), transform, rescale = rescale(value))







Col=c("blue","white","green")
Colors=colorRampPalette(Colors)(10)

heatmap.2(data2, Rowv = FALSE, 
          Colv = FALSE, col = Col ,
          density.info="none",trace="none",
          breaks = Breaks)

library(RColorBrewer)
Colors=brewer.pal(11,"Spectral")
Colors=colorRampPalette(Colors)(50)
heatmap.2(data2, Rowv = FALSE, 
          Colv = FALSE, col = Colors ,
          density.info="none",trace="none")

Breaks=seq(0,100,10)
Colors1=rev(brewer.pal(6,"Spectral"))
heatmap.2(data2, Rowv = FALSE, 
          Colv = FALSE, col = Colors1 ,
          density.info="none",trace="none",
          breaks=Breaks)

Colors=c("blue","white","green")
Colors=colorRampPalette(Colors)(300)
Breaks=seq(0,400,50)
heatmap.2(data2, Rowv = FALSE, 
          Colv = FALSE, col = Colors ,
          density.info="none",trace="none",
          breaks=Breaks)

Colors=c(blue)
Colors=colorRampPalette(Colors)(300)
#Breaks=seq(0,400,50)
heatmap.2(data2, Rowv = FALSE, 
          Colv = FALSE, col = Colors ,
          density.info="none",trace="none")

Breaks=seq(0,400,50)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
heatmap.2(data2, Rowv = FALSE, 
          Colv = FALSE, col = Colors ,
          density.info="none",trace="none",
          breaks = Breaks)
################################################################################################

idk2 <- read.table("~/zebra_stage_RNA/cluster_test/cpm_log.nrm","\t", header=T,row.names=1)

data = idk2[-1,-c(1,2)] View(data2)
cpm_log = data.matrix(data[,1:ncol(data)])

my_palette <- colorRampPalette(c("blue","white","red"))(n =299)
heatmap.2(cpm_log, Rowv = NA, Colv = NA, trace = "none", col = my_palette,
            +           key.xlab = "cpm_log" ,key.ylab = NULL ,main = "cpm_log", xlab = "STAGES", ylab = "log_cpm")



col = colorRampPalette(c("blue","white","yellow"))(n =299)
heatmap.2(cpm_log, Rowv = NA, Colv = NA, trace = "none", col = col,
          +           key.xlab = "cpm" ,key.ylab = NULL , main = "cpm_log", xlab = "STAGES", ylab = "log_cpm")

my_palette <- colorRampPalette(c("blue","white","red"))(n =299)
heatmap.2(cpm_log, Rowv = NA, Colv = NA, trace = "none", col = my_palette,
          key.xlab = "log_cpm", key.ylab = NULL, xlab = "STAGES", ylab = "Log_CPM", main = "Log_cpm")

