library(ggplot2)
setwd("~/zhangyi_rat/rat")
a<-read.table("Enrichment_plot.txt", head=TRUE, sep="\t")
S1<- ggplot(a, aes(x= fold_Enrichment, y=reorder(GO, fold_Enrichment), size=counts,fill=FDR)) + geom_point(shape = 21) +theme_bw() +theme()
#S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab")
#S1<-S1 + theme(axis.title.x =element_text(face="bold", size=20), axis.title.y=element_text(face="bold",size=10))
#S1


S1=S1+ scale_fill_continuous(low = '#d90424', high = '#374a89')+scale_x_continuous(
  labels = scales::number_format(accuracy = 0.1))+ theme(axis.title.y = element_blank())+xlab("fold_Enrichment")+scale_x_continuous(breaks = seq(0, 100, by = 10))
S1