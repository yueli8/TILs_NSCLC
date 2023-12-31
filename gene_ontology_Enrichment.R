library(ggplot2)
setwd("~/gse184053/data/data_finally")
a<-read.table("Enrichment_plot.txt", head=TRUE, sep="\t")
S1<- ggplot(a, aes(x= fold_Enrichment, y=reorder(GO, fold_Enrichment), size=counts,fill=FDR)) + geom_point(shape = 21) +theme_bw() +theme()
S1=S1+ scale_fill_continuous(low = '#d90424', high = '#374a89')+scale_x_continuous(
  labels = scales::number_format(accuracy = 0.1))+ theme(axis.text.y = element_text(size = 16, family = "Arial", color = "black", face = "bold"))+xlab("fold_Enrichment")+scale_x_continuous(breaks = seq(0, 45, by = 10))
S1