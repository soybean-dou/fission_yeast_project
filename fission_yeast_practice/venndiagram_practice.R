library(VennDiagram)
library("ggVennDiagram")


grid.newpage()
draw.single.venn(area = 20, 
                 category = "single - colorful", 
                 lty = "blank", 
                 fill = "steelblue", 
                 alpha = 0.5)
grid.newpage()
draw.pairwise.venn(area1 = 25, 
                   area2 = 20, 
                   cross.area = 10, 
                   category = c("diagram A", "diagram B"))


set.seed(20190708)
genes <- paste("gene",1:1000,sep="")
x <- list(
    A = sample(genes,300), 
    B = sample(genes,525), 
    C = sample(genes,440),
    D = sample(genes,350)
)

ggVennDiagram(x, label_alpha = 0)

cds_sample<-read.csv("DATA/cds_data.csv")
utr_sample<-read.csv("DATA/utr_data.csv")
overlap_cds_data<-read.csv("DATA/overlap_cds_data(1).csv")

x<-list(
    current = overlap_cds_data[,1],
    before = cds_sample[,1]
)
ggVennDiagram(x, label_alpha = 0)


x<-list(
    current = overlap_strand_UTR_all[,1],
    before = utr_sample[,1]
)
ggVennDiagram(x, label_alpha = 0)

x<-list(
    current = overlap_strand_3UTR[,1],
    before = utr_sample[,1]
)
ggVennDiagram(x, label_alpha = 0)