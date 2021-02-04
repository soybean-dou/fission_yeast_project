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
