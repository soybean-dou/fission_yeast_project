
    
overlap_info<-data.frame(nrow(CDS_overlap$ch1),nrow(CDS_overlap$ch2),nrow(CDS_overlap$ch3))


dev.new()

barplot(as.matrix(overlap_info), names.arg = c("1","2","3"),
        col=c('red','green','blue'), border=F, 
        xlab="chromosome", 
        ylab="number of overlapping site", 
        ylim=c(0,80),
        width = 0.8, space=0.8,
        main="Overlap between deletion model and CDS")
