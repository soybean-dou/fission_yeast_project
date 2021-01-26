library(RColorBrewer)
source('find_overlap_func.r')

mRNA_ch1_pos<-subset(gene_data, seqid == "I" & type == "mRNA" & strand=="+")
mRNA_ch1_neg<-subset(gene_data, seqid == "I" & type == "mRNA" & strand=="-")

mRNA_ch2_pos<-subset(gene_data, seqid == "II" & type == "mRNA" & strand=="+")
mRNA_ch2_neg<-subset(gene_data, seqid == "II" & type == "mRNA" & strand=="-")

mRNA_ch3_pos<-subset(gene_data, seqid == "III" & type == "mRNA" & strand=="+")
mRNA_ch3_neg<-subset(gene_data, seqid == "III" & type == "mRNA" & strand=="-")


ch1_tendom_5<-find_overlapping_gene(mRNA_ch1_pos, mRNA_ch1_pos, is_same = T, is_type=T)
ch1_tendom_3<-find_overlapping_gene(mRNA_ch1_neg, mRNA_ch1_neg, is_same = T, is_type=T)

ch1_convergent<-find_overlapping_gene(mRNA_ch1_pos, mRNA_ch1_neg, is_same = F, is_type=T)
ch1_divergent<-find_overlapping_gene(mRNA_ch1_neg, mRNA_ch1_pos, is_same = F, is_type=T)

#-----------------------------------------------------------------------------------------

ch2_tendom_5<-find_overlapping_gene(mRNA_ch2_pos, mRNA_ch2_pos, is_same = T, is_type=T)
ch2_tendom_3<-find_overlapping_gene(mRNA_ch2_neg, mRNA_ch2_neg, is_same = T, is_type=T)

ch2_convergent<-find_overlapping_gene(mRNA_ch2_pos, mRNA_ch2_neg, is_same = F, is_type=T)
ch2_divergent<-find_overlapping_gene(mRNA_ch2_neg, mRNA_ch2_pos, is_same = F, is_type=T)

#-----------------------------------------------------------------------------------------

ch3_tendom_5<-find_overlapping_gene(mRNA_ch3_pos, mRNA_ch3_pos, is_same = T, is_type=T)
ch3_tendom_3<-find_overlapping_gene(mRNA_ch3_neg, mRNA_ch3_neg, is_same = T, is_type=T)

ch3_convergent<-find_overlapping_gene(mRNA_ch3_pos, mRNA_ch3_neg, is_same = F, is_type=T)
ch3_divergent<-find_overlapping_gene(mRNA_ch3_neg, mRNA_ch3_pos, is_same = F, is_type=T)

#-----------------------------------------------------------------------------------------
ch1_overlap_type<-mget(ls(pattern = "^ch1_"))
ch2_overlap_type<-mget(ls(pattern = "^ch2_"))
ch3_overlap_type<-mget(ls(pattern = "^ch3_"))

ch1_overlap_info<-c(nrow(ch1_overlap_type[[1]]),nrow(ch1_overlap_type[[2]]),
                    nrow(ch1_overlap_type[[3]]),nrow(ch1_overlap_type[[4]]))
ch2_overlap_info<-c(nrow(ch2_overlap_type[[1]]),nrow(ch2_overlap_type[[2]]),
                    nrow(ch2_overlap_type[[3]]),nrow(ch2_overlap_type[[4]]))
ch3_overlap_info<-c(nrow(ch3_overlap_type[[1]]),nrow(ch3_overlap_type[[2]]),
                    nrow(ch3_overlap_type[[3]]),nrow(ch3_overlap_type[[4]]))

coul <- brewer.pal(4, "Set2") 
barplot(ch1_overlap_info, names=c("convergent","divergent","3-tendom", "5-tendom"), 
        col=coul, border=F, 
        xlab="overlapping type", 
        ylab="number of gene", 
        main="Overlapping gene type of chromosome 1", ylim=c(0,420))

barplot(ch2_overlap_info, names=c("convergent","divergent","3-tandom", "5-tendom"), 
        col=coul, border=F, 
        xlab="overlapping type", 
        ylab="number of gene", 
        main="Overlapping gene type of chromosome 2", ylim=c(0,360))

barplot(ch3_overlap_info, names=c("convergent","divergent","3-tendom", "5-tendom"), 
        col=coul, border=F, 
        xlab="overlapping type", 
        ylab="number of gene", 
        main="Overlapping gene type of chromosome 3", ylim=c(0,170))
#-----------------------------------------------------------------------------------------
overlap_type_data<-data.frame(chromosome1 = ch1_overlap_info, chromosome2=ch2_overlap_info, chromosome3=ch3_overlap_info)
overlap_type_data

layout(matrix(c(1,1,1,2,2), 1,5))
layout.show(2)

barplot(as.matrix(overlap_type_data), beside=T, col=coul, names.arg = c(1,2,3),
        xlab="choromosome", 
        ylab="number of gene", 
        main="Overlapping gene type of S.pombe", ylim=c(0,400))
legend("topright",c("convergent","divergent","3-tandem", "5-tandem"), fill=coul)

barplot(as.matrix(overlap_type_data), col=coul, horiz = F,
        names.arg = c(1,2,3),
        xlab="choromosome", 
        ylab="number of gene", 
        ylim=c(0,800),
        main="Overlapping gene type of S.pombe\n(merged)", width = 60, space = 0.7)
legend("topright",c("convergent","divergent","3-tandem", "5-tandem"), fill=coul,cex = 0.75)



