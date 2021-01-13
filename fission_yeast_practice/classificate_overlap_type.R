classificate_overlap_type<-function(x,y){
}

library(RColorBrewer)
source('find_overlap_func.r')

mRNA_ch1_pos<-subset(gene_data, seqid == "I" & type == "mRNA" & strand=="+")
mRNA_ch1_neg<-subset(gene_data, seqid == "I" & type == "mRNA" & strand=="-")

mRNA_ch2_pos<-subset(gene_data, seqid == "II" & type == "mRNA" & strand=="+")
mRNA_ch2_neg<-subset(gene_data, seqid == "II" & type == "mRNA" & strand=="-")

mRNA_ch3_pos<-subset(gene_data, seqid == "III" & type == "mRNA" & strand=="+")
mRNA_ch3_neg<-subset(gene_data, seqid == "III" & type == "mRNA" & strand=="-")


ch1_tendom_5<-find_overlapping_gene(mRNA_ch1_pos, mRNA_ch1_pos, is_same = T)
ch1_tendom_3<-find_overlapping_gene(mRNA_ch1_neg, mRNA_ch1_neg, is_same = T)

ch1_convergent<-find_overlapping_gene(mRNA_ch1_pos, mRNA_ch1_neg, is_same = F)
ch1_divergent<-find_overlapping_gene(mRNA_ch1_neg, mRNA_ch1_pos, is_same = F)

#-----------------------------------------------------------------------------------------

ch2_tendom_5<-find_overlapping_gene(mRNA_ch2_pos, mRNA_ch2_pos, is_same = T)
ch2_tendom_3<-find_overlapping_gene(mRNA_ch2_neg, mRNA_ch2_neg, is_same = T)

ch2_convergent<-find_overlapping_gene(mRNA_ch2_pos, mRNA_ch2_neg, is_same = F)
ch2_divergent<-find_overlapping_gene(mRNA_ch2_neg, mRNA_ch2_pos, is_same = F)

#-----------------------------------------------------------------------------------------

ch3_tendom_5<-find_overlapping_gene(mRNA_ch3_pos, mRNA_ch3_pos, is_same = T)
ch3_tendom_3<-find_overlapping_gene(mRNA_ch3_neg, mRNA_ch3_neg, is_same = T)

ch3_convergent<-find_overlapping_gene(mRNA_ch3_pos, mRNA_ch3_neg, is_same = F)
ch3_divergent<-find_overlapping_gene(mRNA_ch3_neg, mRNA_ch3_pos, is_same = F)

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

coul <- brewer.pal(5, "Set2") 
barplot(ch1_overlap_info, names=c("convergent","divergent","3-tendom", "5-tendom"), 
        col=coul, border=F, 
        xlab="overlapping type", 
        ylab="number of gene", 
        main="Overlapping gene type of chromosome 1", ylim=c(0,420))

barplot(ch2_overlap_info, names=c("convergent","divergent","3-tendom", "5-tendom"), 
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


