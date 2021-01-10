classificate_overlap_type<-function(x,y){
}

source('find_overlap_func.r')

mRNA_ch1<-subset(gene_data, seqid == "I" & type == "mRNA")
mRNA_ch2<-subset(gene_data, seqid == "II" & type == "mRNA")
mRNA_ch3<-subset(gene_data, seqid == "III" & type == "mRNA")

mRNA_ch1_pos_strand<-subset(mRNA_ch1, strand=="+")
mRNA_ch1_neg_strand<-subset(mRNA_ch1, strand=="-")

tendom_5<-find_overlapping_gene(mRNA_ch1_pos_strand, mRNA_ch1_pos_strand, is_same = T)
tendom_3<-find_overlapping_gene(mRNA_ch1_neg_strand, mRNA_ch1_neg_strand, is_same = T)

convergent<-find_overlapping_gene(mRNA_ch1_pos_strand, mRNA_ch1_neg_strand, is_same = F)
divergent<-find_overlapping_gene(mRNA_ch1_neg_strand, mRNA_ch1_pos_strand, is_same = F)
