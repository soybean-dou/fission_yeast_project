##fission yeast overlapping gene data
##2020-11-12

# install.packages("ape")
# install.packages("stringr")
# install.packages("writexl")
# install.packages("progress")
# install.packages("dplyr")
# install.packages("seqinr")

library(ape)
library(stringr)
library(dplyr)
library(writexl)
library(progress)
library(seqinr)
source('find_overlap_func.r')


#----------------------fission yeast data--------------------------
gene_data<-read.gff("Schizosaccharomyces_pombe_all_chromosomes.gff3")
GSprimer<-read.table("GSprimerPos.tsv", sep = "\t")
Bprimer<-read.table("bPrimerPos2.tsv", sep = "\t")
SERIALprimer<-read.table("serialPrimerPos2.tsv")

colnames(GSprimer)<-c("ID","c-id","up_match","u_ch","up_start","up_end","up_strand","down_match","d_ch","down_start","down_end","down_strand")
colnames(Bprimer)<-c("ID","c-id","match1","1_chr","start1_up","end1_down","strand1",
                      "match2","2_chr","start2_up","end2_down","strand2",
                      "match3","3_chr","start3_up","end3_down","strand3",
                      "match4","4_chr","start4_up","end4_down","strand4",
                      "match5","5_chr","start5_up","end5_down","strand5",
                      "match6","6_chr","start6_up","end6_down","strand6")
#-------------------------------------------------------------
#------------------find CDS site------------------------------

CDS_ch1<-gene_data[0,]; CDS_ch2<-gene_data[0,]; CDS_ch3<-gene_data[0,]; CDS_mit<-gene_data[0,];

for(i in 1:41795){
    if(gene_data[i,3]=="CDS" & 
       !str_detect(gene_data[i,9],'NCRNA') &!str_detect(gene_data[i,9],'TRNA')
       &!str_detect(gene_data[i,9],'RRNA') &!str_detect(gene_data[i,9],'SNORNA')
       &!str_detect(gene_data[i,9],'SNRNA')){
        if(gene_data[i,1]=="I")
            CDS_ch1<-rbind(CDS_ch1,gene_data[i,])
        else if(gene_data[i,1]=="II")
            CDS_ch2<-rbind(CDS_ch2,gene_data[i,])
        else if(gene_data[i,1]=="III")
            CDS_ch3<-rbind(CDS_ch3,gene_data[i,])
        else if(gene_data[i,1]=="mitochondrial")
            CDS_mit<-rbind(CDS_mit,gene_data[i,])
    }
}

overlap_cds_ch1<-find_overlapping_gene(CDS_ch1,CDS_ch1)
overlap_cds_ch2<-find_overlapping_gene(CDS_ch2,CDS_ch2)
overlap_cds_ch3<-find_overlapping_gene(CDS_ch3,CDS_ch3)
overlap_cds_mit<-find_overlapping_gene(CDS_mit,CDS_mit)

CDS_overlap<-mget(ls(pattern = "overlap_cds_"))

#-------------------------------------------------------------
#--------------find CDS-UTR overlapping site-----------------

five_UTR_ch1<-subset(gene_data, seqid == "I" & type =="five_prime_UTR")
five_UTR_ch2<-subset(gene_data, seqid == "II" & type =="five_prime_UTR")
five_UTR_ch3<-subset(gene_data, seqid == "III" & type =="five_prime_UTR")
five_UTR_mit<-subset(gene_data, seqid == "mitochondrial" & type =="five_prime_UTR")

three_UTR_ch1<-subset(gene_data, seqid == "I" & type =="three_prime_UTR")
three_UTR_ch2<-subset(gene_data, seqid == "II" & type =="three_prime_UTR")
three_UTR_ch3<-subset(gene_data, seqid == "III" & type =="three_prime_UTR")
three_UTR_mit<-subset(gene_data, seqid == "mitochondrial" & type =="three_prime_UTR")

overlap_fUTR_ch1<-find_overlapping_gene(CDS_ch1,five_UTR_ch1)
overlap_fUTR_ch2<-find_overlapping_gene(CDS_ch2,five_UTR_ch2)
overlap_fUTR_ch3<-find_overlapping_gene(CDS_ch3,five_UTR_ch3)
overlap_fUTR_mit<-find_overlapping_gene(CDS_mit,five_UTR_mit)

overlap_tUTR_ch1<-find_overlapping_gene(CDS_ch1,three_UTR_ch1)
overlap_tUTR_ch2<-find_overlapping_gene(CDS_ch2,three_UTR_ch2)
overlap_tUTR_ch3<-find_overlapping_gene(CDS_ch3,three_UTR_ch3)
overlap_tUTR_mit<-find_overlapping_gene(CDS_mit,three_UTR_mit)

utr5_overlap<-mget(ls(pattern = "overlap_fUTR_"))
utr3_overlap<-mget(ls(pattern = "overlap_tUTR_"))


###########------save data-----###############
write_xlsx(overlap_cds_ch1,"over_cds_ch1.xlsx")
write_xlsx(overlap_cds_ch2,"over_cds_ch2.xlsx")
write_xlsx(overlap_cds_ch3,"over_cds_ch3.xlsx")
write_xlsx(overlap_cds_mit,"over_cds_mit.xlsx")

write_xlsx(overlap_fUTR_ch1,"over_5utr_ch1.xlsx")
write_xlsx(overlap_fUTR_ch2,"over_5utr_ch2.xlsx")
write_xlsx(overlap_fUTR_ch3,"over_5utr_ch3.xlsx")
write_xlsx(overlap_fUTR_mit,"over_5utr_mit.xlsx")

write_xlsx(overlap_tUTR_ch1,"over_3utr_ch1.xlsx")
write_xlsx(overlap_tUTR_ch2,"over_3utr_ch2.xlsx")
write_xlsx(overlap_tUTR_ch3,"over_3utr_ch3.xlsx")
write_xlsx(overlap_tUTR_mit,"over_3utr_mit.xlsx")
###########------save data-----###############

GS_ch1<-(as.data.frame
         (t(as.data.frame(apply(subset(GSprimer, u_ch == "I"), 
                                1,reverse_negative_strand)))))

GS_ch1 <- type.convert(GS_ch1, as.is = TRUE)


GS_ch2<-(as.data.frame
         (t(as.data.frame(apply(subset(GSprimer, u_ch == "II"),
                                1,reverse_negative_strand)))))
GS_ch2 <- type.convert(GS_ch2, as.is = TRUE)


GS_ch3<-(as.data.frame
         (t(as.data.frame(apply(subset(GSprimer, u_ch == "III"), 
                                1,reverse_negative_strand)))))
GS_ch3 <- type.convert(GS_ch3, as.is = TRUE)

overlap_strand_gs1<-find_overlapping_gene(GS_ch1,CDS_ch1,x_start_name = "up_start", 
                                          x_end_name = "down_end", is_strain = TRUE)

overlap_strand_gs1<-delete_duplicate(overlap_strand_gs1)

overlap_strand_gs2<-find_overlapping_gene(GS_ch2,CDS_ch2, x_start_name = "up_start", 
                                          x_end_name = "down_end", is_strain = TRUE)

overlap_strand_gs3<-find_overlapping_gene(GS_ch3,CDS_ch3)







