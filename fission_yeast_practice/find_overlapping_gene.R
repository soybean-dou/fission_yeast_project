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
GSprimer<-read.table("GSPrimerPos.tsv", sep = "\t")
Bprimer<-read.table("bPrimerPos2.tsv", sep = "\t")
SERIALprimer<-read.table("serialPrimerPos2.tsv")

colnames(GSprimer)<-c("ID","c-id","match1","chr1","start1","end1","strand1",
                      "match2","chr2","start2","end2","strand2")

colnames(Bprimer)<-c("ID","c-id","match1","chr1","start1","end1","strand1",
                      "match2","chr2","start2","end2","strand2",
                      "match3","chr3","start3","end3","strand3",
                      "match4","chr4","start4","end4","strand4",
                      "match5","chr5","start5","end5","strand5",
                      "match6","chr6","start6","end6","strand6")
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

overlap_cds_ch1<-find_overlapping_gene(CDS_ch1,CDS_ch1,is_same = TRUE)
overlap_cds_ch2<-find_overlapping_gene(CDS_ch2,CDS_ch2,is_same = TRUE)
overlap_cds_ch3<-find_overlapping_gene(CDS_ch3,CDS_ch3,is_same = TRUE)
overlap_cds_mit<-find_overlapping_gene(CDS_mit,CDS_mit,is_same = TRUE)

CDS_site<-mget(ls(pattern = "CDS_"))
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


five_UTR<-mget(ls(pattern = "five_UTR_"))
three_UTR<-mget(ls(pattern = "three_UTR_"))
utr5_overlap<-mget(ls(pattern = "overlap_fUTR_"))
utr3_overlap<-mget(ls(pattern = "overlap_tUTR_"))


###########------save data-----###############
for(i in 1:3){
    file_name <- paste0("CDS_overlap_", toString(i),".xlsx")
    write_xlsx(CDS_overlap[i],file_name)
}

###########------GS primer data overlap site finding code-----###############

GS_ch1<-type.convert((as.data.frame(t(as.data.frame(
                    apply(subset(GSprimer, chr1 == "I"), 
                          1,reverse_negative_strand,
                          up_start="start1", up_end="end1", 
                          down_start="start2", down_end="end2",
                          up_strand="strand1", down_strand="strand2"))))), 
                    as.is = TRUE)

GS_ch2<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(GSprimer, chr1 == "II"), 
          1,reverse_negative_strand,
          up_start="start1", up_end="end1", 
          down_start="start2", down_end="end2",
          up_strand="strand1", down_strand="strand2"))))), 
    as.is = TRUE)

GS_ch3<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(GSprimer, chr1 == "III"), 
          1,reverse_negative_strand,
          up_start="start1", up_end="end1", 
          down_start="start2", down_end="end2",
          up_strand="strand1", down_strand="strand2"))))), 
    as.is = TRUE)

overlap_strand_gs1<-find_overlapping_gene(GS_ch1,CDS_site$CDS_ch1,
                                          x_start_name = "start1", x_end_name = "end2", 
                                          is_strain = TRUE)

overlap_strand_gs2<-find_overlapping_gene(GS_ch2,CDS_site$CDS_ch2,
                                          x_start_name = "start1", x_end_name = "end2", 
                                          is_strain = TRUE)

overlap_strand_gs3<-find_overlapping_gene(GS_ch3,CDS_site$CDS_ch3,
                                          x_start_name = "start1", x_end_name = "end2", 
                                          is_strain = TRUE)

GS<-mget(ls(pattern = "GS_"))
GS_overlap<-mget(ls(pattern = "overlap_strand_gs"))
write_xlsx(overlap_strand_gs1,"GS_overlap_1.xlsx")

for(i in 1:3){
    file_name <- paste0("GS_overlap_", toString(i),".xlsx")
    write_xlsx(GS_overlap[i],file_name)
}
###########------block PCR data overlap site finding code-----###############

#--prepare data : negative strand modification

B_ch1<-type.convert((as.data.frame(t(as.data.frame(
        apply(subset(Bprimer, chr1 == "I"), 
                1,reverse_negative_strand,
                up_start="start3", up_end="end4", 
                down_start="start4", down_end="end4",
                up_strand="strand3", down_strand="strand4"))))), as.is = TRUE)

B_ch2<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(Bprimer, chr1 == "II"), 
          1,reverse_negative_strand,
          up_start="start1", up_end="end1", 
          down_start="start2", down_end="end2",
          up_strand="strand1", down_strand="strand2"))))), as.is = TRUE)

B_ch3<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(Bprimer, chr1 == "III"), 
          1,reverse_negative_strand,
          up_start="start1", up_end="end1", 
          down_start="start2", down_end="end2",
          up_strand="strand1", down_strand="strand2"))))), as.is = TRUE)

#--finding overlapping site
overlap_strand_B1<-find_overlapping_gene(B_ch1,CDS_site$CDS_ch1,
                                          x_start_name = "start3", x_end_name = "end4", 
                                          is_strain = TRUE)

overlap_strand_B2<-find_overlapping_gene(B_ch2,CDS_site$CDS_ch2,
                                         x_start_name = "start3", x_end_name = "end4", 
                                         is_strain = TRUE)

overlap_strand_B3<-find_overlapping_gene(B_ch3,CDS_site$CDS_ch3,
                                         x_start_name = "start3", x_end_name = "end4", 
                                         is_strain = TRUE)

Block<-mget(ls(pattern = "B_"))
B_overlap<-mget(ls(pattern = "overlap_strand_B"))

for(i in 1:3){
    file_name <- paste0("B_overlap_", toString(i),".xlsx")
    write_xlsx(B_overlap[i],file_name)
}
