##fission yeast overlapping gene data
##2020-11-12

# install.packages("ape")
# install.packages("stringr")
# install.packages("writexl")
# install.packages("progress")
# install.packages("dplyr")
# install.packages("seqinr")
# install.packages("readxl")

library(ape)
library(stringr)
library(dplyr)
library(writexl)
library(xlsx)
library(readxl)
library(progress)
library(seqinr)
library(RColorBrewer)
source('find_overlap_func.r')


#----------------------fission yeast data--------------------------
gene_data<-read.gff("Schizosaccharomyces_pombe_all_chromosomes.gff3")

ch1<-read.table("DATA/chromosome1.cds.coords.tsv")
colnames(ch1)<-c("ID","start","end","strand")

ch2<-read.table("DATA/chromosome2.cds.coords.tsv")
colnames(ch2)<-c("ID","start","end","strand")

ch3<-read.table("DATA/chromosome3.cds.coords.tsv")
colnames(ch3)<-c("ID","start","end","strand")

CDS<-mget(ls(pattern = "ch"))

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
colnames(SERIALprimer)<-c("ID","c-id","sequence1","match1","chr1","start1","end1","strand1",
                      "sequence2","match2","chr2","start2","end2","strand2",
                      "match3","chr3","start3","end3","strand3",
                      "match4","chr4","start4","end4","strand4")


Bprimer<-transform(Bprimer, start5 = as.numeric(start5), end5 = as.numeric(end5), 
                   start6 = as.numeric(start6),end6 = as.numeric(end6))
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

CDS<-mget(ls(pattern = "CDS_"))
CDS_overlap<-mget(ls(pattern = "overlap_cds_"))

#-------------------------------------------------------------
#--------------find CDS-UTR overlapping site-----------------

ch1<-subset(gene_data, seqid == "I" & type =="five_prime_UTR")
ch2<-subset(gene_data, seqid == "II" & type =="five_prime_UTR")
ch3<-subset(gene_data, seqid == "III" & type =="five_prime_UTR")
mit<-subset(gene_data, seqid == "mitochondrial" & type =="five_prime_UTR")

five_UTR<-mget(ls(pattern = "ch"))

ch1<-subset(gene_data, seqid == "I" & type =="three_prime_UTR")
ch2<-subset(gene_data, seqid == "II" & type =="three_prime_UTR")
ch3<-subset(gene_data, seqid == "III" & type =="three_prime_UTR")
mit<-subset(gene_data, seqid == "mitochondrial" & type =="three_prime_UTR")

three_UTR<-mget(ls(pattern = "ch"))

overlap_fUTR_ch1<-find_overlapping_gene(CDS_ch1,five_UTR_ch1)
overlap_fUTR_ch2<-find_overlapping_gene(CDS_ch2,five_UTR_ch2)
overlap_fUTR_ch3<-find_overlapping_gene(CDS_ch3,five_UTR_ch3)
overlap_fUTR_mit<-find_overlapping_gene(CDS_mit,five_UTR_mit)

utr5_overlap<-mget(ls(pattern = "overlap_fUTR_"))

overlap_tUTR_ch1<-find_overlapping_gene(CDS_ch1,three_UTR_ch1)
overlap_tUTR_ch2<-find_overlapping_gene(CDS_ch2,three_UTR_ch2)
overlap_tUTR_ch3<-find_overlapping_gene(CDS_ch3,three_UTR_ch3)
overlap_tUTR_mit<-find_overlapping_gene(CDS_mit,three_UTR_mit)

utr3_overlap<-mget(ls(pattern = "overlap_tUTR_"))

###########------save data-----###############
for(i in 1:3){
    file_name <- paste0("CDS_overlap_", toString(i),".xlsx")
    write_xlsx(CDS_overlap[i],file_name)
}

#--prepare data : primer data modification

ch1<-type.convert((as.data.frame(t(as.data.frame(
                    apply(subset(GSprimer, chr1 == "I"), 
                          1,reverse_negative_strand,
                          up_start="start1", up_end="end1", 
                          down_start="start2", down_end="end2",
                          up_strand="strand1", down_strand="strand2"))))), 
                    as.is = TRUE)

ch2<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(GSprimer, chr1 == "II"), 
          1,reverse_negative_strand,
          up_start="start1", up_end="end1", 
          down_start="start2", down_end="end2",
          up_strand="strand1", down_strand="strand2"))))), 
    as.is = TRUE)

ch3<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(GSprimer, chr1 == "III"), 
          1,reverse_negative_strand,
          up_start="start1", up_end="end1", 
          down_start="start2", down_end="end2",
          up_strand="strand1", down_strand="strand2"))))), 
    as.is = TRUE)

GS<-mget(ls(pattern = "ch"))

ch1<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(Bprimer, chr1 == "I"), 
          1,reverse_negative_strand,
          up_start="start3", up_end="end3", 
          down_start="start4", down_end="end4",
          up_strand="strand3", down_strand="strand4"))))), as.is = TRUE)


ch1 <- na.omit(B_ch1)

ch2<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(Bprimer, chr1 == "II"), 
          1,reverse_negative_strand,
          up_start="start3", up_end="end3", 
          down_start="start4", down_end="end4",
          up_strand="strand3", down_strand="strand4"))))), as.is = TRUE)

ch3<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(Bprimer, chr1 == "III"), 
          1,reverse_negative_strand,
          up_start="start3", up_end="end3", 
          down_start="start4", down_end="end4",
          up_strand="strand3", down_strand="strand4"))))), as.is = TRUE)

Block<-mget(ls(pattern = "ch"))

ch1<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(SERIALprimer, chr1 == "I"), 
          1,reverse_negative_strand,
          up_start="start1", up_end="end1", 
          down_start="start2", down_end="end2",
          up_strand="strand1", down_strand="strand2"))))), as.is = TRUE)

ch2<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(SERIALprimer, chr1 == "II"), 
          1,reverse_negative_strand,
          up_start="start1", up_end="end1", 
          down_start="start2", down_end="end2",
          up_strand="strand1", down_strand="strand2"))))), as.is = TRUE)

ch3<-type.convert((as.data.frame(t(as.data.frame(
    apply(subset(SERIALprimer, chr1 == "III"), 
          1,reverse_negative_strand,
          up_start="start1", up_end="end1", 
          down_start="start2", down_end="end2",
          up_strand="strand1", down_strand="strand2"))))), as.is = TRUE)

Serial<-mget(ls(pattern = "ch"))


###########------overlap site finding code-----###############

ch1<-find_overlapping_gene(GS$ch1, CDS$ch1,
                            x_start_name = "start1", x_end_name = "end2", 
                            is_strain = TRUE)

ch2<-find_overlapping_gene(GS$ch2,CDS$ch2,
                            x_start_name = "start1", x_end_name = "end2", 
                            is_strain = TRUE)

ch3<-find_overlapping_gene(GS$ch3,CDS$ch3,
                          x_start_name = "start1", x_end_name = "end2", 
                          is_strain = TRUE)

GS_CDS<-mget(ls(pattern = "ch"))


ch1<-find_overlapping_gene(Block$ch1,CDS$ch1,
                            x_start_name = "start3", x_end_name = "end4", 
                            is_strain = TRUE, b_primer=T)

ch2<-find_overlapping_gene(Block$ch2,CDS$ch2,
                            x_start_name = "start3", x_end_name = "end4", 
                            is_strain = TRUE, b_primer=T)

ch3<-find_overlapping_gene(Block$ch3,CDS$ch3,
                            x_start_name = "start3", x_end_name = "end4", 
                            is_strain = TRUE, b_primer=T)

B_CDS<-mget(ls(pattern = "ch"))

ch1<-find_overlapping_gene(Serial$ch1,CDS$ch1,
                                         x_start_name = "start1", x_end_name = "end2", 
                                         is_strain = TRUE)

ch2<-find_overlapping_gene(Serial$ch2,CDS$ch2,
                                         x_start_name = "start1", x_end_name = "end2", 
                                         is_strain = TRUE)

ch3<-find_overlapping_gene(Serial$ch3,CDS$ch3,
                                         x_start_name = "start1", x_end_name = "end2", 
                                         is_strain = TRUE)

S_CDS<-mget(ls(pattern = "ch1"))


overlap_strand_data<-c(nrow(GS_overlap[[1]])+nrow(B_overlap[[1]])+nrow(S_overlap[[1]]),
                            nrow(GS_overlap[[2]])+nrow(B_overlap[[2]])+nrow(S_overlap[[2]]),
                            nrow(GS_overlap[[3]])+nrow(B_overlap[[3]])+nrow(S_overlap[[3]]))

overlap_strand_data<-as.data.frame(overlap_strand_data)
overlap_strand_data

layout(matrix(c(1,2),1,2))
layout.show(1)

coul <- brewer.pal(3, "Set1") 
barplot(overlap_strand_data, names.arg = c("1","2","3"),
        col=coul, border=F, 
        xlab="chromosome", 
        ylab="number of overlapping site", 
        ylim=c(0,80),
        width = 0.8, space=0.8,
        main="Overlap between deletion model and CDS")

overlap_strand_data_mg<-data.frame(chromosome1=c(overlap_strand_data[1], nrow(GS_overlap[[1]]),nrow(B_overlap[[1]]),nrow(S_overlap[[1]])), 
           chromosome2=c(overlap_strand_data[2], nrow(GS_overlap[[2]]),nrow(B_overlap[[2]]),nrow(S_overlap[[2]])),
           chromosome3=c(overlap_strand_data[3],nrow(GS_overlap[[3]]),nrow(B_overlap[[3]]),nrow(S_overlap[[3]])))
overlap_strand_data_mg



coul <- brewer.pal(4, "Set1") 
barplot(as.matrix(overlap_strand_data_mg), beside=T,names.arg = c("1","2","3"),
        col=coul, 
        xlab="chromosome", 
        ylab="number of overlapping site", ylim=c(0,90),
        main="Overlap between \ndeletion model and CDS")
legend("topright",c("all","GS","Block","SERIAL"), fill=coul,cex = 0.75)

#--------------------------------------
overlap_matrix_ch1<-data.matrix(overlap_matrix$overlap_matrix_ch1)
overlap_heatmap <- heatmap(overlap_matrix_ch1, Rowv=NA, Colv=NA, col=brewer.pal(9, "Blues"), scale="row", margin=c(5,10))


#--------------------------------------
ch1<-find_overlapping_gene(GS$ch1, five_UTR$ch1,
                           x_start_name = "start1", x_end_name = "end2", 
                          is_strain = TRUE)

ch2<-find_overlapping_gene(GS$ch2,five_UTR$ch2,
                          x_start_name = "start1", x_end_name = "end2", 
                          is_strain = TRUE)

ch3<-find_overlapping_gene(GS$ch3,five_UTR$ch3,
                            x_start_name = "start1", x_end_name = "end2", 
                            is_strain = TRUE)

GS_5utr<-mget(ls(pattern = "ch"))

ch1<-find_overlapping_gene(Block$ch1,five_UTR$ch1,
                            x_start_name = "start3", x_end_name = "end4", 
                            is_strain = TRUE)

ch2<-find_overlapping_gene(Block$ch2,five_UTR$ch2,
                            x_start_name = "start3", x_end_name = "end4", 
                            is_strain = TRUE)

ch3<-find_overlapping_gene(Block$ch3,five_UTR$ch3,
                            x_start_name = "start3", x_end_name = "end4", 
                            is_strain = TRUE)

B_5utr<-mget(ls(pattern = "ch"))

ch1<-find_overlapping_gene(Serial$ch1,five_UTR$ch1,
                                         x_start_name = "start1", x_end_name = "end2", 
                                         is_strain = TRUE)

ch2<-find_overlapping_gene(Serial$ch2,five_UTR$ch2,
                                         x_start_name = "start1", x_end_name = "end2", 
                                         is_strain = TRUE)

ch3<-find_overlapping_gene(Serial$ch3,five_UTR$ch3,
                                         x_start_name = "start1", x_end_name = "end2", 
                                         is_strain = TRUE)

S_5utr<-mget(ls(pattern = "ch"))
#----------------
overlap_strand_UTR<-mget(ls(pattern = "overlap_strand_UTR"))
overlap_strand_UTR_ch1<-rbind(overlap_strand_UTR_B1,overlap_strand_UTR_gs1,overlap_strand_UTR_S1)
overlap_strand_UTR_ch2<-rbind(overlap_strand_UTR_B2,overlap_strand_UTR_gs2,overlap_strand_UTR_S2)
overlap_strand_UTR_ch3<-rbind(overlap_strand_UTR_B3,overlap_strand_UTR_gs3,overlap_strand_UTR_S3)
overlap_strand_UTR_all<-rbind(overlap_strand_UTR_ch1,overlap_strand_UTR_ch2,overlap_strand_UTR_ch3)


#--------------------------------------
ch1<-find_overlapping_gene(GS$ch1, three_UTR$ch1,
                          x_start_name = "start1", x_end_name = "end2", 
                                              is_strain = TRUE)

ch2<-find_overlapping_gene(GS$ch2,three_UTR$ch2,
                          x_start_name = "start1", x_end_name = "end2", 
                          is_strain = TRUE)

ch3<-find_overlapping_gene(GS$ch3,three_UTR$ch3,
                          x_start_name = "start1", x_end_name = "end2", 
                          is_strain = TRUE)

GS_3utr<-mget(ls(pattern = "ch"))

ch1<-find_overlapping_gene(Block$ch1,three_UTR$ch1,
                             x_start_name = "start3", x_end_name = "end4", 
                             is_strain = TRUE)

ch2<-find_overlapping_gene(Block$ch2,three_UTR$ch2,
                             x_start_name = "start3", x_end_name = "end4", 
                             is_strain = TRUE)

ch3<-find_overlapping_gene(Block$ch3,three_UTR$ch3,
                             x_start_name = "start3", x_end_name = "end4", 
                             is_strain = TRUE)

B_3utr<-mget(ls(pattern = "ch"))

ch1<-find_overlapping_gene(Serial$ch1,three_UTR$ch1,
                             x_start_name = "start1", x_end_name = "end2", 
                             is_strain = TRUE)

ch2<-find_overlapping_gene(Serial$ch2,three_UTR$ch2,
                             x_start_name = "start1", x_end_name = "end2", 
                             is_strain = TRUE)

ch3<-find_overlapping_gene(Serial$ch3,three_UTR$ch3,
                             x_start_name = "start1", x_end_name = "end2", 
                             is_strain = TRUE)

S_3utr<-mget(ls(pattern = "ch"))
#----------------
ch1<-subset(gene_data, seqid == "I" & type =="gene")
ch1<-as.data.frame(t(as.data.frame(apply(ch1, 1,modify_id))))

ch2<-subset(gene_data, seqid == "II" & type =="gene")
ch2<-as.data.frame(t(as.data.frame(apply(ch2, 1,modify_id))))

ch3<-subset(gene_data, seqid == "III" & type =="gene")
ch3<-as.data.frame(t(as.data.frame(apply(ch3, 1,modify_id))))

gene<-mget(ls(pattern = "ch"))

five_UTR$ch1<-as.data.frame(t(as.data.frame(apply(five_UTR$ch1, 1,modify_id))))
five_UTR$ch2<-as.data.frame(t(as.data.frame(apply(five_UTR$ch2, 1,modify_id))))
five_UTR$ch3<-as.data.frame(t(as.data.frame(apply(five_UTR$ch3, 1,modify_id))))

percent<-data.frame(1:nrow(utr5_overlap$ch1))

ch1<-data.frame()
ch2<-data.frame()
ch3<-data.frame()

utr5_overlap_<-mget(ls(pattern ="ch"))

j=0

for(i in 1:nrow(utr5_overlap$ch1)){
    df1<-subset(gene$ch1, attributes == utr5_overlap$ch1[i,2])
    #df2<-subset(five_UTR$ch1, attributes == utr5_overlap$ch1[i,2])
    
    percent<-(as.numeric(df1$end)-as.numeric(df1$start))
    ch1<-rbind(ch1,cbind(utr5_overlap$ch1[i,1:4],percent,as.numeric(utr5_overlap$ch1[i,4])/percent))
}

for(i in 1:nrow(utr5_overlap$ch2)){
    df1<-subset(gene$ch2, attributes == utr5_overlap$ch2[i,2])
    #df2<-subset(five_UTR$ch1, attributes == utr5_overlap$ch1[i,2])
    
    percent<-(as.numeric(df1$end)-as.numeric(df1$start))
    ch2<-rbind(ch2,cbind(utr5_overlap$ch2[i,1:4],percent,as.numeric(utr5_overlap$ch2[i,4])/percent))
}

for(i in 1:nrow(utr5_overlap$ch3)){
    df1<-subset(gene$ch3, attributes == utr5_overlap$ch3[i,2])
    #df2<-subset(five_UTR$ch1, attributes == utr5_overlap$ch1[i,2])
    
    percent<-(as.numeric(df1$end)-as.numeric(df1$start))
    ch3<-rbind(ch3,cbind(utr5_overlap$ch3[i,1:4],percent,as.numeric(utr5_overlap$ch3[i,4])/percent))
}

utr5_overlap<-mget(ls(pattern = "ch"))


i<-1
while(i<=nrow(five_UTR$ch1)){
    print(i)
    j<-1
    while(j<=nrow(five_UTR$ch1)){
        if(i==j){
            next;
        }
        if(five_UTR$ch1[i,'attributes'] == five_UTR$ch1[j,'attributes']){
            five_UTR$ch1[i,'start']<-min(five_UTR$ch1[i,"start"],five_UTR$ch1[j,"start"])
            five_UTR$ch1[i,'end']<-max(five_UTR$ch1[i,"end"],five_UTR$ch1[j,"end"])
            five_UTR$ch1<-five_UTR$ch1[-j,]
        }
        j<-j+1
    }
    i<-i+1
}

list_1<-unlist(df3$`as.numeric(utr5_overlap$ch1[i, 4])/percent`)
list_2<-normalize(list_1)

for(i in 1:nrow(df3)){
    df3$`as.numeric(utr5_overlap$ch1[i, 4])/percent`[i]<-list_2[i]
}


#-----------------------------------------------------------
ch1<-subset(gene_data, type == "ncRNA" & seqid == "I")
ch2<-subset(gene_data, type == "ncRNA" & seqid == "II")
ch3<-subset(gene_data, type == "ncRNA" & seqid == "III")

ncRNA<-mget(ls(pattern = "ch"))

ch1<-find_overlapping_gene(GS$ch1, ncRNA$ch1,
                           x_start_name = "start1", x_end_name = "end2", 
                           is_strain = TRUE)

ch2<-find_overlapping_gene(GS$ch2, ncRNA$ch2,
                           x_start_name = "start1", x_end_name = "end2", 
                           is_strain = TRUE)

ch3<-find_overlapping_gene(GS$ch3, ncRNA$ch3,
                           x_start_name = "start1", x_end_name = "end2", 
                           is_strain = TRUE)

GS_ncRNA<-mget(ls(pattern = "ch"))


ch1<-find_overlapping_gene(Block$ch1, ncRNA$ch1,
                           x_start_name = "start3", x_end_name = "end4", 
                           is_strain = TRUE)

ch2<-find_overlapping_gene(Block$ch2, ncRNA$ch2,
                           x_start_name = "start3", x_end_name = "end4", 
                           is_strain = TRUE)

ch3<-find_overlapping_gene(Block$ch3, ncRNA$ch3,
                           x_start_name = "start3", x_end_name = "end4", 
                           is_strain = TRUE)

B_ncRNA<-mget(ls(pattern = "ch"))

ch1<-find_overlapping_gene(Serial$ch1, ncRNA$ch1,
                           x_start_name = "start3", x_end_name = "end4", 
                           is_strain = TRUE)

ch2<-find_overlapping_gene(Serial$ch2, ncRNA$ch2,
                           x_start_name = "start3", x_end_name = "end4", 
                           is_strain = TRUE)

ch3<-find_overlapping_gene(Serial$ch3, ncRNA$ch3,
                           x_start_name = "start3", x_end_name = "end4", 
                           is_strain = TRUE)

S_ncRNA<-mget(ls(pattern = "ch"))


#finding function of overlapping gene
library(readxl)

GO_annotation<-read_excel("DATA/gene_association.tsv.xlsx")

