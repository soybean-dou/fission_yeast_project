#finding overlaping gene(gene_data-gene_data)



library(ape)
library(stringr)
library(dplyr)
library(writexl)
library(progress)
library(seqinr)
library(RColorBrewer)
source('find_overlap_func.r')

gene_data<-read.gff("Schizosaccharomyces_pombe_all_chromosomes.gff3")

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

for(i in 1:3){
    file_name <- paste0("CDS_overlap_", toString(i),".tsv")
    write_xlsx(CDS_overlap[i],file_name)
    
    file_name <- paste0("CDS_", toString(i),".tsv")
    write_xlsx(CDS_site[i],file_name)
}
