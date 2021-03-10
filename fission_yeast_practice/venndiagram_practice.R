library(VennDiagram)
library("ggVennDiagram")


cds_sample<-read.csv("DATA/cds_data.CSV")
utr_sample<-read.csv("DATA/utr_data.CSV")

overlap_strain_cds<-rbind(GS_overlap$overlap_strand_gs1,B_overlap$overlap_strand_B1,S_overlap$overlap_strand_S1,
                          GS_overlap$overlap_strand_gs2,B_overlap$overlap_strand_B2,S_overlap$overlap_strand_S3,
                          GS_overlap$overlap_strand_gs3,B_overlap$overlap_strand_B3,S_overlap$overlap_strand_S3)

colnames(overlap_strain_cds)<- c("Strain","CDS")

write_xlsx(overlap_strain_cds,"overlap_strain_cds.xlsx")

x<-list(
    current = overlap_strain_cds[,1],
    before = cds_sample[,1]
)
ggVennDiagram(x)

overlap_strain_5utr<-rbind(overlap_strand_UTR_gs1,overlap_strand_UTR_B1,overlap_strand_UTR_S1,
                           overlap_strand_UTR_gs2,overlap_strand_UTR_B2,overlap_strand_UTR_S2,
                           overlap_strand_UTR_gs3,overlap_strand_UTR_B3,overlap_strand_UTR_S3)

x<-list(
    current = overlap_strain_5utr[,1],
    before = utr_sample[,1]
)
ggVennDiagram(x, label_alpha = 0)

overlap_strand_3utr<-rbind(overlap_strand_UTR_B1,overlap_strand_UTR_gs1,overlap_strand_UTR_S1,
                           overlap_strand_UTR_B2,overlap_strand_UTR_gs2,overlap_strand_UTR_S2,
                           overlap_strand_UTR_B3,overlap_strand_UTR_gs3,overlap_strand_UTR_S3)

x<-list(
    current = overlap_strand_3utr[,1],
    before = utr_sample[,1]
)
ggVennDiagram(x, label_alpha = 0)
