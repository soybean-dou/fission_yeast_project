library(xlsx)

gene_id<-read.table("gene_IDs_names.tsv", sep = "\t")
gene_id<-gene_id$V1

strain_id<-read.xlsx2('strain_id.xlsx',sheetIndex=1,header=TRUE)
strain_id<-strain_id$Systemic.ID

overlap_matrix<-setNames(data.frame(matrix(ncol = 6998, nrow = 4914)), gene_id)
rownames(overlap_matrix)<-strain_id
