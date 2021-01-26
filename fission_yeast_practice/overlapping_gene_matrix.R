library(xlsx)
library(highcharter)
library(tidyr)
library(dplyr)

gene_id<-read.table("gene_IDs_names.tsv", sep = "\t")
gene_id_ch1 <- gene_id[1:2344,1]
gene_id_ch2 <- gene_id[1:1880,2]
gene_id_ch3 <- gene_id[1:927,3]

strain_id<-read.xlsx2('strain_id.xlsx',sheetIndex=1,header=F)
strain_id_ch1 <- strain_id[1:2232,1]
strain_id_ch2 <- strain_id[1:1808,2]
strain_id_ch3 <- strain_id[1:874,3]

#overlap_matrix<-setNames(data.frame(matrix(ncol = 6998, nrow = 4914)), gene_id)

overlap_matrix_ch1<-setNames(data.frame(matrix(ncol = 2344, nrow = 2232)), gene_id_ch1)
overlap_matrix_ch2<-setNames(data.frame(matrix(ncol = 1880, nrow = 1808)), gene_id_ch2)
overlap_matrix_ch3<-setNames(data.frame(matrix(ncol = 927, nrow = 874)), gene_id_ch3)

rownames(overlap_matrix_ch1)<-strain_id_ch1
rownames(overlap_matrix_ch2)<-strain_id_ch2
rownames(overlap_matrix_ch3)<-strain_id_ch3

overlap_matrix<-mget(ls(pattern = "overlap_matrix_"))
gene_id_list<-mget(ls(pattern = "gene_id_"))
strain_id_list<-mget(ls(pattern = "strain_id_"))

sum(!is.na(overlap_matrix$overlap_matrix))

for(i in 1:3){
    for(j in 1:nrow(GS_overlap[[i]])){
            overlap_matrix[[i]][GS_overlap[[i]][j,1], GS_overlap[[i]][j,2]]<-1
    }
}

for(i in 1:2){
    for(j in 1:nrow(S_overlap[[i]])){
        overlap_matrix[[i]][S_overlap[[i]][j,1], S_overlap[[i]][j,2]]<-1
    }
}

for(i in 1:3){
    for(j in 1:nrow(B_overlap[[i]])){
        overlap_matrix[[i]][B_overlap[[i]][j,1], B_overlap[[i]][j,2]]<-1
    }
}

for(i in 1:3){
    pb <- progress_bar$new(total = length(strain_id_list[[i]]))
    for(j in 1:length(strain_id_list[[i]])){
        overlap_matrix[[i]][strain_id_list[[i]][j], strain_id_list[[i]][j]]<-1
        pb$tick()
    }
}
sum(!is.na(overlap_matrix$overlap_matrix_ch1))
overlap_matrix$overlap_matrix_ch1[is.na(overlap_matrix$overlap_matrix_ch1)] <- 0
overlap_matrix$overlap_matrix_ch2[is.na(overlap_matrix$overlap_matrix_ch2)] <- 0
overlap_matrix$overlap_matrix_ch3[is.na(overlap_matrix$overlap_matrix_ch3)] <- 0
write_xlsx(overlap_matrix$overlap_matrix_ch1,"overlap_matrix_ch1.xlsx")
heatmap(overlap_matrix$overlap_matrix_ch1)

hchart(overlap_matrix$overlap_matrix_ch1, "heatmap") %>% 
    hc_colorAxis(stops = color_stops(2, c("red","green")))

# gene_tmp<-subset(gene_data,seqid == "I" & type == "mRNA")
# gene_temp<-gene_tmp$attributes
# gene_temp
# 
# extract_id<-function(x){
#     id<-unlist(strsplit(x[1], ";"))[1]
#     id_split<-unlist(strsplit(id,".",fixed = TRUE))
#     result<-paste(id_split[1],id_split[2],sep = ".")
#     x[1]<-gsub("ID=", "", result)
#     return(x)
# }
# 
# gene_temp<-as.data.frame(apply(as.data.frame(gene_temp), 1,extract_id))
# gene_temp <- type.convert(gene_temp, as.is = TRUE)
