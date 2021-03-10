
library(pheatmap)

ch1<-rbind(GS_CDS$ch1,B_CDS$ch1,S_CDS$ch1)

ch2<-rbind(GS_CDS$ch2,B_CDS$ch2,S_CDS$ch2)

ch3<-rbind(GS_CDS$ch3,B_CDS$ch3,S_CDS$ch3)

CDS_overlap<-mget(ls(pattern = "ch"))

for(k in 1:3){
    CDS_overlap[[k]]<-CDS_overlap[[k]][order(CDS_overlap[[k]][1]),]
    
    strain_name<-CDS_overlap[[k]]$V1
    
    i=1;
    j=1;
    as.numeric(length(strain_name))
    while(i<as.numeric(length(strain_name))){
        j=1;
        while(j<length(strain_name)){
            if(i==j){
                j<-j+1
                next
            }
            if(strain_name[i]==strain_name[j]){
                strain_name<-strain_name[-j]
            }
            print(j)
            j<-j+1
        }
        print(i)
        i<-i+1
    }
    
    cds_matrix<-matrix(ncol=nrow(CDS_overlap[[k]]),nrow=length(strain_name))
    
    rownames(cds_matrix)<-strain_name
    colnames(cds_matrix)<-CDS_overlap[[k]]$V2
    
    for(h in 1:nrow(CDS_overlap[[k]])){
        cds_matrix[CDS_overlap[[k]][h,1], CDS_overlap[[k]][h,2]]<-as.numeric(CDS_overlap[[k]][h,4])
    }
    
    p <- pheatmap(cds_matrix, color = colorRampPalette(c("yellow", "red"))(100),
                  cluster_row = FALSE, cluster_cols = F, main = paste("Strain-CDS overlap percent in chromosome",as.character(k)), na_col = "white")
}









cds_matrix
cds_matrix[is.na(cds_matrix)]<-as.numeric(0)

heatmap(cds_matrix)



