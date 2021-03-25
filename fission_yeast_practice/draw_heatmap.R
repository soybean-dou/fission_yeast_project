
library(pheatmap)

ch1<-rbind(GS_CDS$ch1,B_CDS$ch1,S_CDS$ch1)

ch2<-rbind(GS_CDS$ch2,B_CDS$ch2,S_CDS$ch2)

ch3<-rbind(GS_CDS$ch3,B_CDS$ch3,S_CDS$ch3)

CDS_overlap<-mget(ls(pattern = "ch"))

for(k in 1:3){
    CDS_overlap[[k]]<-CDS_overlap[[k]][order(CDS_overlap[[k]][1]),]
    
    strain_name<-CDS_overlap[[k]]$V1
    cds_name<-CDS_overlap[[k]]$V2
    
    i=1;
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
            j<-j+1
        }
        i<-i+1
    }
    
    i=1;
    while(i<as.numeric(length(cds_name))){
        j=1;
        while(j<length(cds_name)){
            if(i==j){
                j<-j+1
                next
            }
            if(cds_name[i]==cds_name[j]){
                cds_name<-cds_name[-j]
            }
            j<-j+1
        }
        i<-i+1
    }
    
    cds_matrix<-matrix(ncol=length(cds_name),nrow=length(strain_name))
    
    rownames(cds_matrix)<-strain_name
    colnames(cds_matrix)<-cds_name
    
    for(h in 1:nrow(CDS_overlap[[k]])){
        cds_matrix[CDS_overlap[[k]][h,1], CDS_overlap[[k]][h,2]]<-as.numeric(CDS_overlap[[k]][h,6])
    }
    
    png(filename=paste0("strain_cds_ch",k,".png"),width=1000,height=1000,unit="px",bg="transparent")
    p <- pheatmap(cds_matrix, color = colorRampPalette(c("yellow", "red"))(100),
                  cluster_row = FALSE, cluster_cols = F, main = paste("Strain-CDS overlap percent in chromosome",as.character(k)), na_col = "white")
    dev.off()
    
}

ch1<-rbind(GS_5utr$ch1,B_5utr$ch1,S_5utr$ch1)

ch2<-rbind(GS_5utr$ch2,B_5utr$ch2,S_5utr$ch2)

ch3<-rbind(GS_5utr$ch3,B_5utr$ch3,S_5utr$ch3)

utr5_overlap<-mget(ls(pattern = "ch"))

for(k in 1:1){
    utr5_overlap[[k]]<-utr5_overlap[[k]][order(utr5_overlap[[k]][1]),]
    
    strain_name<-utr5_overlap[[k]]$V1
    utr5_name<-utr5_overlap[[k]]$V2
    
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
            j<-j+1
        }
        i<-i+1
    }
    
    i=1;
    while(i<as.numeric(length(utr5_name))){
        j=1;
        while(j<length(utr5_name)){
            if(i==j){
                j<-j+1
                next
            }
            if(utr5_name[i]==utr5_name[j]){
                utr5_name<-utr5_name[-j]
            }
            j<-j+1
        }
        i<-i+1
    }
    
    utr5_matrix<-matrix(ncol=length(utr5_name),nrow=length(strain_name))
    
    utr5_name
    strain_name
    
    rownames(utr5_matrix)<-strain_name
    colnames(utr5_matrix)<-utr5_name
    
    for(h in 1:nrow(utr5_overlap[[k]])){
        utr5_matrix[utr5_overlap[[k]][h,1], utr5_overlap[[k]][h,2]]<-as.numeric(df3[h,6])
    }
    
    png(filename=paste0("strain_5utr_ch",k,".png"),width=1200,height=800,unit="px",bg="transparent")
    p <- pheatmap(utr5_matrix, color = colorRampPalette(c("yellow", "red"))(100),
                  cluster_row = FALSE, cluster_cols = F, main = paste("Strain-5UTR overlap percent in chromosome",as.character(k)), na_col = "white")
    dev.off()

}


ch1<-rbind(GS_3utr$ch1,B_3utr$ch1,S_3utr$ch1)

ch2<-rbind(GS_3utr$ch2,B_3utr$ch2,S_3utr$ch2)

ch3<-rbind(GS_3utr$ch3,B_3utr$ch3,S_3utr$ch3)

utr3_overlap<-mget(ls(pattern = "ch"))

for(k in 1:3){
    utr3_overlap[[k]]<-utr3_overlap[[k]][order(utr3_overlap[[k]][1]),]
    
    strain_name<-utr3_overlap[[k]]$V1
    utr3_name<-utr3_overlap[[k]]$V2
    
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
            j<-j+1
        }
        i<-i+1
    }
    
    i=1;
    while(i<as.numeric(length(utr3_name))){
        j=1;
        while(j<length(utr3_name)){
            if(i==j){
                j<-j+1
                next
            }
            if(utr3_name[i]==utr3_name[j]){
                utr3_name<-utr3_name[-j]
            }
            j<-j+1
        }
        i<-i+1
    }
    
    utr3_matrix<-matrix(ncol=length(utr3_name),nrow=length(strain_name))
    
    strain_name
    
    rownames(utr3_matrix)<-strain_name
    colnames(utr3_matrix)<-utr3_name
    
    for(h in 1:nrow(utr3_overlap[[k]])){
        utr3_matrix[utr3_overlap[[k]][h,1], utr3_overlap[[k]][h,2]]<-as.numeric(utr3_overlap[[k]][h,4])
    }
    
    png(filename=paste0("strain_3utr_ch",k,".png"),width=1000,height=1000,unit="px",bg="transparent")
    p <- pheatmap(utr3_matrix, color = colorRampPalette(c("yellow", "red"))(100),
                  cluster_row = FALSE, cluster_cols = F, main = paste("Strain-3UTR overlap percent in chromosome",as.character(k)), na_col = "white")
    dev.off()
}

ch1<-rbind(GS_ncRNA$ch1,B_ncRNA$ch1,S_ncRNA$ch1)

ch2<-rbind(GS_ncRNA$ch2,B_ncRNA$ch2,S_ncRNA$ch2)

ch3<-rbind(GS_ncRNA$ch3,B_ncRNA$ch3,S_ncRNA$ch3)

ncRNA_overlap<-mget(ls(pattern = "ch"))

for(k in 1:3){
    ncRNA_overlap[[k]]<-ncRNA_overlap[[k]][order(ncRNA_overlap[[k]][1]),]
    
    strain_name<-ncRNA_overlap[[k]]$V1
    ncRNA_name<-ncRNA_overlap[[k]]$V2
    
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
            j<-j+1
        }
        i<-i+1
    }
    
    i=1;
    while(i<as.numeric(length(ncRNA_name))){
        j=1;
        while(j<length(ncRNA_name)){
            if(i==j){
                j<-j+1
                next
            }
            if(ncRNA_name[i]==ncRNA_name[j]){
                ncRNA_name<-ncRNA_name[-j]
            }
            j<-j+1
        }
        i<-i+1
    }
    
    ncRNA_matrix<-matrix(ncol=length(ncRNA_name),nrow=length(strain_name))
    
    strain_name
    
    rownames(ncRNA_matrix)<-strain_name
    colnames(ncRNA_matrix)<-ncRNA_name
    
    for(h in 1:nrow(ncRNA_overlap[[k]])){
        ncRNA_matrix[ncRNA_overlap[[k]][h,1], ncRNA_overlap[[k]][h,2]]<-as.numeric(ncRNA_overlap[[k]][h,4])
    }
    
    png(filename=paste0("strain_ncRNA_ch",k,".png"),width=1200,height=800,unit="px",bg="transparent")
    p <- pheatmap(ncRNA_matrix, color = colorRampPalette(c("yellow", "red"))(100),
                  cluster_row = FALSE, cluster_cols = F, main = paste("Strain-ncRNA overlap percent in chromosome",as.character(k)), 
                  na_col = "gray98", fontsize_row = 5, fontsize_col = 5)
    dev.off()
}

dev.new()
