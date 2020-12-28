delete_duplicate_factor<-function(x){
    i=1
    while(i<=nrow(x)){
        ID_1<-unlist(strsplit(x[i,1], ":"))[1]
        ID_3<-unlist(strsplit(x[i,2], ":"))[1]
        ID_1_comp<-unlist(strsplit(ID_1,".",fixed = TRUE))
        ID_3_comp<-unlist(strsplit(ID_3,".",fixed = TRUE))
        if(ID_1_comp[1]==ID_3_comp[1] & ID_1_comp[2]==ID_3_comp[2]){
            x<-x[-i,]
            next
        }
        j=1
        while(j<=nrow(x)){
            ID_2<-unlist(strsplit(x[j,1], ":"))[1]
            if(j==i){ j<-j+1; next; }
            if(ID_1 == ID_2){
                x<-x[-j,]
            }
            j<-j+1
        }
        i<-i+1
    }
    return(x)
}

find_overlapping_gene<-function(x,y, x_start_name='start', x_end_name='end',
                                y_start_name='start', y_end_name='end',
                                is_same=FALSE, is_strain=FALSE){
    if(is_strain==TRUE){
        x<-as.data.frame(t(as.data.frame(apply(x, 1,reverse_uptag_neg))))
        x <- type.convert(x, as.is = TRUE)
    }
    x<-x[order(x$start1),]
    y<-y[order(y$start),]
    len_1<-nrow(x)-1
    len_2<-nrow(y)-1
    pb <- progress_bar$new(total = len_1)
    result<-data.frame()
    for(i in 1:len_1){
        for(j in 1:len_2){
            if(is_same==TRUE & i==j)
                next;
            if(x[i,x_end_name]<y[j,y_start_name])
                break;
            if((x[i,x_start_name]<=y[j,y_start_name])&(x[i,x_end_name]>=y[j,y_start_name])){
                result<-rbind(result,cbind(x[i,'ID'],y[j,'attributes'],y[j,y_start_name]))
            }
            else if((x[i,x_start_name]<=y[j,y_end_name])&(x[i,x_end_name]>=y[j,y_end_name])){
                result<-rbind(result,cbind(x[i,'ID'],y[j,'attributes'],y[j,y_start_name]))
            }
            else if((x[i,x_start_name]>y[j,y_start_name])&(x[i,x_end_name]<y[j,y_end_name])){
                result<-rbind(result,cbind(x[i,'ID'],y[j,'attributes'],y[j,y_start_name]))
            }
        }
        pb$tick()
    }
    if(is_strain == TRUE){
        result<-as.data.frame(t(as.data.frame(apply(result, 1,modify_id))))
        result <- type.convert(result, as.is = TRUE)
        #result <- delete_duplicate(result)
    }
    return(result)
}

reverse_negative_strand<-function(x, up_start, up_end, down_start, down_end,
                                  up_strand, down_strand){
    if(x[up_strand] == "-"){
        swap(x[up_start],x[up_end])
    }
    else if(x[down_strand] == "-"){
        swap(x[down_start],x[down_end])
    }
    return(x)
}


reverse_uptag_neg<-function(x){
    if(x["strand1"] == "-"){
        swap(x["start1"],x["start2"])
        swap(x["end1"],x["end2"])
    }
    return(x)
}

modify_id<-function(x){
    id<-unlist(strsplit(x[2], ":"))[1]
    id_split<-unlist(strsplit(id,".",fixed = TRUE))
    result<-paste(id_split[1],id_split[2],sep = ".")
    x[2]<-gsub("ID=", "", result)
    return(x)
}


delete_duplicate<-function(x){
    i<-1
    while(i<nrow(x)){
        if(x[i,1]==x[i,2]){
            x<-x[-i,]
            next
        }
        j<-1
        while(j<nrow(x)){
            test1<-x[i,]
            test2<-x[j,]
            if(x[i,1]==x[j,1] & i!=j & x[i,2]==x[j,2]){
                x<-x[-j,]
                next
            }
            j<-j+1
        }
        i<-i+1
    }
    return(x)
}
