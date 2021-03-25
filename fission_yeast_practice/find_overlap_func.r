find_overlapping_gene<-function(x,y, x_start_name='start', x_end_name='end',
                                y_start_name='start', y_end_name='end',
                                is_same=FALSE, is_strain=FALSE, is_type=F, b_primer=F){
    if(is_strain){
        ID<-"ID"
    }else{ ID <-"attributes" }
    if(is_strain==TRUE){
        x<-as.data.frame(t(as.data.frame(apply(x, 1, reverse_uptag_neg,x_name=x_start_name))))
        x <- type.convert(x, as.is = TRUE)
    }
    # x<-x[order(x$start1),]
    # y<-y[order(y$start),]
    len_1<-nrow(x)-1
    len_2<-nrow(y)-1
    pb <- progress_bar$new(total = len_1)
    result<-data.frame()
    
    for(i in 1:len_1){
        for(j in 1:len_2){
            if(is_same==TRUE & i==j)
                next;
            # if(x[i,x_end_name]<y[j,y_start_name])
            #     break;
            if(is_type==T){
                if((x[i,x_start_name]<=y[j,y_start_name])&(x[i,x_end_name]>=y[j,y_start_name])){
                    result<-rbind(result,cbind(x[i,ID],y[j,'attributes'],y[j,y_start_name]))
                }
            }
            else{
                overlap_percent<-NULL
                if((x[i,x_start_name]<=y[j,y_start_name])&(x[i,x_end_name]>=y[j,y_start_name])){
                    if(y[j,y_end_name]>=x[i,x_end_name]){
                        overlap_percent <- (x[i,x_end_name]-y[j,y_start_name])/(y[j,y_end_name]-y[j,y_start_name])
                    }
                    else if(y[j,y_end_name]<x[i,x_end_name]){
                        overlap_percent <- (y[j,y_end_name]-y[j,y_start_name])/(y[j,y_end_name]-y[j,y_start_name])
                    }
                    #col1<-cbind(x[i,ID],y[j,'attributes'],y[j,y_start_name],overlap_percent)
                    result<-rbind(result,cbind(x[i,ID],y[j,'attributes'],y[j,y_start_name],overlap_percent))
                }
                else if((x[i,x_start_name]<=y[j,y_end_name])&(x[i,x_end_name]>=y[j,y_end_name])){
                    if(y[j,y_start_name]<x[i,x_start_name]){
                        overlap_percent <- (y[j,y_end_name]-x[i,x_start_name])/(y[j,y_end_name]-y[j,y_start_name])
                    }
                    else if(y[j,y_start_name]>x[i,x_start_name]){
                        overlap_percent <- (y[j,y_end_name]-y[j,y_start_name])/(y[j,y_end_name]-y[j,y_start_name])
                    }
                    result<-rbind(result,cbind(x[i,ID],y[j,'attributes'],y[j,y_start_name],overlap_percent))
                }
                else if((x[i,x_start_name]>y[j,y_start_name])&(x[i,x_end_name]<y[j,y_end_name])){
                    overlap_percent <- (x[i,x_end_name]-x[i,x_start_name])/(y[j,y_end_name]-y[j,y_start_name])
                    result<-rbind(result,cbind(x[i,ID],y[j,'attributes'],y[j,y_start_name],overlap_percent))
                }
                if(!is.null(overlap_percent)){
                        if(overlap_percent<0){
                        a0<-x[i,ID]
                        a1<-x[i,x_start_name]
                        a2<-x[i,x_end_name]
                        b1<-y[j,y_start_name]
                        b2<-y[j,y_end_name]
                    }
                }
            }
        }
        pb$tick()
    }
    if(is_strain == TRUE){
        result <-as.data.frame(t(as.data.frame(apply(result, 1,modify_id))))
        result <- delete_duplicate(result)
    }else{
        result<-as.data.frame(t(as.data.frame(apply(result, 1,modify_id,is_strain=F))))
        result <- delete_duplicate(result)
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


reverse_uptag_neg<-function(x,x_name){
    #is_block<-T
    if(x_name=="start1"){
        if(x["strand1"] == "-"){
            swap(x["start1"],x["start2"])
            swap(x["end1"],x["end2"])
        }
    }
    else{    
        if(x["start3"]>x["start4"]){
            swap(x["start3"],x["start4"])
            swap(x["end3"],x["end4"])
        }
    }
    return(x)
}

modify_id<-function(x,is_strain=T){
    if(is_strain==F){
        id<-unlist(strsplit(x[1], ":"))[1]
        id_split<-unlist(strsplit(id,".",fixed = TRUE))
        result<-paste(id_split[1],id_split[2],sep = ".")
        x[1]<-gsub("ID=", "", result)
    }
    id<-unlist(strsplit(x[2], ";"))[1]
    id_split<-unlist(strsplit(id,".",fixed = TRUE))
    result<-paste(id_split[1],id_split[2],sep = ".")
    x[2]<-gsub("ID=", "", result)
    return(x)
}


delete_duplicate<-function(x){
    i<-1
    while(i<nrow(x)+1){
        # if(x[i,1]=="SPAC1039.10"){
        #     print("debug\n")
        # }
        if(x[i,1]==x[i,2]){
            x<-x[-i,]
            next
        }
        j<-1
        while(j<nrow(x)+1){
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


find_overlapping_gene_m2<-function(x,y, x_start_name='start', x_end_name='end',
                                y_start_name='start', y_end_name='end',
                                is_same=FALSE, is_strain=FALSE, is_type=F, b_primer=F){
    if(is_strain){
        ID<-"ID"
    }else{ ID <-"attributes" }
    if(is_strain==TRUE){
        x<-as.data.frame(t(as.data.frame(apply(x, 1, reverse_uptag_neg,x_name=x_start_name))))
        x <- type.convert(x, as.is = TRUE)
    }
    # x<-x[order(x$start1),]
    # y<-y[order(y$start),]
    len_1<-nrow(x)-1
    len_2<-nrow(y)-1
    pb <- progress_bar$new(total = len_1)
    result<-data.frame()
    
    for(i in 1:len_1){
        for(j in 1:len_2){
            if(is_same==TRUE & i==j)
                next;
            # if(x[i,x_end_name]<y[j,y_start_name])
            #     break;
            if(is_type==T){
                if((x[i,x_start_name]<=y[j,y_start_name])&(x[i,x_end_name]>=y[j,y_start_name])){
                    result<-rbind(result,cbind(x[i,ID],y[j,'attributes'],y[j,y_start_name]))
                }
            }
            else{
                overlap_percent<-NULL
                if((x[i,x_start_name]<=y[j,y_start_name])&(x[i,x_end_name]>=y[j,y_start_name])){
                    if(y[j,y_end_name]>=x[i,x_end_name]){
                        overlap_percent <- (x[i,x_end_name]-y[j,y_start_name])/(y[j,y_end_name]-y[j,y_start_name])
                    }
                    else if(y[j,y_end_name]<x[i,x_end_name]){
                        overlap_percent <- (y[j,y_end_name]-y[j,y_start_name])/(y[j,y_end_name]-y[j,y_start_name])
                    }
                    #col1<-cbind(x[i,ID],y[j,'attributes'],y[j,y_start_name],overlap_percent)
                    result<-rbind(result,cbind(x[i,ID],y[j,'attributes'],y[j,y_start_name],overlap_percent))
                }
                else if((x[i,x_start_name]<=y[j,y_end_name])&(x[i,x_end_name]>=y[j,y_end_name])){
                    if(y[j,y_start_name]<x[i,x_start_name]){
                        overlap_percent <- (y[j,y_end_name]-x[i,x_start_name])/(y[j,y_end_name]-y[j,y_start_name])
                    }
                    else if(y[j,y_start_name]>x[i,x_start_name]){
                        overlap_percent <- (y[j,y_end_name]-y[j,y_start_name])/(y[j,y_end_name]-y[j,y_start_name])
                    }
                    result<-rbind(result,cbind(x[i,ID],y[j,'attributes'],y[j,y_start_name],overlap_percent))
                }
                else if((x[i,x_start_name]>y[j,y_start_name])&(x[i,x_end_name]<y[j,y_end_name])){
                    overlap_percent <- (x[i,x_end_name]-x[i,x_start_name])/(y[j,y_end_name]-y[j,y_start_name])
                    result<-rbind(result,cbind(x[i,ID],y[j,'attributes'],y[j,y_start_name],overlap_percent))
                }
                if(!is.null(overlap_percent)){
                    if(overlap_percent<0){
                        a0<-x[i,ID]
                        a1<-x[i,x_start_name]
                        a2<-x[i,x_end_name]
                        b1<-y[j,y_start_name]
                        b2<-y[j,y_end_name]
                    }
                }
            }
        }
        pb$tick()
    }
    if(is_strain == TRUE){
        result <-as.data.frame(t(as.data.frame(apply(result, 1,modify_id))))
        result <- delete_duplicate(result)
    }else{
        result<-as.data.frame(t(as.data.frame(apply(result, 1,modify_id,is_strain=F))))
        result <- delete_duplicate(result)
    }
    return(result)
}


normalize <- function(x) {
    return((x-min(x))/(max(x)-min(x)))
}
