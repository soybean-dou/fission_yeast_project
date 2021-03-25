i<-1
while(i<=nrow(five_UTR$ch1)){
    print(i)
    j<-1
    while(j<=nrow(five_UTR$ch1)){
        if(i==j){
            j<-j+1
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

