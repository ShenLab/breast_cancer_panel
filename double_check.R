Mendelian_check <- function(){
    ## pedigree information predicted by Primus with second degree information
    pedf <- "family/Pedigree_sec.txt"
    ped <- read.delim(pedf)
    subs <- ped[,3]==0 & ped[,4]==0
    ped <- ped[!subs,]
    ped <- ped[,c(2,3,4)]
    
    ## pedigree information predicted by Primus with first degree information
    pedf <- "family/Pedigree_fir.txt"
    ped_1 <- read.delim(pedf)
    subs <- ped_1[,3]==0 & ped_1[,4]==0
    ped_1 <- ped_1[!subs,]
    ped_1 <- ped_1[,c(2,3,4)]
    
    tmp <- intersect(ped_1[,1],ped[,1])
    ped <- ped[!(ped[,1] %in% tmp),]
    
    preped <- rbind(ped,ped_1)
    
    ## manual checked ped
    pedf1 <- "family/Famcheckedids.csv"
    ped1 <- read.csv(pedf1)
    subs <- is.na(ped1[,5]) & is.na(ped1[,6])
    ped1 <- ped1[!subs,]
    ped1 <- ped1[,c("Subject_ID","PID_ID","MID_ID")]
    
    pedf2 <- "family/PriorFam_Qiang.csv"
    pedf3 <- "family/PriorFam2_Qiang.csv"
    ped2 <- read.csv(pedf2)
    ped3 <- read.csv(pedf3)
    ped2 <- rbind(ped2,ped3)
    subs <- is.na(ped2[,"PID.1"]) & is.na(ped2[,"MID.1"])
    ped2 <- ped2[!subs,]
    ped2 <- ped2[,c("Subject_ID","PID.1","MID.1")]
    colnames(ped2) <- c("Subject_ID","PID_ID","MID_ID")
    manped <- rbind(ped1,ped2)
    
    tmp <- intersect(manped[,1],preped[,1])
    preped <- preped[!(preped[,1] %in% tmp),]
    colnames(preped) <- c("Subject_ID","PID_ID","MID_ID")
    
    fullped <- rbind(manped,preped)
    fullped[is.na(fullped)] <- 0
    fullped1 <- fullped[fullped[,3]!=0 & fullped[,2]!=0,]
    
    list(fullped=fullped,fullped1=fullped1)
}

lookupcheck <- function(){
    re <- Mendelian_check()
    fullped <- re$fullped
    fullped1 <- re$fullped1
    ## all checked mutations
    mutaf <- "lookup.tsv"
    mutas <- read.delim(mutaf)
    chvars <- paste(mutas[,"Chromosome"],mutas[,"Position"],mutas[,"REF"],mutas[,"ALT"],sep="_")
    a <- cbind(chvars,mutas[,"Subject_ID"],"","","","")
    colnames(a) <- c("variant","Subject_ID","oneParent","Menerror","twoParent","Menerror2")
    a[a[,2] %in% fullped[,1],3] <- TRUE
    a[a[,2] %in% fullped1[,1],5] <- TRUE
    
    ## all mutations in families
    load("alllist_9_30")
    allV <- paste(alllist[,1],alllist[,2],alllist[,4],alllist[,5],sep="_")
    varT <- allV
    
    subs <- which(a[,3]==TRUE)
    for(i in subs){
        pars <- fullped[fullped[,1]==a[i,2],2:3]
        onep <- which(alllist[,"Subject_ID"] %in% pars)
        a[i,4] <- a[i,1] %in% varT[onep]
    }
    
    subs <- which(a[,5]==TRUE)
    for(i in subs){
        pars <- fullped1[fullped1[,1]==a[i,2],2:3]
        onep <- which(alllist[,"Subject_ID"] %in% pars)
        a[i,6] <- a[i,1] %in% varT[onep]
    }
    source("~/.Rprofile")
    qwt(a,file="MendelianError_lookup.txt",flag=2)
    
}

onecheck <- function(){
    re <- Mendelian_check()
    fullped <- re$fullped
    fullped1 <- re$fullped1
    ## all checked mutations
    mutaf <- "filtered_Qiang/mutationsinfo.txt"
    mutas <- read.delim(mutaf)
    chvars <- paste(mutas[,"Chromosome"],mutas[,"Position"],mutas[,"REF"],mutas[,"ALT"],sep="_")
    a <- cbind(chvars,mutas[,"Subject_ID"],"","","","")
    colnames(a) <- c("variant","Subject_ID","oneParent","Menerror","twoParent","Menerror2")
    a[a[,2] %in% fullped[,1],3] <- TRUE
    a[a[,2] %in% fullped1[,1],5] <- TRUE
    
    ## all mutations in families
    load("alllist_9_10")
    allV <- paste(alllist[,1],alllist[,2],alllist[,4],alllist[,5],sep="_")
    varT <- allV
    
    subs <- which(a[,3]==TRUE)
    for(i in subs){
        pars <- fullped[fullped[,1]==a[i,2],2:3]
        onep <- which(alllist[,"Subject_ID"] %in% pars)
        a[i,4] <- a[i,1] %in% varT[onep]
    }
    
    subs <- which(a[,5]==TRUE)
    for(i in subs){
        pars <- fullped1[fullped1[,1]==a[i,2],2:3]
        onep <- which(alllist[,"Subject_ID"] %in% pars)
        a[i,6] <- a[i,1] %in% varT[onep]
    }
    source("~/.Rprofile")
    qwt(a,file="MendelianError.txt",flag=2)

}

DP4_check <- function(){

    dp4f <- "filtered_Qiang/gene_mutatons_dp4_9_24.txt"
    varT <- read.delim(dp4f,sep="\t")
    dpfil <- rep(TRUE,dim(varT)[1])
    
    for(i in 1:dim(varT)[1]){
        onedp <- varT[i,"alldp4"]
        tmp <- unlist(strsplit(onedp,";"))
        onevar <- paste(varT[i,"REF"],varT[i,"ALT"],sep="_")
        nj <- length(tmp)
        onech <- TRUE
        for(j in 1:nj){
            tmp1 <- unlist(strsplit(tmp[j],":"))
            nk <- length(tmp1)
            if(nk >= 5){onedp4 <- tmp1[5];}else{onedp4 <- tmp1[3];}
            onetmp <- unlist(strsplit(onedp4,","))
            onech=onech & (onetmp[3]==0 & onetmp[4]==0)
        }
        if(is.na(onech)) print(i)
        dpfil[i] <- onech
    }
    dpfil[is.na(dpfil)] <- FALSE
    varT <- cbind(varT,dpfil)
    
    source("~/.Rprofile")
    qwt(varT,file="DP4error.txt",flag=2)
    
}

vcf_check_WRN_CBL0 <- function(){
    varT <- read.delim("gene_mutations_dp4.txt",sep="\t")
    
    #subs <- which(varT[,3]==30945376)
    subs <- which(varT[,3]==119149355)
    a <- matrix(,6,5)
    tmp <- varT[subs,c("Subject_IDs","GT","AD","GT.AD.DP.GQ.PL","alldp4")]
    for(i in 1:5){
        if(i==2 | i==3){
            a0 <- unlist(strsplit(tmp[1,i],";"))
            a1 <- sapply(1:length(a0),function(k) unlist(strsplit(a0[k],":"))[2])
            a[,i] <- a1
        }else if(i==4){
            a0 <- unlist(strsplit(tmp[1,i],";"))
            a1 <- sapply(1:length(a0),function(k) substr(a0[k],8,nchar(a0[k])))
            a[,i] <- a1
        }else if(i==5){
            a0 <- unlist(strsplit(tmp[1,i],";"))
            a1 <- sapply(1:length(a0),function(k)  paste(unlist(strsplit(a0[k],":"))[4],unlist(strsplit(a0[k],":"))[5],sep=":" )     )
            a[,i] <- a1        
        }else{
            a[,i] <- unlist(strsplit(tmp[1,i],","))
        }
    }
    
    write.csv(a,file="WRN_CBL.csv",row.names=FALSE)
}

vcf_check_WRN_CBL <- function(){
    onef <- "data/CBLvcf.txt"
    #onef <- "data/WRNvcf.txt"
    a <- read.delim(onef,sep="\t",header=FALSE)
    ad <- sapply(10:1366,function(i) unlist(strsplit(a[1,i],":"))[2])
    ge <- sapply(10:1366,function(i) unlist(strsplit(a[1,i],":"))[1])
    alt  <- sapply(1:length(ad), function(i) unlist(strsplit(ad[i],",")) )
    alt <- t(alt)
    mode(alt) <- "numeric"
    subs <- sapply(10:1366,function(i) grepl("0/0:",a[1,i]) | grepl("\\./\\.:",a[1,i]) )
    
    alt0 <- alt[subs,]
    
    table(alt0[,2])
    table(alt0[,3])
    
}

famSKAT_check <- function(){

    alllist <- genos$alllist
    gener <- read.delim("famSKATresult/Gene_FALSE_All_ALL_1.txt")
    gener <- gener[order(gener[,6]),]
    qwt(gener,file="famSKAT_gene.txt",flag=2)
    subg <- c()
    for (i in 1:dim(gener)[1]){
        onel <- which(alllist[,"Gene"]==gener[i,1])
        subg <- c(subg,onel)
        }
    genelist <- alllist[subg,]
    qwt(genelist,file="famSKAT_genelist.txt",flag=2)
    
    varr <- read.delim("famSKATresult/Variant_FALSE_All_1.txt")
    varr <- varr[order(varr[,6]),]
    qwt(varr,file="famSKAT_vars.txt",flag=2)
    
    varV <- paste(alllist[,1],alllist[,2],alllist[,4],alllist[,5],sep="_")
    tmp <- do.call(c,sapply(1:dim(varr)[1],function(i) which(varV==varr[i,2]) ))
    varlist <- alllist[tmp,]
    qwt(varlist,file="famSKAT_varlist.txt",flag=2)
    
}
