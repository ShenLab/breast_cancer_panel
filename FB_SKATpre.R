FB_SKATpre <- function(){
    source("SKAT_ana.R")
    source("Faminfo.R")
    pedf <- "data/ALL_pedigree.csv"
    peds <- read.csv(pedf)
    
    pheno <- pheno_all()
    
    idi <- as.vector(pheno[,2])
    fullkins <- getKinshipMatrix(1)
    allped <- colnames(fullkins)
    idi <- intersect(idi,allped)
    id <- pheno[match(idi,pheno[,2]),3]
    colnames(fullkins)[match(idi,colnames(fullkins))] <- id
    rownames(fullkins)[match(idi,rownames(fullkins))] <- id
    
    Z1 <- genotype_wts_FBSKAT(id,FALSE)
    
    phe <- as.vector(pheno[match(id,pheno[,3]),"BreastCancer"])
    phe[phe=="Yes"] <- 1
    phe[phe=="No"] <- 0
    phe <- as.numeric(phe)
    
    sex <- as.vector(pheno[match(id,pheno[,3]),"Sex"])
    sex[sex=="Male"] <- 1
    sex[sex=="Female"] <- 2
    sex <- as.numeric(sex)
    
    peds <- peds[match(idi,peds[,2]),]
    peds[,2] <- id
    ip <- intersect(peds[,3],pheno[,2])
    ip0 <- peds[,3]
    peds[,3] <- 9
    im <- intersect(peds[,4],pheno[,2])
    im0 <- peds[,4]
    peds[,4] <- 9
    
    peds[match(ip,ip0),3] <- pheno[match(ip,pheno[,2]),3]
    peds[match(im,im0),4] <- pheno[match(im,pheno[,2]),3]
    
    peds <- cbind(peds,sex,phe,Z1)
    write.table(peds,file="family/pedigree_FBSKAT.ped",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    
    
}


genotype_wts_FBSKAT <- function(id,sig,dirstr="famSKATresult/"){
    
    genos <- genotype_wts(id,sig)
    #load(paste(dirstr,"genos_",sig,sep=""))
    alllist <- genos$alllist
    Z <- genos$Z
    fres <- genos$fres
    vars <- paste(alllist[,1],alllist[,2],alllist[,4],alllist[,5],sep="_")
    var_g <- cbind(alllist[,"Gene"],vars)
    
    Z1 <- c()
    Z1n <- c()
    genes <- unique(alllist[,"Gene"])
    
    k <- 0
    gst <- rep(0,length(genes))
    for(i in 1:length(genes)){
        onevar <- unique(vars[var_g[,1]==genes[i]])
        Z1 <- cbind(Z1,Z[,onevar])
        Z1n <- c(Z1n,onevar)
        st <- k+1 
        ed <- k+length(onevar)
        #st <- st+6
        #ed <- ed+6
        gst[i] <- paste(genes[i],st,ed,sep="\t")
        k <- ed
    }
    colnames(Z1) <- Z1n
    write.table(gst,file="family/gene.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
    f <- fres[colnames(Z1)]
    write.table(f,file="family/weights.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
    varpass <- alllist[match(Z1n,vars),]
    write.table(varpass,file="family/var_pass.txt",row.names=FALSE,quote=FALSE)
    
    Z1
}
