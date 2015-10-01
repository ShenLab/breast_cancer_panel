FB_SKATpre <- function(){
    source("SKAT_ana.R")
    source("Faminfo.R")
    pedf <- "family/Pedigree_sec.txt"
    peds <- read.delim(pedf,sep="\t")
    
    pheno <- pheno_all()
    id <- as.vector(pheno[,3])
    fullkins <- getKinshipMatrix(1)
    allped <- colnames(fullkins)
    id <- intersect(id,allped)
    Z1 <- genotype_wts_FBSKAT(id,FALSE)
    
    phe <- as.vector(pheno[match(id,pheno[,3]),"BreastCancer"])
    phe[phe=="Yes"] <- 1
    phe[phe=="No"] <- 0
    phe <- as.numeric(phe)
    
    sex <- as.vector(pheno[match(id,pheno[,3]),"Sex"])
    sex[sex=="Male"] <- 1
    sex[sex=="Female"] <- 2
    sex <- as.numeric(sex)
    
    peds <- peds[match(id,peds[,2]),]
    peds[peds[,3]==0,3] <- 9
    peds[peds[,4]==0,4] <- 9
    
    peds <- cbind(peds[,1:4],sex,phe,Z1)
    write.table(peds,file="family/pedigree_FBSKAT.ped",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    
    
}


genotype_wts_FBSKAT <- function(id,sig,dirstr="famSKATresult/"){
    
    #genos <- genotype_wts(id,sig)
    load(paste(dirstr,"genos_",sig,sep=""))
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
        gst[i] <- paste(genes[i],st,ed,sep="\t")
        k <- ed
    }
    colnames(Z1) <- Z1n
    write.table(gst,file="family/gene.ped",row.names=FALSE,col.names=FALSE,quote=FALSE)
    f <- fres[colnames(Z1)]
    write.table(f,file="family/weights.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
    
    Z1
}
