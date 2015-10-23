parallelfamSKAT <- function(){
    source("FB_SKATpre.R")
    library(parallel)
    
    mclapply(1:20,function(kk) run_FBSKAT(kk,FALSE),mc.cores = 20)
    mclapply(1:20,function(kk) run_FBSKAT(kk,TRUE),mc.cores = 20)
    
    ## bash parallel run
    ##parallel './runFBSKAT.sh $i' ::: {1..40} 
}

run_FBSKAT <- function(kk,sig){
    
    dirstr <- "FBSKATresult/"
    vartype <- c("LGD","D-mis","indels","LGD+D-mis","ALL")
    poptype <- c("Jewish","Hispanic","JH","All")
    
    source("SKAT_ana.R")
    source("Faminfo.R")
    pedf <- "data/ALL_pedigree.csv"
    peds <- read.csv(pedf)
    pheno <- pheno_all()
    pheno[,3] <- gsub("222357, 222966","222357",pheno[,3])
    load("alllist_10_20")
    listsubs <- unique(alllist[,"Subject_ID"])
    pheno <- pheno[pheno[,3] %in% listsubs,]
    
    idi <- as.vector(pheno[,2])
    #fullkins <- getKinshipMatrix(1)
    #allped <- colnames(fullkins)
    allped <- union(peds[,2],union(peds[,3],peds[,4]))
    idi <- intersect(idi,allped)
    id <- pheno[match(idi,pheno[,2]),3]
    #colnames(fullkins)[match(idi,colnames(fullkins))] <- id
    #rownames(fullkins)[match(idi,rownames(fullkins))] <- id
    
    phe <- as.vector(pheno[match(id,pheno[,3]),"BreastCancer"])
    phe[phe=="Yes"] <- 1
    phe[phe=="No"] <- 0
    phe <- as.numeric(phe)
    
    sex <- as.vector(pheno[match(id,pheno[,3]),"Sex"])
    sex[sex=="Male"] <- 1
    sex[sex=="Female"] <- 2
    sex <- as.numeric(sex)
    
    peds0 <- matrix(9,length(id),4) ## 9 for miss 
    peds0[,1] <- pheno[match(id,pheno[,3]),1]
    peds0[,2] <- id
    for(i in 1:length(id)){
        if(any(peds[,2]==idi[i])){
        fid <- peds[which(peds[,2]==idi[i]),3]
        mid <- peds[which(peds[,2]==idi[i]),4]
        
        if(fid %in% pheno[,2]){
            peds0[i,3] <- pheno[pheno[,2]==fid,3]
        }
        if(mid %in% pheno[,2]){
            peds0[i,4] <- pheno[pheno[,2]==mid,3]
        }
        
        }
    } 
    peds <- peds0
    
    ### run with one flag, sig and flag
    fig <- floor((kk-1)/4) + 1
    pop <- ifelse(kk %% 4==0,4,kk %% 4)
    genos <- genotype_wts(id,sig)
    
    oneresult <- subFBSKAT(genos,fig,pop,peds,sex,phe,id)
    pedsub <- oneresult$peds
    sexsub <- oneresult$sex
    phesub <- oneresult$phe
    Z1sub <- oneresult$Z1
    ## try to use genotype data
    #Z1sub[Z1sub > 0] <- 2
    #Z1sub[Z1sub == 0] <- 1
    
    gstsub <- oneresult$gst
    fres <- oneresult$fres
    wtssub <- fres #dbeta(fres,1,25)
    varsub <- oneresult$onelist
    
    peds <- cbind(pedsub,sexsub,phesub,Z1sub)
    ## order a trio need to be in order child, parent1, parent2
    peds <- changepeds(peds)
    
    write.table(peds,file=paste(dirstr,"pedigree_FBSKAT_",vartype[fig],"_",poptype[pop],"_",sig,".ped",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    write.table(gstsub,file=paste(dirstr,"gene_FBSKAT_",vartype[fig],"_",poptype[pop],"_",sig,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
    if(fig==5){
        tmp <- cbind(colnames(Z1sub),1:length(wtssub),1:length(wtssub))
        write.table(tmp,file=paste(dirstr,"variant_FBSKAT_",vartype[fig],"_",poptype[pop],"_",sig,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
    write.table(wtssub,file=paste(dirstr,"wts_FBSKAT_",vartype[fig],"_",poptype[pop],"_",sig,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(rep(1,length(wtssub)),file=paste(dirstr,"wts1_FBSKAT_",vartype[fig],"_",poptype[pop],"_",sig,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(varsub,file=paste(dirstr,"var_FBSKAT_",vartype[fig],"_",poptype[pop],"_",sig,".txt",sep=""),row.names=FALSE,quote=FALSE)
    write.table(rep(1,length(wtssub)),file=paste(dirstr,"var1_FBSKAT_",vartype[fig],"_",poptype[pop],"_",sig,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE)
    
}

changepeds <- function(peds){
    
    n <- dim(peds)[2]
    misf <- rep(9,n)
    misf[5] <- 1
    mism <- rep(9,n)
    mism[5] <- 2
    peds0 <- c()
    for(i in 1:dim(peds)[1]){
        if(!( (peds[i,2] %in% union(peds[,3],peds[,4])) & peds[i,3]==9 & peds[i,4]==9)){
            onetrio <- peds[i,]
            if(peds[i,3]==9){
                tmpf <- misf
                tmpf[1] <- peds[i,1]
            }else{
                tmpf <- peds[peds[,2]==peds[i,3],]
            }
            
            if(peds[i,4]==9){
                tmpm <- mism
                tmpm[1] <- peds[i,1]
            }else{
                tmpm <- peds[peds[,2]==peds[i,4],]
            }
            
            onetrio <- rbind(onetrio,tmpf,tmpm)
            peds0 <- rbind(peds0,onetrio)
        }
    }
    
    peds0
}

subFBSKAT <- function(genos,fig,pop,peds,sex,phe,id){
    Z <- genos$Z
    fres <- genos$fres
    onelist <- genos$alllist
    
    ## sub columns corresponding to variants
    source("SKAT_ana.R")
    onelist <- subSKAT(onelist,fig,pop)
    vars <- paste(onelist[,1],onelist[,2],onelist[,4],onelist[,5],sep="_")
    univar <- unique(vars)
    subs <- match(univar,colnames(Z))
    Z <- Z[,univar]
    fres <- fres[subs]
    var_g <- cbind(onelist[,"Gene"],vars)
    
    Z1 <- c()
    Z1n <- c()
    genes <- unique(onelist[,"Gene"])
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
    fres <- fres[colnames(Z1)]
    onelist <- onelist[match(Z1n,vars),]  ## just use the first one
    
    ## sub rows corresponding to subjects
    bc.pop <- read.delim("data/WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    Jp <-  bc.pop[bc.pop[,4] %in% "J",3]
    Hp <-  bc.pop[bc.pop[,4] %in% "H",3]
    coln <- ifelse(any(grepl("SubID",colnames(onelist))), "SubID","Subject_ID")
    
    if(pop==1){ onep = Jp; }
    if(pop==2){ onep = Hp; }
    if(pop==3){ onep = c(Jp,Hp); }
    if(pop==4){ onep = id;}
    ## cannot change the orders
    Z1 <- Z1[rownames(Z1) %in% onep,]
    idsub <- intersect(onep,id)
    phe <- phe[id %in% idsub]
    peds <- peds[id %in% idsub,]
    sex <- sex[id %in% idsub]
    
    list(peds=peds,sex=sex,phe=phe,Z1=Z1,gst=gst,fres=fres,onelist=onelist)
    
}

qqplot_FBSKAT <- function(){
    ## qq plots for pvalues in variant and gene level
    library(Haplin)
    dirstr <- "FBSKATresult/FBout/"
    dirstr1 <- "FBSKATresult/qqout/"
    vartype <- c("LGD","D-mis","indels","LGD+D-mis","ALL")
    poptype <- c("Jewish","Hispanic","JH","All")
    for(flag in 1){
        for(sig in c(FALSE,TRUE)){
            for(pop in 1:4){
                for(fig in 1:5){
                    gfile <- paste(dirstr,"FBSKATr_",vartype[fig],"_",poptype[pop],"_",sig,"0.050000_0.010000_0.txt",sep="")
                    if(file.exists(gfile)){
                        x=read.table(gfile);
                        y1=1-pchisq(x[,4],x[,5]);
                        y2=1-pchisq(x[,6],x[,7]);
                        #y1[is.na(y1)] <- 1
                        #y2[is.na(y2)] <- 1
                        names(y1) <- x[,1]
                        names(y2) <- x[,1]
                        # del nan
                        y1 <- y1[!is.na(y1)]
                        y2 <- y2[!is.na(y2)]
                        pdf(file=paste(dirstr1,"FBSKATqq_",vartype[fig],"_",poptype[pop],"_",sig,".pdf",sep=""),width=20,height=10)
                        par(mfrow=c(1,1))
                        pQQ(y1, nlabs = sum(y1<0.05), conf = 0.95, mark = 0.05)
                        #pQQ(y2, nlabs = sum(y2<0.05), conf = 0.95, mark = 0.05)
                        dev.off()
                    }
                }
            }
        }
    }
    
}

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
