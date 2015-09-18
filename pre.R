## one case from one family data 
caseonefam <- function(){
    pheno <- read.csv("WES BCFR phenotypic data.csv")
    
    ## excluded outlier samples and undetermined and mismatched sex samples
    outliers <- read.delim("Potential_Outliers.tsv")[,1]
    sexcheck <- read.delim("CUMC_Regeneron.mismatches.sexcheck.tsv")[,1]
    outlfam <- read.delim("Potential_Problem_families.tsv",sep=" ")[,1]
    
    outliers <- union(outliers,sexcheck)
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% outliers),]
    pheno <- pheno[!(pheno[,"FAMILYID"] %in% outlfam),]
    
    ## one case in one family
    famid <- unique(pheno[,1])
    indid <- unique(pheno[,1])
    
    nsubj <- dim(pheno)[1]
    nfam <- length(famid)
    subs <- rep(FALSE,nsubj)
    for(i in 1:nfam){
        tmpsub <- which(pheno[,1] %in% famid[i])
        casesub <- which(pheno[tmpsub,"BreastCancer"] == "Yes")
        onesub <- which.min(pheno[tmpsub[casesub],"LiveAge"])
        subs[tmpsub[casesub[onesub]]] <- TRUE
    }
    
    phenoin <- pheno[subs,1:15]
    
    phenoin
    
}

## contorl phenotype
controlpheno <- function(){

    pheno <- read.csv("WES BCFR phenotypic data.csv")
    
    ## excluded outlier samples and undetermined and mismatched sex samples
    outliers <- read.delim("Potential_Outliers.tsv")[,1]
    sexcheck <- read.delim("CUMC_Regeneron.mismatches.sexcheck.tsv")[,1]
    outlfam <- read.delim("Potential_Problem_families.tsv",sep=" ")[,1]
    
    outliers <- union(outliers,sexcheck)
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% outliers),]
    pheno <- pheno[!(pheno[,"FAMILYID"] %in% outlfam),]
    
    subs <- pheno[,"Sex"]== "Female" & pheno[,"BreastCancer"]=="No"
    pheno <- pheno[subs,]
    subs <- sapply(1:dim(pheno)[1],function(i) unlist(strsplit(pheno[i,"BIRTHDT"],"/"))[3]<= 44 )
    pheno <- pheno[subs,]

    canfil <- read.csv("WES_CaseControl_PossibleControls_OtherCancer.csv")
    tmpid <- canfil[canfil[,"Control.Status"]=="N","Subject_ID"]    
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% tmpid),]
    
    pheno[is.na(pheno)] <- ""
    write.csv(pheno,file="controls_Qiang.csv",row.names=FALSE)
    
    pheno
}

## case SKAT analysis 
caseSKAT <- function(sig,fig,pop){    
    source("Faminfo.R")
    load(paste("case3_single",sig,sep=""))
    caselist <- case3
    caselist <- subSKAT(caselist,fig,pop)
    
    ## MAF (minor-allele frequency): 
    freN <- c("ExAC.nfe.freq","ExAC.afr.freq")
    fres <- sapply(1:dim(caselist)[1], function(i){
        tmp <- unlist(strsplit(caselist[i,"INFO"],";"))
        tmp2 <- unlist(strsplit(tmp,"="))
        #print(as.numeric(tmp2[match(freN,tmp2)+1]))
        min(as.numeric(tmp2[match(freN,tmp2)+1]))
    })
    fres[is.na(fres)] <- 0
    vars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
    names(fres) <- vars
    vars <- unique(vars)
    fres <- fres[match(vars,names(fres))]
    #save(fres,file="MAF")
    
    ### genotype matrix in SKAT
    caseids <- unique(caselist[,"Subject_ID"])
    n.case <- length(caseids)
    n.var <- length(vars)
    Z <- matrix(0,n.case,n.var,dimnames=list(caseids,vars))
    for(i in 1:n.case){
        tmp <- caselist[caselist[,"Subject_ID"]==caseids[i],]
        svar <- paste(tmp[,1],tmp[,2],tmp[,4],tmp[,5],sep="_")
        geo <- rep(2,length(svar))
        geo[tmp[,"GT"]=="0/0"] <- 0
        geo[tmp[,"GT"]=="0/1"] <- 1
        Z[i,match(svar,vars)] <- geo
    }
    #save(Z,file="genotype")
    
    list(fres=fres,Z=Z,caselist=caselist)
}

## control SKAT analysis
controlSKAT <- function(sig,fig,pop){
#     source("pre.R")
#     control <- controlpheno()
# 
#     control_ID <- control[,"Subject_ID"]
#     contf <- paste(control_ID,".tsv",sep="")
#     #path="/ifs/scratch/c2b2/ys_lab/yshen/WENDY/BreastCancer/Regeneron/Filtering_for_Qiang/"
#     path="/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_for_Qiang/"
#     files <-  list.files(path=path,pattern=".tsv$")
#     contf <- intersect(contf,files)
#     
#     contlist <- c()
#     for(i in 1:length(contf)){
#         tmp <- paste(path,contf[i],sep="")
#         oner <- read.delim(tmp)
#         oner <- cbind(oner,gsub(".tsv","",contf[i]))
#         colnames(oner)[c(23,24,29,30)] <- c("GT","AD","Subject_INFO","SubID")
#         contlist <- rbind(contlist,oner)
#     }
#     save(contlist,file="contlist_8_6")
    
    source("Faminfo.R")
    load(paste("cont3_single",sig,sep=""))
    contlist <- cont3
    contlist <- subSKAT(contlist,fig,pop)
    
    ## MAF (minor-allele frequency): 
    freN <- c("ExAC.nfe.freq","ExAC.afr.freq")
    fres <- sapply(1:dim(contlist)[1], function(i){
        tmp <- unlist(strsplit(contlist[i,"INFO"],";"))
        tmp2 <- unlist(strsplit(tmp,"="))
        min(as.numeric(tmp2[match(freN,tmp2)+1]))
    })
    fres[is.na(fres)] <- 0
    vars <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
    names(fres) <- vars
    vars <- unique(vars)
    fres <- fres[match(vars,names(fres))]
    fres1 <- fres
    #save(fres1,file="MAF_cont")
    
    ### genotype matrix in SKAT
    contids <- unique(contlist[,"SubID"])
    n.cont <- length(contids)
    n.var <- length(vars)
    Z1 <- matrix(0,n.cont,n.var,dimnames=list(contids,vars))
    for(i in 1:n.cont){
        tmp <- contlist[contlist[,"SubID"]==contids[i],]
        svar <- paste(tmp[,1],tmp[,2],tmp[,4],tmp[,5],sep="_")
        geo <- rep(2,length(svar))
        geo[tmp[,"GT"]=="0/0"] <- 0
        geo[tmp[,"GT"]=="0/1"] <- 1
        Z1[i,match(svar,vars)] <- geo
    }
    #save(Z1,file="genotype_cont")
    
    list(fres1=fres1,Z1=Z1,contlist=contlist)
}

subSKAT <- function(onelist,fig,pop){

    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss","none")  
    oneind <- nchar(onelist[,"REF"]) != nchar(onelist[,"ALT"])
    
    #  fig=1; LOF only; fig=2; D-mis only; fig=3; indels only; fig=4; LOF and D-mis; fig=5; LOF, D-mis, indels
    if(fig==1){ onelist <- onelist[onelist[,"VariantClass"] %in% lof,];}
    if(fig==2){ onelist <- onelist[(onelist[,"VariantClass"] %in% mis) & !oneind,];}
    if(fig==3){ onelist <- onelist[(onelist[,"VariantClass"] %in% mis) & oneind,];}
    if(fig==4){ onelist <- onelist[onelist[,"VariantClass"] %in% lof | ((onelist[,"VariantClass"] %in% mis) & !oneind),];}
    #if(fig==5){ onelist <- onelist;}
    
    # # hp and jp: p=1 jp only; p=2; hp only; p=3 jp and hp; p=4 jp, hp and unknown;
    ### populations
    bc.pop <- read.delim("WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    
    Jp <-  bc.pop[bc.pop[,4] %in% "J",3]
    Hp <-  bc.pop[bc.pop[,4] %in% "H",3]
    coln <- ifelse(any(grepl("SubID",colnames(onelist))), "SubID","Subject_ID")
    if(pop==1){onelist <- onelist[onelist[,coln] %in% Jp,];}
    if(pop==2){onelist <- onelist[onelist[,coln] %in% Hp,];}
    if(pop==3){onelist <- onelist[onelist[,coln] %in% c(Jp,Hp),]}
    
    onelist
}

## single variants SKAT analysis
runSKAT <- function(sig,fig,pop){
    library(SKAT)
    
    pheno <- read.csv("WES BCFR phenotypic data.csv")
    r1 <- caseSKAT(sig,fig,pop)
    r2 <- controlSKAT(sig,fig,pop)
    
    Z <- r1$Z; fres <- r1$fres;
    Z1 <- r2$Z1; fres1 <- r2$fres1;
    
    sams <- union(rownames(Z),rownames(Z1))
    vars <- union(colnames(Z),colnames(Z1))
    
    fres[fres=="."] <- 0
    fres1[fres1=="."] <- 0
    fres <- as.numeric(fres)
    fres1 <- as.numeric(fres1)
    fres[is.na(fres)] <- 0
    fres1[is.na(fres1)] <- 0
    mafv <- rep(0,length(vars))
    mafv[match(colnames(Z),vars)] <- fres
    mafv[match(colnames(Z1),vars)] <- fres1
    olapr <- intersect(colnames(Z),colnames(Z1))
    mafv[match(olapr,vars)] <- apply(cbind(fres[match(olapr,colnames(Z))],fres1[match(olapr,colnames(Z1))]),1,min)
    
    wts <- dbeta(mafv,1,25)
    
    G <- matrix(0,length(sams),length(vars),dimnames=list(sams,vars))
    G[rownames(Z),colnames(Z)] <- Z
    G[rownames(Z1),colnames(Z1)] <- Z1
    
    X <- matrix(unlist(pheno[match(sams,pheno[,"Subject_ID"]),c("Sex","BIRTHDT")]),ncol=2)
    X[X[,1]=="Female",1] <- 0
    X[X[,1]=="Male",1] <- 1
    ###X[113,] <- unlist(pheno[512,c("Sex","LiveAge")]) #???
    X[,2] <- sapply(1:dim(X)[1],function(i) 114 - as.numeric(unlist(strsplit(X[i,2],"/"))[3]) )
    X <- as.numeric(X)
    X <- matrix(X,ncol=2)
    y <- rep(0,length(sams))
    y[match(rownames(Z1),sams)] <- 1
    
    #save(wts,file="weights")
    #save(G,file="Genotype")
    #save(X,file="X")
    #save(y,file="y")
    
    list(G=G,X=X,y=y,wts=wts,caselist=r1$caselist,contlist=r2$contlist)
    
    #### run SKAT
#     obj <- SKAT_Null_Model(y ~ X, out_type="D")
#     ## weighted
#     p0 <- SKAT(G, obj, kernel = "linear.weighted", weights=wts)
#     save(p0,file="SKATr0")
#     pV <- SKAT(G, obj)
#     save(pV,file="SKATr")
#     
#     # weighted rare variants
#     pV <- rep(1,dim(G)[2])
#     for(i in 1:dim(G)[2]){
#         pV[i] <- SKAT(as.matrix(G[,i],ncol=1), obj, kernel = "linear.weighted", weights=wts[i])$p.value
#     }
#     
#     save(pV,file="singleSKATr")
    
    ## weighted combined and rare variants
#     p0 <- SKAT_CommonRare(G, obj)$p.value
#     save(p0,file="SKATr_RC")
# 
#     pV <- rep(1,dim(G)[2])
#     for(i in 1:dim(G)[2]){
#         pV[i] <- SKAT_CommonRare(as.matrix(G[,i],ncol=1), obj)$p.value
#     }
#     
#     save(pV,file="singleSKATr_RC")
    
    
}

## combined variants and combined test of burden test and SKAT
comSKAT <- function(){
    library(SKAT)
    load("Genotype")
    load("weights")
    load("X")
    load("y")
    load("singleSKATr")
    
    obj <- SKAT_Null_Model(y ~ X, out_type="D")
    
    ## combined p.value < 0.05
    subs <- pV < 0.05
    p1 <- SKAT(G[,subs], obj, kernel = "linear.weighted", weights=wts[subs])$p.value ## 2.645224e-12
    
    ## combined p.value < 0.01
    subs <- pV < 0.01
    p2 <- SKAT(G[,subs], obj, kernel = "linear.weighted", weights=wts[subs])$p.value ## 6.001993e-10
    
    ## combined test of burden test and SKAT
    subs <- pV < 0.05
    p3 <- SKAT(G[,subs], obj, kernel = "linear.weighted", weights=wts[subs], method="optimal.adj") ## 7.626381e-18
    ## the optimal rho is 1, that is, the smallest p.value is based on the burden test, that means all the variants influence the phenotype in the same direction and with the same magnitude of effect. (Maybe because of LOF and Dmis)
    
    ## most related genes 
    source("~/.Rprofile")
    source("Faminfo.R")
    sig=FALSE
    load(paste("case3_single",sig,sep=""))
    caselist <- case3
    casevar <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
    
    subs <- pV < 0.05 ## 3970
    sivars <- colnames(G)[subs]
    
    #sum(sivars %in% casevar) ## 1811
    gesubs <- which(casevar %in% sivars)
    sigenes <- unique(caselist[gesubs,"Gene"])
    qwt(sigenes,file="effectgenes.txt")
    
    ## top genes
    subs <- pV < 0.01 ## 159
    sivars <- colnames(G)[subs]

    gesubs <- which(casevar %in% sivars)
    sigenes <- unique(caselist[gesubs,"Gene"])
    qwt(sigenes,file="topgenes.txt")
    
    ## top 50 varriants
    subs <- pV < 0.01
    sivars <- colnames(G)[subs]
    gesubs <- which(casevar %in% sivars)
    topvar <- caselist[gesubs,] 
    write.csv(topvar,file="topvars.csv",row.names=FALSE)
    
}

parallelSKAT <- function(){
    library(parallel)
    mclapply(1,function(kk) allSKAT(kk),mc.cores = 1)
}

allSKAT <- function(kk){
    source("pre.R")
    
    #  fig=1; LOF only; fig=2; D-mis only; fig=3; indels only; fig=4; LOF and D-mis; fig=5; LOF, D-mis, indels
    # hp and jp: p=1 jp only; p=2; hp only; p=3 jp and hp; p=4 jp, hp and unknown;
    if(kk < 21){
        fig <- floor((kk-1)/4) + 1
        pop <- ifelse(kk %% 4==0,4,kk %% 4)
        
        sig=FALSE
        r1 <- runSKAT(sig,fig,pop)
        comSKAT_gene(r1,sig,fig,pop,"SKATresult/")
        
        sig=TRUE
        r2 <- runSKAT(sig,fig,pop)
        comSKAT_gene(r2,sig,fig,pop,"SKATresult/")
    }
    
    if(kk==21){
        fig=5;pop=4;
        sig=FALSE
        r1 <- runSKAT(sig,fig,pop)
        comSKAT_variant(r1,sig,"SKATresult/")
        
        sig=TRUE
        r2 <- runSKAT(sig,fig,pop)
        comSKAT_variant(r2,sig,"SKATresult/")
    }
    
}

comSKAT_gene <- function(r1,sig,fig,pop,dirstr="SKATresult/"){
    library(SKAT)
    source("~/.Rprofile")
    source("Faminfo.R")
    G <- r1$G;wts <- r1$wts; X <- r1$X; y <- r1$y;caselist <- r1$caselist
    obj <- SKAT_Null_Model(y ~ X, out_type="D")
    casevar <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
    
    genes <- unique(caselist[,"Gene"])
    n.gene <- length(genes)
    print(n.gene)
    pM <- matrix(,n.gene,6)
    caselist[,c("SKAT_w","SKAT","SKAT-O_w","rho_w","SKAT-O","rho")] <- c("","","","","","")
    for(i in 1:n.gene){
        onevar <- casevar[which(caselist[,"Gene"] %in% genes[i])]
        subs <- which(colnames(G) %in% onevar)
        oneG <- as.matrix(G[,subs])
        pM[i,1] <- (SKAT(oneG, obj, kernel = "linear.weighted", weights=wts[subs])$p.value)
        pM[i,2] <- (SKAT(oneG, obj)$p.value)
        
        ## combined test of burden test and SKAT
        a1 <- SKAT(oneG, obj, kernel = "linear.weighted", weights=wts[subs], method="optimal.adj")
        a2 <- SKAT(oneG, obj, method="optimal.adj")
        pM[i,3] <- a1$p.value
        pM[i,4] <- ifelse(is.null(a1$param$rho_est), -1, a1$param$rho_est)
        pM[i,5] <- a2$p.value
        pM[i,6] <- ifelse(is.null(a2$param$rho_est), -1, a2$param$rho_est)
        
        n <- sum(caselist[,"Gene"]==genes[i])
        caselist[caselist[,"Gene"]==genes[i],c("SKAT_w","SKAT","SKAT-O_w","rho_w","SKAT-O","rho")] <- matrix(pM[i,],n,6,byrow=TRUE)
    }
    
    ## most related genes 
    source("~/.Rprofile")
    qwt(caselist,file=paste(dirstr,"caselist_SKAT_gene_",sig,"_",fig,"_",pop,".txt",sep=""),flag=2)
    
}

comSKAT_variant <- function(r1,sig,dirstr="SKATresult/"){
    library(SKAT)
    G <- r1$G;wts <- r1$wts; X <- r1$X; y <- r1$y;
    obj <- SKAT_Null_Model(y ~ X, out_type="D")
    
    source("~/.Rprofile")
    source("Faminfo.R")
    load(paste("case3_single",sig,sep=""))
    caselist <- case3
    casevar <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
    
    vars <- unique(casevar)
    n.var <- length(vars)
    pM <- matrix(,n.var,6)
    
    #Q_p = (1-rho)Q_s + rho*Q_B
    #Q_s is a test statistic of SKAT
    #Q_B is a score test statistic of the burden test
    caselist[,c("SKAT_w","SKAT","SKAT-O_w","rho_w","SKAT-O","rho")] <- c("","","","","","")
    for(i in 1:n.var){
        subs <- which(colnames(G)==vars[i])
        oneG <- as.matrix(G[,subs])
        pM[i,1] <- (SKAT(oneG, obj, kernel = "linear.weighted", weights=wts[subs])$p.value)
        pM[i,2] <- (SKAT(oneG, obj)$p.value)
        
        ## combined test of burden test and SKAT
        a1 <- SKAT(oneG, obj, kernel = "linear.weighted", weights=wts[subs], method="optimal.adj")
        a2 <- SKAT(oneG, obj, method="optimal.adj")
        pM[i,3] <- a1$p.value
        pM[i,4] <- ifelse(is.null(a1$param$rho_est), -1, a1$param$rho_est)
        pM[i,5] <- a2$p.value
        pM[i,6] <- ifelse(is.null(a2$param$rho_est), -1, a2$param$rho_est)
        
        n <- sum(casevar==vars[i])
        caselist[casevar==vars[i],c("SKAT_w","SKAT","SKAT-O_w","rho_w","SKAT-O","rho")] <- matrix(pM[i,],n,6,byrow=TRUE)
    }
    
    ## most related genes 
    source("~/.Rprofile")
    qwt(caselist,file=paste(dirstr,"caselist_SKAT_variant_",sig,".txt",sep=""),flag=2)
    
}

# burden test -------------------------------------------------------------
burdent <- function(){
    
    ## case not filtered by matched controls
    source("pre.R")
    cases <- caseonefam()
    case_ID <- cases[,"Subject_ID"]
    path="/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_for_Qiang/"
    casef <- paste(case_ID,".tsv",sep="")
       
    files <-  list.files(path=path,pattern=".tsv$")
    casef <- intersect(casef,files)
    casef <- setdiff(casef,"280002.tsv")#!!!!
    
    onlycase <- c()
    for(i in 1:length(casef)){
        tmp <- paste(path,casef[i],sep="")
        oner <- read.delim(tmp)
        oner <- cbind(oner,gsub(".tsv","",casef[i]))
        colnames(oner)[c(23,24,29,30)] <- c("GT","AD","Subject_INFO","Subject_ID")
        onlycase <- rbind(onlycase,oner)
    }
    save(onlycase,file="caselist_9_15")
    
    ### case filtered by match controls
    source("pre.R")
    path <- "/ifs/scratch/c2b2/ys_lab/yshen/WENDY/BreastCancer/Regeneron/CaseControl_Filtering"
    
    control <- controlpheno()
    control_ID <- control[,"Subject_ID"]
    control_ID <- paste("X",control_ID,".GT",sep="")
    
    cases <- caseonefam()
    case_ID <- cases[,"Subject_ID"]
    case_f <- cases[,1]
    casefile <- unlist(strsplit(paste(path,"/Fam_",case_f,".tsv",sep="",collapse = " ")," "))
    
    files <-  list.files(path=path,pattern=".tsv$",full.names=TRUE)
    subs <- which(casefile %in% files)
    casefile <- casefile[subs]
    case_ID <- case_ID[subs]
        
    ncase <- length(casefile)
    caselist <- c()
    cols <- c("Chromosome","Position","ID","REF","ALT","Gene","VariantFunction","VariantClass","AAchange","AlleleFrequency.ExAC","AlleleFrequency.1KG","AlleleFrequency.ESP","MetaSVM","SIFTprediction","PP2prediction","MAprediction","MTprediction","GERP..","CADDscore","SegmentalDuplication","PredictionSummary","VariantCallQuality","AlternateAlleles","MetaSVMScore","FILTER","INFO")
    
    for(i in 1:ncase){
        tmp <- read.delim(casefile[i])
        if(dim(tmp)[1] > 0){
            coln <- c(paste("X",case_ID[i],".GT",sep=""),paste("X",case_ID[i],".AD",sep=""),paste("X",case_ID[i],sep=""))
            subs <- tmp[,coln[1]]!="0/0" & tmp[,coln[1]]!="\\./\\."
            sub1 <- case_INFO(tmp[,"INFO"]) ## add filters the same with fileters in each sample
            if(sum(control_ID %in% colnames(tmp)) > 0){
                onec <- intersect(control_ID,colnames(tmp))
                sub2 <- rep(TRUE,length(subs))
                for(kk in 1:length(onec)){
                    subtmp <- tmp[,onec[kk]]=="0/0" | tmp[,onec[kk]]=="\\./\\."
                    sub2 <- sub2 & subtmp
                }
                subs <- subs & sub1 & sub2
            }else{
                subs <- subs & sub1
            }
            if(sum(subs)>0){
                tmp <- cbind(tmp[subs,cols],tmp[subs,coln],case_ID[i])
                if(any(is.na(tmp[,"Gene"]))) print(i)
                colnames(tmp)[27:30] <- c("GT","AD","INFO_VCF","Subject_ID")
                caselist <- rbind(caselist,tmp)
            }
        }
    }
    
    save(caselist,file="caselist_8_26")
    
    ##==============================================================================

}

# all population frequency
popvariant <- function(){
    pheno <- read.csv("WES BCFR phenotypic data.csv")
    
    ## excluded outlier samples and undetermined and mismatched sex samples
    outliers <- read.delim("Potential_Outliers.tsv")[,1]
    sexcheck <- read.delim("CUMC_Regeneron.mismatches.sexcheck.tsv")[,1]
    outlfam <- read.delim("Potential_Problem_families.tsv",sep=" ")[,1]
    
    outliers <- union(outliers,sexcheck)
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% outliers),]
    pheno <- pheno[!(pheno[,"FAMILYID"] %in% outlfam),]
    allf <- paste(pheno[,3],".tsv",sep="")
    allf <- gsub("222357, 222966.tsv","222357.tsv",allf)
    
    path="/ifs/scratch/c2b2/ys_lab/yshen/WENDY/BreastCancer/Regeneron/Filtering_for_Qiang_with_Synonymous/"
    files <-  list.files(path=path,pattern=".tsv$")
    allf <- intersect(allf,files)
    #allf <- setdiff(allf,"280002.tsv") ## double check 
    
    alllist <- c()
    for(i in 1:length(allf)){
        tmp <- paste(path,allf[i],sep="")
        oner <- read.delim(tmp)
        oner <- cbind(oner,gsub(".tsv","",allf[i]))
        colnames(oner)[c(24,25,30,45)] <- c("GT","AD","Subject_INFO","Subject_ID")
        alllist <- rbind(alllist,oner)
    }
    save(alllist,file="alllist_9_10")

    allV <- paste(alllist[,1],alllist[,2],alllist[,4],alllist[,5],sep="_")
    varT <- table(allV)
    save(varT,file="varT_9_10")
    
    #======================================================
    ### add synonmous variants
    load("caselist_9_15")
    caseid <- unique(onlycase[,"Subject_ID"])
    #load("caseid")
    load("contid")      
    path="/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_for_Qiang_with_Synonymous/"
    casef <- paste(caseid,".tsv",sep="")
    contf <- paste(contid,".tsv",sep="")
    syn <- c("","")
    
    casesy <- c()
    for(i in 1:length(casef)){
        tmp <- paste(path,casef[i],sep="")
        oner <- read.delim(tmp)
        oner <- cbind(oner,gsub(".tsv","",casef[i]))
        colnames(oner)[c(24,25,30,45)] <- c("GT","AD","Subject_INFO","Subject_ID")
        casesy <- rbind(casesy,oner)
    }
    
    casesy <- casesy[casesy[,"VariantClass"] %in% "synonymousSNV",]
    save(casesy,file="casesy_9_15")
    
    contsy <- c()
    for(i in 1:length(contf)){
        tmp <- paste(path,contf[i],sep="")
        oner <- read.delim(tmp)
        oner <- cbind(oner,gsub(".tsv","",contf[i]))
        colnames(oner)[c(24,25,30,45)] <- c("GT","AD","Subject_INFO","Subject_ID")
        contsy <- rbind(contsy,oner)
    }
    contsy <- contsy[contsy[,"VariantClass"] %in% "synonymousSNV",]
    save(contsy,file="contsy_9_15")

}

## INFO filtering 
case_INFO <- function(INFOs){
    #1000 genomes alternate allele frequency maximum: 0.01
    #ESP alternate allele frequency maximum: 0.01
    #Within VCF allele frequency maximum: 0.05
    
    cutN <- c(0.01,0.01,0.05)
    freN <- c("1KGfreq","ESPfreq","AF")
    sub1 <- sapply(1:length(INFOs),function(i) oneInfo(INFOs[i],freN,cutN))
    
    sub1
}

oneInfo <- function(Info,freN,cutN){

    tmp1 <- unlist(strsplit(Info,";"))
    tmp2 <- unlist(strsplit(tmp1,"="))
    
    sub1 = TRUE
    for(i in 1:length(freN)){
        if( is.na(match(freN[i],tmp2)) ){tmp <- TRUE;
        }else{
            va <- tmp2[match(freN[i],tmp2)+1];
            if(va==""){ tmp <- TRUE;
                        }else{
            va <- gsub("\\.,","0,",va)
            va <- gsub(",\\.",",0",va)
            va <- as.numeric(unlist(strsplit(va,",")))
            tmp <- max(va) < cutN[i]
                        }
        }
        sub1 <- sub1 & tmp    
    }
    
    sub1
}

# igvplot -----------------------------------------------------------------
igvplot <- function(){
    source("pre.R")

    pheno <- read.csv("WES BCFR phenotypic data.csv")
    ## excluded outlier samples and undetermined and mismatched sex samples
    outliers <- read.delim("Potential_Outliers.tsv")[,1]
    sexcheck <- read.delim("CUMC_Regeneron.mismatches.sexcheck.tsv")[,1]
    outlfam <- read.delim("Potential_Problem_families.tsv",sep=" ")[,1]
    
    outliers <- union(outliers,sexcheck)
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% outliers),]
    pheno <- pheno[!(pheno[,"FAMILYID"] %in% outlfam),]
    
#     load("Burden_caselist")
#     caselist <-caseL
#     a1 <- nchar(caselist[,"REF"])
#     a2 <- nchar(caselist[,"ALT"])
#     subs <- a1!=a2
#     indels <- caselist[subs,]
#     save(indels,file="indels_8_13")
    
    load("indels_8_13")
    con <- file("indels_IGV_8_13.txt","w")
    for(i in 1:dim(indels)[1]){
        famid <- pheno[pheno[,3]==indels[i,"Subject_ID"],1]
        Ss <- pheno[pheno[,1]==famid,3]
        if(length(Ss)>1) Ss <- paste(Ss,sep="",collapse=",")
        writeLines(paste(indels[i,1],indels[i,2],Ss,sep="\t"),con)
    }
    close(con)

    load("case3")
    a1 <- nchar(case3[,"REF"])
    a2 <- nchar(case3[,"ALT"])
    subs <- a1!=a2
    indels <- case3[subs,]
    save(indels,file="indels_8_24")
    
}

