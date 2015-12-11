phenoinfo <- function(){
    phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv"
    variantpath <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"
    
    files <- list.files(path=variantpath,pattern=".tsv$")
    subjects <- gsub(".AllVariants.tsv","",files)
    pheno <- read.csv(phenofile)
    pheno <- pheno[pheno[,3] %in% subjects,]
    
    pheno
}

getHispanicCases <- function(){
    source("src.R")
    source("misc.R")
    pheno <- phenoinfo()
    HIcases <- pheno[pheno[,"HISPFAM"]=="H" & pheno[,"AJFAM"]!="J",3]
    qwt(HIcases,file="../data/HispanicCases548.txt")
}

FamilyDis <- function(pheno,familyf,casef){

    ### case and non-case
    tmp <- table(pheno[,1])
    tmp1 <- sort(unique(tmp))
    famdis <- matrix(0,length(tmp1),4)
    teM <- matrix(0,length(tmp1),5)
    for(i in 1:length(tmp1)){
        tmpfam <- names(tmp)[tmp==tmp1[i]]
        teM[i,1] <- sum(pheno[,1] %in% tmpfam)
        teM[i,2] <- sum(pheno[,1] %in% tmpfam & pheno[,"BreastCancer"]=="Yes" & pheno[,"Sex"]=="Female")
        teM[i,3] <- sum(pheno[,1] %in% tmpfam & pheno[,"BreastCancer"]=="Yes" & pheno[,"Sex"]=="Male")
        teM[i,4] <- sum(pheno[,1] %in% tmpfam & pheno[,"BreastCancer"]=="No" & pheno[,"Sex"]=="Female")
        teM[i,5] <- sum(pheno[,1] %in% tmpfam & pheno[,"BreastCancer"]=="No" & pheno[,"Sex"]=="Male")
        
        a <- sum(tmp==tmp1[i]) / sum(pheno[,1] %in% tmpfam) 
        for(j in 1:4){
            famdis[i,j] <- a* teM[i,j+1]
        }
    }
    testr <- paste(teM[,1],teM[,2],teM[,3],teM[,4],teM[,5],sep=",")
    testr <- paste("(",testr,")",sep="")
    
    pdf(file=familyf,width=12,height=10)
    par(mai=c(2,2,1,1))
    mp <- barplot(t(famdis[,c(4,3,2,1)]),ylim=c(0,max(rowSums(famdis))+15), space=0.4, col=c("green","cyan","blue","red") ,cex.axis=1.6,xlab="Number of family members",ylab="Number of families",cex.lab=2,legend = c("No Breast Cancer(Male)","No Breast Cancer(Female)","Breast Cancer(Male)","Breast Cancer(Female)"),main="Family Distribution",cex.main=2 )
    axis(1, at=mp, labels=tmp1,cex.axis=2)
    text(x=mp,y=rowSums(famdis),labels=testr,pos=3,cex=0.8)
    text(x=mp[5]+5,y=max(rowSums(famdis))/2+10,labels="(#subject, #case(F), #case(M), #non-case(F), #non-case(M))",cex=1)
    dev.off()   
    
    ## the number of cases in each family
    fams <- unique(pheno[,1])
    numcases <- rep(0,length(fams))
    for(i in 1:length(fams)){
        numcases[i] <- sum(pheno[pheno[,1] %in% fams[i],"BreastCancer"]=="Yes")
    }
    
    tmp <- table(pheno[,1])
    ncol <- max(numcases)
    nrow <- sort(unique(tmp))
    barM <- matrix(0,length(nrow),ncol+1)
    testr <- 1:dim(barM)[1]
    for(i in 1:dim(barM)[1]){
        for(j in 1:dim(barM)[2]){
            onefams <- names(tmp)[tmp==nrow[i]]
            barM[i,j] <- sum(numcases[fams %in% onefams]==j-1)
        }
        testr[i] <- paste(barM[i,1:min(c(nrow[i]+1,dim(barM)[2]))],sep="",collapse=",")
    }
    testr <- paste("(",testr,")",sep="")
    
    pdf(file=casef,width=12,height=10)
    par(mai=c(2,2,1,1))
    mp <- barplot(t(barM),ylim=c(0,max(rowSums(barM))+5), space=0.6, col=c("black","cyan","blue","red","yellow","green") ,cex.axis=1.6,xlab="Number of family members",ylab="",cex.lab=2,legend = c("no case","1 case","2 case","3 case","4 case","5 case"),main="Number of cases in each family",cex.main=2)
    axis(1, at=mp, labels=nrow,cex.axis=2)
    text(x=mp,y=rowSums(barM),labels=testr,pos=3,cex=0.8)
    dev.off()
    
}

PopulationDis <- function(){
    source("src.R")
    pheno <- phenoinfo()
    familyf <- "../data/phenotype/Family.pdf"
    casef <- "../data/phenotype/case_perFamily.pdf"
    FamilyDis(pheno,familyf,casef)
    
    AJcasefile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/tablesandlists/All_AJ_samples.list"
    AJs <- unlist(read.table(AJcasefile))
    AJpheno <- pheno[pheno[,3] %in% AJs,]
    FamilyDis(AJpheno,"../data/phenotype/JewishFamily.pdf","../data/phenotype/Jewishcase_perFamily.pdf")
    
    Hpheno <- pheno[pheno[,"HISPFAM"] %in% "H" & pheno[,"AJFAM"]!="J",]
    FamilyDis(Hpheno,"../data/phenotype/HisFamily.pdf","../data/phenotype/Hiscase_perFamily.pdf")

}

SubjectBatches <- function(){
    source("misc.R")
    bampaths <- "/home/local/ARCS/yshen/mnt/BigData/WENDY/BreastCancer/Regeneron/bams"
    bams <- list.files(path=bampaths,pattern=".bam$")
    bams <- basename(bams)
    
    bats <- cbind(bams,"","")
    n <- length(bams)
    for(i in 1:n){
        if(grepl("BRC",bams[i])){
            bats[i,2] <- unlist(strsplit(bams[i],"\\."))[1]
            bats[i,3] <- "BRC"
        }
        if(grepl("^COL",bams[i])){
            bats[i,2] <- unlist(strsplit(bams[i],"_"))[3]
            bats[i,3] <- "Regeneron"
        }
        if(grepl("^M_",bams[i])){
            bats[i,2] <- gsub("M_","",unlist(strsplit(bams[i],"\\."))[1])
            bats[i,3] <- "CGC"
        }
        if(grepl("^Sample_",bams[i])){
            bats[i,2] <- unlist(strsplit(bams[i],"_"))[3]
            bats[i,3] <- "Yale"
        }
    }
    colnames(bats) <- c("bamfile","Subject_ID","Batches")
    qwt(bats,file="../data/phenotype/Sample_batches.txt")
    
}

pseudoControl <- function(){
    ## did not run this updated version
    pheno <- phenoinfo() ## update
    
    ## old version
    pheno <- read.csv("../data/phenotype/WES BCFR phenotypic data.csv")
    ## excluded outlier samples and undetermined and mismatched sex samples
    outliers <- read.delim("data/Potential_Outliers.tsv")[,1]
    sexcheck <- read.delim("data/CUMC_Regeneron.mismatches.sexcheck.tsv")[,1]
    outlfam <- read.delim("data/Potential_Problem_families.tsv",sep=" ")[,1]
    outliers <- union(outliers,sexcheck)
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% outliers),]
    pheno <- pheno[!(pheno[,"FAMILYID"] %in% outlfam),]

    ## pseudo controls definition
    subs <- pheno[,"Sex"]== "Female" & pheno[,"BreastCancer"]=="No"
    pheno <- pheno[subs,]
    subs <- sapply(1:dim(pheno)[1],function(i) as.numeric(unlist(strsplit(pheno[i,"BIRTHDT"],"/"))[3])<= 44 )
    pheno <- pheno[subs,]
    canfil <- read.csv("../data/phenotype/WES_CaseControl_PossibleControls_OtherCancer.csv")
    tmpid <- canfil[canfil[,"Control.Status"]=="N","Subject_ID"]    
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% tmpid),]
    pheno[is.na(pheno)] <- ""
    write.csv(pheno,file="../data/phenotype/controls_Qiang.csv",row.names=FALSE)
    
}

genelists <- function(){
## tumor suppressors 
ts2  <- unlist(read.table("../../genelist/Tumor_supressor/TS_Cosmic.txt")) ## tumor supressors
ts3 <- unlist(read.table("../data/hotspots/TScell_filtered.txt"))
ts <- union(ts2,ts3)
write.table(ts,file="../data/hotspots/TS_collected.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

ts <- gsub("Sep-05","SEPT5",ts);ts <- gsub("Sep-06","SEPT6",ts);ts <- gsub("Sep-09","SEPT9",ts);
qwt(ts,"../data/hotspots/Tumor_suppressors_11_11.txt")

}

singleVarianttest <- function(){
    
    caselistf <- "../data/Rdata/AJcaselist_11_9"
    contlistf <- "../data/Rdata/AJcontlist_11_9"
    indexf <- "../data/AJindexcases265.txt"
    indexcases <- unlist(read.table(indexf))
    mis <- "nonsynonymousSNV"
    Ecut=0.01
    pheno <- phenoinfo()
    cases <- pheno[pheno[,"BreastCancer"]=="Yes",3]
    
    load(caselistf)
    caselist <- onelist
    indexlist <- caselist[caselist[,"Subject_ID"] %in% indexcases, ]
    caselist <- caselist[caselist[,"Subject_ID"] %in% cases, ]
    rm(onelist)
    load(contlistf)
    contlist <- onelist
    rm(onelist)
    
    
    caselist <- variant_filtering(caselist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf="",alleleFrefile=NULL,popcut=0.05)
    caselist <- caselist[caselist[,"filtered"], ]
    
    indexlist <- variant_filtering(indexlist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf="",alleleFrefile=NULL,popcut=0.05)
    indexlist <- indexlist[indexlist[,"filtered"], ]
    
    contlist <- variant_filtering(contlist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf="",alleleFrefile=NULL,popcut=0.05)
    contlist <- contlist[contlist[,"filtered"], ]
    
    vburdenfile <- "../resultf/variant_level.burden.txt"
    variantTable <- burden_test(indexlist,contlist,flag=3,sig=FALSE)
    qwt(variantTable,file=vburdenfile,flag=2) 
    
    vburdenfile <- "../resultf/variant_level.burden_affected_control.txt"
    vT <- burden_test(caselist,contlist,flag=3,sig=FALSE)
    qwt(vT,file=vburdenfile,flag=2)

}

variantDis <- function(){
    source("misc.R")
    source("src.R")
    load("../data/Rdata/AJcaselist_11_9")
    phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv" ## phenotype file to get index cases only
    outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"  ## outlier samples
    caselist <- onelist
    rm(onelist)
    indexcases <- getindexcase(phenofile)
    caselist <- caselist[caselist[,"Subject_ID"] %in% indexcases, ] ## only index cases
    outliers <- unlist(read.table(outlierfile)) ## exclude outlier subjects
    outliers <- sapply(1:length(outliers),function(i) {tmp=unlist(strsplit(outliers[i],"_"));setdiff(gsub("\\D","",tmp),c("",paste("00",1:9,sep="")));})
    caselist <- caselist[!(caselist[,"Subject_ID"] %in% outliers), ]
    load("../data/Rdata/AJcontlist_11_9")
    contlist <- onelist
    varTypes <- variantDis_one(caselist,contlist)
    
    rm(caselist)
    rm(contlist)
    
    load("../resultf/caselist_singleton_0.01")
    load("../resultf/contlist_singleton_0.01")
    varTypes1 <- variantDis_one(caselist,contlist)
    
    #=======================================
    outputpath <- "../resultf/"
    qwt(varTypes,file=paste(outputpath,"Jewish_case_control_variants.txt",sep=""),flag=2)
    qwt(varTypes1,file=paste(outputpath,"Jewish_case_control_variants_filtered_singleton.txt",sep=""),flag=2)
 
    density_plots(varTypes[,3:7],varTypes[,1],"../resultf/plots_all/")
    density_plots(varTypes1[,3:7],varTypes1[,1],"../resultf/plots_singleton_0.01/")
    
}

variantDis_one <- function(caselist,contlist){

    VariantClass <- c(".","frameshiftdeletion","frameshiftinsertion","none","nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV","stopgain","stoploss","synonymousSNV","unknown")
    lof <- c("frameshiftdeletion","frameshiftinsertion","none","stopgain","stoploss",".")
    mis <- "nonsynonymousSNV"
    indel <- c("nonframeshiftdeletion","nonframeshiftinsertion")
    syn <- "synonymousSNV"
    unknown <- c("unknown")
    
    cases <- unique(caselist[,"Subject_ID"])
    conts <- unique(contlist[,"Subject_ID"])
    n.case <- length(cases)
    n.cont <- length(conts)
    varTypes <- matrix(0,n.case+n.cont,7)
    
    for(i in 1:n.case){
        onecase <- caselist[caselist[,"Subject_ID"]==cases[i], ]
        varTypes[i,3:7] <- c(sum(onecase[,"VariantClass"] %in% lof), sum(onecase[,"VariantClass"] %in% mis), sum(onecase[,"VariantClass"] %in% indel),sum(onecase[,"VariantClass"] %in% syn),sum(onecase[,"VariantClass"] %in% unknown))
    }
    for(i in 1:n.cont){
        onecont <- contlist[contlist[,"Subject_ID"]==conts[i], ]
        varTypes[i+n.case,3:7] <- c(sum(onecont[,"VariantClass"] %in% lof), sum(onecont[,"VariantClass"] %in% mis), sum(onecase[,"VariantClass"] %in% indel), sum(onecont[,"VariantClass"] %in% syn),sum(onecont[,"VariantClass"] %in% unknown))
    }
    
    varTypes[,1] <- "control"
    varTypes[1:n.case,1] <- "case"
    varTypes[1:n.case,2] <- cases
    varTypes[(n.case+1):(n.case+n.cont),2] <- conts
    colnames(varTypes) <- c("case/control","Subject_ID","#LOF","#MIS","#indel","#synonymous","#unknown")
    
    varTypes  
}

getCohortVariantlist <- function(){

    namestr=".AllVariants.tsv"
    path <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"
    samf <- list.files(path=path,pattern=".tsv$")
    
    onelist <- c()
    for(i in 1:length(samf)){
        tmp <- paste(path,samf[i],sep="")
        oner <- read.delim(tmp,check.names=FALSE)
        subj <- gsub(namestr,"",samf[i])
        oner <- cbind(oner,subj)
        cols <- colnames(oner)
        colsub <- c(which(grepl(paste(subj,".GT",sep=""),cols) | grepl(paste(toupper(subj),".GT",sep=""),cols)),which(grepl(paste(subj,".AD",sep=""),cols) | grepl(paste(toupper(subj),".AD",sep=""),cols)),which(subj==cols | paste("X",subj,sep="")==cols | toupper(subj)==cols), dim(oner)[2])
        colnames(oner)[colsub] <- c("GT","AD","Subject_INFO","Subject_ID")
        onelist <- rbind(onelist,oner)
    }
    
    save(onelist,file="../resultf/BreastCancer_VariantList_11_12")
    
}

getLargeFamily <- function(){
    source("misc.R")
    source("src.R")
    mis <- "nonsynonymousSNV"
    pheno <- phenoinfo()
    nfam <- table(pheno[,1])
    LargeFam <- names(nfam[nfam>=7])
    Lpheno <- pheno[pheno[,1] %in% LargeFam,]
    subjects <- unique(Lpheno[,3])
    load("../resultf/BreastCancer_VariantList_11_12")
    onelist <- onelist[onelist[,"Subject_ID"] %in% subjects,]
    
    outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"  ## outlier samples
    outliers <- unlist(read.table(outlierfile))
    outliers <- sapply(1:length(outliers),function(i) {tmp=unlist(strsplit(outliers[i],"_"));setdiff(gsub("\\D","",tmp),c("",paste("00",1:9,sep="")));})
    
    onelist <- onelist[!(onelist[,"Subject_ID"] %in% outliers), ]
    onelist <- variant_filtering(onelist,mis,Ecut=0.01,segd=0.95,pp2=TRUE,hotf="",alleleFrefile=NULL,popcut=0.05)
    onelist <- onelist[onelist[,"filtered"], ]
    
    qwt(Lpheno,file="../resultf/LargeFamilyPhenotype.txt",flag=2)
    qwt(onelist,file="../resultf/LargeFamilyVariants.txt",flag=2)
    
}

variantLargeFamily <- function(){
    # low in control data, high in family cases and low in unaffected memebers
    pedigree <- read.csv("/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/ALL_pedigree.csv")
    onelist <- read.delim("../resultf/LargeFamilyVariants.txt")
    pheno <- read.delim("../resultf/LargeFamilyPhenotype.txt")
    
    subAJ <- pheno[pheno[,"AJFAM"]=="J",3]
    contlistf <- "../data/Rdata/AJcontlist_11_9"
    wfile="../resultf/LargeFamilyJewish.txt"
    Variantcounting(pheno,subAJ,onelist,contlistf,wfile)
    
    subAJ <- pheno[pheno[,"HISPFAM"]=="H",3]
    contlistf <- "../data/Rdata/HIcontlist_11_20"
    wfile="../resultf/LargeFamilyHispanic.txt"
    Variantcounting(pheno,subAJ,onelist,contlistf,wfile)

}

Variantcounting <- function(pheno,subAJ,onelist,contlistf,wfile){
    source("misc.R")
    subAJf <- subAJ[pheno[match(subAJ,pheno[,3]),"BreastCancer"]=="Yes"]
    subAJu <- subAJ[pheno[match(subAJ,pheno[,3]),"BreastCancer"]=="No"]
    caseAJf <- onelist[onelist[,"Subject_ID"] %in% subAJf,]
    caseAJu <- onelist[onelist[,"Subject_ID"] %in% subAJu,]
    load(contlistf)
    contlist <- onelist
    
    fvar <- paste(caseAJf[,1],caseAJf[,2],caseAJf[,4],caseAJf[,5],sep="_")
    uvar <- paste(caseAJu[,1],caseAJu[,2],caseAJu[,4],caseAJu[,5],sep="_")
    convar <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
    
    ffvar <- table(fvar)
    fuvar <- table(uvar)
    fcvar <- table(convar)
    
    utmp <- setdiff(names(fuvar),names(ffvar))
    
    atmp <- setdiff(union(names(ffvar),names(fuvar)),names(fcvar))
    fcvar[atmp] <- 0 
    atmp <- setdiff(names(ffvar),names(fuvar))
    fuvar[atmp] <- 0
    
    
    unifvar <- unique(fvar)
    s <- (ffvar[unifvar]+1)/(fuvar[unifvar]+2*fcvar[unifvar]+1)
    s1 <- - fuvar[utmp] 
    aa <- cbind(c(unifvar,utmp),c(s,s1))  
    aa <- cbind(aa,0,0,0)
    n1 <- length(unifvar)
    n2 <- length(c(unifvar,utmp))
    aa[1:n1,3] <- ffvar[unifvar]
    aa[1:n1,4] <- fuvar[unifvar]
    aa[1:n1,5] <- fcvar[unifvar]
    aa[(n1+1):n2,4] <- fuvar[utmp]
    aa[(n1+1):n2,5] <- fcvar[utmp]

    wvar <- rbind(caseAJf,caseAJu)
    vars <- paste(wvar[,1],wvar[,2],wvar[,4],wvar[,5],sep="_")
    wvar <- cbind(aa[match(vars,aa[,1]),2:5],wvar)
    wvar <- wvar[order(-as.numeric(wvar[,1])),]
    colnames(wvar)[1:4] <- c("score","affect","unaffect","control")
    wvar <- wvar[,!(colnames(wvar) %in% c("filtered","ExACfreq","VCFPASS","noneSegmentalDup","meta.SVM_PP2","hotspot","alleleFre"))]
    wvar <- wvar[,c(dim(wvar)[2],1:(dim(wvar)[2]-1))]
    wvar <- cbind(pheno[match(wvar[,"Subject_ID"],pheno[,3]),c("FAMILYID","BreastCancer")],wvar)
    colnames(wvar)[1] <- "FamilyID"
    
    qwt(wvar,file=wfile,flag=2)

}

doublecheck <- function(){
    source("misc.R")
    source("src.R")
    varfile <- "../resultf/burdentest_FALSE_0.01_HMM_hotspots_11_12/Panel_genes_variantlist.txt"
    onevar <- read.delim(varfile)
    onevar <- onevar[onevar[,"Gene"]=="MSH3",]
    
    ## write indels
    oneids <- c()
    for(i in 2:5){
        tmp <- unlist(strsplit(onevar[1,i],":"))
        oneids <- union(oneids,tmp)
    }
    onevarf <- matrix(0,length(oneids),3)
    onevarf[,1] <- onevar[1,"Chromosome"]
    onevarf[,2] <- onevar[1,"Position"]
    onevarf[,3] <- oneids
    qwt(onevarf,file="../single_check/MSH3_IGV.txt")
        
    pheno <- phenoinfo()
    onephe <- pheno[pheno[,3] %in% oneids,]
    qwt(onephe,file="../single_check/MSH3_pheno.txt",flag=2)
    
    pedigree <- read.csv("/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/ALL_pedigree.csv") # ALLadd_pedigree.csv
    indids <- pheno[pheno[,3] %in% oneids,2]
    oneped <- pedigree[pedigree[,2] %in% indids | pedigree[,3] %in% indids | pedigree[,4] %in% indids,]
    oneped <- cbind(oneped,oneped[,2:4])
    oneped[oneped[,5] %in% pheno[,2],5] <- pheno[match(oneped[oneped[,5] %in% pheno[,2],5],pheno[,2]), 3]
    oneped[oneped[,6] %in% pheno[,2],6] <- pheno[match(oneped[oneped[,6] %in% pheno[,2],6],pheno[,2]), 3]
    oneped[oneped[,7] %in% pheno[,2],7] <- pheno[match(oneped[oneped[,7] %in% pheno[,2],7],pheno[,2]), 3]
    qwt(oneped,file="../single_check/MSH3_pedigree.txt",flag=2)
    
    
    pheno <- phenoinfo()
    affs <- pheno[pheno[,"BreastCancer"]=="Yes",3]
    uffs <- pheno[pheno[,"BreastCancer"]=="No",3]
    AJcasefile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/tablesandlists/All_AJ_samples.list"
    HIcasefile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/HispanicCases548.txt"
    Ajs <- unlist(read.table(AJcasefile))
    His <- unlist(read.table(HIcasefile))
    
    fams <- pheno[match(oneids,pheno[,3]),1]
    allmems <- pheno[pheno[,1] %in% fams,3]
    tyvar <- read.delim("../single_check/MSH3_detail.txt",header=FALSE)
    tys <- c("ins","substitute","complex","substitute/del")
    numT <- matrix(0,4,4)
    conl <- setdiff(allmems,oneids)
    conNum <- c(sum(conl %in% intersect(affs,Ajs)),sum(conl %in% intersect(uffs,Ajs)),sum(conl %in% intersect(affs,His)),sum(conl %in% intersect(uffs,His)))
    for (i in 1:4){
        numT[i,] <- c(sum(oneids %in% intersect(intersect(affs,Ajs),tyvar[tyvar[,4] %in% tys[i],3])),sum(oneids %in% intersect(intersect(uffs,Ajs),tyvar[tyvar[,4] %in% tys[i],3])),sum(oneids %in% intersect(intersect(affs,His),tyvar[tyvar[,4] %in% tys[i],3])),sum(oneids %in% intersect(intersect(uffs,His),tyvar[tyvar[,4] %in% tys[i],3])) )
    }
    
}

depth_variant <- function(){

	AJindexf <- "/home/local/ARCS/qh2159/breast_cancer/variants/AJdepth/indexAJ.ldepth.mean"
	AJcontf <- "/home/local/ARCS/qh2159/breast_cancer/variants/AJdepth/AJ.ldepth.mean"

	HIindexf <- "/home/local/ARCS/qh2159/breast_cancer/variants/HIdepth/indexHI.ldepth.mean"
	HIcontf <- "/home/local/ARCS/qh2159/breast_cancer/variants/HIdepth/HI.ldepth.mean"
	
	onedepthcompare(AJindexf,AJcontf)
	onedepthcompare(HIindexf,HIcontf)
}

onedepthcompare <- function(AJindexf,AJcontf){

  	indexf <- read.delim(AJindexf)
        contf  <- read.delim(AJcontf)
        caseV <- paste(indexf[,1],indexf[,2],sep="_")
        contV <- paste(contf[,1],contf[,2],sep="_")


	interV <- intersect(caseV,contV)
	unicase <- indexf[!(caseV %in% interV),]
	unicont <- contf[!(contV %in% interV),]


        chr=c(1:22,"X")
        for(i in 1:length(chr)){

                a1= unicase[unicase[,1]==chr[i],2]
                a2= unicont[unicont[,1]==chr[i],2]

                a1M <- matrix(a1,ncol=length(a2),nrow=length(a1),byrow=FALSE)
                a2M <- matrix(a2,ncol=length(a2),nrow=length(a1),byrow=TRUE)
                print(sum(abs(a1M-a2M)<=30))

        }
}

filtered_indexcases  <- function(){

source("misc.R")

load("../data/Rdata/AJcaselist_11_9")
phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv"
indexcases <- getindexcase(phenofile)
subjs <- unique(onelist[,"Subject_ID"])
aa <- intersect(subjs,indexcases)

outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"  ## outlier samples
outliers <- unlist(read.table(outlierfile))
outliers <- sapply(1:length(outliers),function(i) {tmp=unlist(strsplit(outliers[i],"_"));setdiff(gsub("\\D","",tmp),c("",paste("00",1:9,sep="")));})
aa <- setdiff(aa,outliers)

BRCA1_2pathogenicfile <- "../data/phenotype/BRCA1_2.txt" 
tmp <- read.delim(BRCA1_2pathogenicfile)
pathogenic_sample <- tmp[tmp[,1] %in% c("likely pathogenic","pathogenic"),"Subject_ID"]
aa <- setdiff(aa,pathogenic_sample)
qwt(aa,file="../data/AJindexcases265.txt")

load("../data/Rdata/HIcaselist_11_20")
subjs <- unique(onelist[,"Subject_ID"])
aa <- intersect(subjs,indexcases)
aa <- setdiff(aa,outliers)
aa <- setdiff(aa,pathogenic_sample)
qwt(aa,file="../data/HIindexcases138.txt")

}

oneVariantSta <- function(){

source("misc.R")
outputpath="../resultf/"
casestaf <- "../data/BR_20151209.hardfiltered.stats_hq.tsv"
contstaf <- "../data/contAJ_20151209.hardfiltered.stats_hq.tsv"
indexcase <- unlist(read.table("../data/AJindexcases265.txt"))
conts <- unlist(read.table("/home/local/ARCS/qh2159/breast_cancer/variants/data/AJs_557.txt"))

VariantSta(casestaf,contstaf,indexcase,conts,paste(outputpath,"variantSta_overlap/",sep=""))

}

interested_genes <- function(){

genes <- c("MSH3","PARP4","PTPRF","ARID1B","POT1","SETD2","FBXW7","HOXD11")
varlist1 <- "../resultf/burdentest_FALSE_0.01_HMM_hotspots_11_12/Panel_genes_variantlist.txt"
varlist2 <- "../resultf/burdentest_FALSE_0.01_HMM_hotspots_11_12/Variantlist.txt"

var1 <- read.delim(varlist1)
var2 <- read.delim(varlist2)

vars <- rbind(var1,var2)
testg <- c()
testv <- c()
for(i in 1:length(genes)){
	tmps <- which(vars[,"Gene"]==genes[i])
	testv <- c(testv,vars[tmps,"Variant"])
	testg <- c(testg,rep(genes[i],length(tmps)))
}

phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv"
load("../resultf/BreastCancer_VariantList_11_12")
pheno <- phenoinfo()
onevar <- paste(onelist[,1],onelist[,2],onelist[,4],onelist[,5],sep="_")

for(i in 1:length(testv)){
	oneVariant(testv[i],testg[i],onelist,onevar,pheno)
}

}

oneVariant <- function(var,gene,onelist,onevar,pheno){
	pops <- paste(pheno[,"AJFAM"],pheno[,"HISPFAM"],sep="")
	samids <- onelist[onevar==var,"Subject_ID"]
	famids <- pheno[match(samids,pheno[,3]),1]
	indids <- pheno[match(samids,pheno[,3]),2]
	popids <- pops[match(samids,pheno[,3])]
	affids <- pheno[match(samids,pheno[,3]),"BreastCancer"]

	oneC <- matrix(0,length(famids),4)
	for(i in 1:length(famids)){
		onefams <- pheno[pheno[,1]==famids[i],3]
		oneC[i,1] <- sum((samids %in% onefams) & affids=="Yes")
		oneC[i,2] <- sum((samids %in% onefams) & affids=="No")
		oneC[i,3] <- sum(pheno[match(onefams,pheno[,3]),"BreastCancer"]=="Yes")
		oneC[i,4] <- sum(pheno[match(onefams,pheno[,3]),"BreastCancer"]=="No")
	}

	oneT <- cbind(famids,indids,samids,popids,affids,oneC)
	oneT[is.na(oneT)] <- ""
	colnames(oneT) <- c("FAMILYID","INDIVIDUALID","SUBJECTID","POP","AFFSTATUS","#case-with","#control-with","#CASE","#CONTROL")
	qwt(oneT,file=paste("../single_check/",gene,"_",var,".txt",sep=""),flag=2)

	## output IGV indels
	tmp <- unlist(strsplit(var,"_"))
	indf <- cbind(tmp[1],tmp[2],samids)
	qwt(indf,file=paste("../single_check/",gene,"_",var,"_IGVs.txt",sep=""))
}

