phenoinfo <- function(){
    phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv"
    variantpath <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"
    
    files <- list.files(path=variantpath,pattern=".tsv$")
    subjects <- gsub(".AllVariants.tsv","",files)
    pheno <- read.csv(phenofile)
    pheno <- pheno[pheno[,3] %in% subjects,]
    
    pheno
}

## label unknown samples by admixture model
labelUnknown <- function(){


}

## samples and phenotype informations
pheno_samples <- function(){
    outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"  
    outliers <- read.delim(outlierfile,header=FALSE)[,2]
    dups <- unlist(read.table("/home/local/ARCS/qh2159/breast_cancer/Panel/data/CGC_YALE_Rege_repeatedSamples.txt"))
    
    length(intersect(outliers,dups))
    
    phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv"
    variantpath <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"
    files <- list.files(path=variantpath,pattern=".tsv$")
    subjects <- gsub(".AllVariants.tsv","",files)
    pheno <- read.csv(phenofile)
    
    length(subjects)
    dim(pheno)[1]
    
    pheno <- pheno[pheno[,3] %in% subjects,]
    dim(pheno)[1]
    sum(outliers %in% pheno[,3])
    
    pop <- paste(pheno[,4],pheno[,5],sep="")
    pop[pop=="JH"] <- "J"
    table(pop)
    
    AJBRfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/tablesandlists/All_AJ_samples.list"
    HIBRfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/HispanicCases548.txt"
    
    AJs <- unlist(read.table(AJBRfile))
    table(pop[pheno[,3] %in% AJs])
    
    HIs <- unlist(read.table(HIBRfile))
    length(setdiff(pheno[,3],c(AJs,HIs)))
    
    phenofile <- "../data/phenotype/WES_BCFR_phenotypic_data-19062015.txt"
    pheno1 <- read.delim(phenofile)
    pop <- paste(pheno1[,4],pheno1[,5],sep="")
    pop[pop=="JH"] <- "J"
    table(pop)
    ## two pheno file table give the same anotation for pop
}

## pseudo-contorl phenotype
pseudoControl <- function(){
    pheno <- phenoinfo()
    subs <- pheno[,"Sex"]== "Female" & pheno[,"BreastCancer"]=="No"
    pheno <- pheno[subs,]
    ### using both Birth data and also liveage
    subsb <- sapply(1:dim(pheno)[1],function(i) as.numeric(unlist(strsplit(pheno[i,"BIRTHDT"],"/"))[3]) <= 44 )
    subs <- sapply(1:dim(pheno)[1],function(i) if(is.na(pheno[i,"LiveAge"])){subsb[i];}else{as.numeric(pheno[i,"LiveAge"]) >= 70;})
    pheno <- pheno[subs,]
    
    canfil <- read.csv("../data/phenotype/WES_CaseControl_PossibleControls_OtherCancer.csv")
    tmpid <- canfil[canfil[,"Control.Status"]=="N","Subject_ID"]    
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% tmpid),]
    
    ## excluding outliers
    outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"  ## outlier samples
    outliers <- unlist(read.table(outlierfile))
    outliers <- sapply(1:length(outliers),function(i) {tmp=unlist(strsplit(outliers[i],"_"));setdiff(gsub("\\D","",tmp),c("",paste("00",1:9,sep="")));})
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% outliers),]
    
    pheno[is.na(pheno)] <- ""
    write.csv(pheno,file="../data/phenotype/controls_Qiang_1_5_2016.csv",row.names=FALSE)
    
    pheno
}

pseudoControl_old <- function(){
    ## old version discard
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
    mp <- barplot(t(barM),ylim=c(0,max(rowSums(barM))+5), space=0.6, col=c("black","cyan","blue","red","yellow","green") ,cex.axis=1.6,xlab="Number of family members",ylab="Number of families",cex.lab=2 ,main="Number of cases in each family",cex.main=2)
    legend("topright",legend=c("0 case","1 case","2 case","3 case","4 case","5 case"), col=c("black","cyan","blue","red","yellow","green"),cex=2.2,pch=15)
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

singleVariantTest_Jewish_Hispanic <- function(){
    source("misc.R")
    source("src.R")
    Cohortfile <- "../data/Rdata/BreastCancer_VariantList_11_12"
    pseudoCont <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/controls_Qiang_1_5_2016.csv"
    
    ### Jewish single variant test
    caselistf <- "../data/Rdata/AJcaselist_12_17"
    contlistf <- "../data/Rdata/AJcontlist_12_17"
    indexf <- "../data/AJindexcases265.txt"
    vburdenfile <- "../resultf/variant_level.burden.txt"
    singleVariantTest(Cohortfile,pseudoCont,caselistf,contlistf,indexf,vburdenfile)
    
    ### Hispanic single variant test
    caselistf <- "../data/Rdata/HIcaselist_1_5"
    contlistf <- "../data/Rdata/HIcontlist_1_5"
    indexf <- "../data/HIindexcases138.txt"
    vburdenfile <- "../resultf/HISP_variant_level.burden.txt"
    n.case=164
    group="HISP548"
    singleVariantTest(Cohortfile,pseudoCont,caselistf,contlistf,indexf,vburdenfile,n.case,group)
}

singleVariantTest  <- function(Cohortfile,pseudoCont,caselistf,contlistf,indexf,vburdenfile,n.case=336,group="AJ716"){
    
    indexcases <- unlist(read.table(indexf))
    mis <- "nonsynonymousSNV"
    Ecut=0.01
    pheno <- phenoinfo()
    cases <- pheno[pheno[,"BreastCancer"]=="Yes",3]
    nons <- pheno[pheno[,"BreastCancer"]=="No",3]
    
    load(caselistf)
    caselist <- onelist
    indexlist <- caselist[caselist[,"Subject_ID"] %in% indexcases, ]
    nonslist <- caselist[caselist[,"Subject_ID"] %in% nons, ]
    caselist <- caselist[caselist[,"Subject_ID"] %in% setdiff(cases,indexcases), ]
    rm(onelist)
    load(contlistf)
    contlist <- onelist
    rm(onelist)
    
    indexlist <- variant_filtering(indexlist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf="",alleleFrefile=NULL,popcut=0.05)
    indexlist <- indexlist[indexlist[,"filtered"], ]
    contlist <- variant_filtering(contlist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf="",alleleFrefile=NULL,popcut=0.05)
    contlist <- contlist[contlist[,"filtered"], ]
    
    ## single variant burden test in index cases and controls
    variantTable <- burden_test(indexlist,contlist,flag=3,sig=FALSE)
    qwt(variantTable,file=vburdenfile,flag=2) 

    ## add annotation information 
    indexvars <- paste(indexlist[,1],indexlist[,2],indexlist[,4],indexlist[,5],sep="_") 
    burdenf <- read.delim(vburdenfile)
    burdenlist <- cbind(burdenf,indexlist[match(burdenf[,2],indexvars),])
    burdenlist <- burdenlist[,!(colnames(burdenlist) %in% "Subject_ID")]
    
    ## add index cases, non-index cases, non-index controls and control ID for each variants
    casevars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_") 
    contvars <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
    nonsvars <- paste(nonslist[,1],nonslist[,2],nonslist[,4],nonslist[,5],sep="_")
    id3set <- sapply(1:dim(burdenlist)[1],function(i){
	c(paste(unique(indexlist[indexvars == burdenlist[i,2],"Subject_ID"]),sep="",collapse="_"),paste(unique(caselist[casevars == burdenlist[i,2],"Subject_ID"]),sep="",collapse="_"),paste(unique(nonslist[nonsvars == burdenlist[i,2],"Subject_ID"]),sep="",collapse="_"),paste(unique(contlist[contvars == burdenlist[i,2],"Subject_ID"]),sep="",collapse="_"))	
	})
    id3set <- t(id3set)
    colnames(id3set) <- c("index_case","non-index_case","non_case_cohort","control")
    burdenlist <- cbind(burdenlist,id3set)
    qwt(burdenlist,file=gsub(".burden.txt","_burden_anno.txt",vburdenfile),flag=2)
    
    ## add case, non-case frequency in our breast cancer cohort
    n.case <- length(unique(indexlist[,"Subject_ID"])) + length(unique(caselist[,"Subject_ID"]))
    n.noncase <- length(unique(nonslist[,"Subject_ID"]))
    burdenlist <- read.delim(gsub(".burden.txt","_burden_anno.txt",vburdenfile))
    fre <- sapply(1:dim(burdenlist)[1],function(i) c( length(unlist(strsplit(burdenlist[i,"index_case"],"_"))) + length(unlist(strsplit(burdenlist[i,"non.index_case"],"_"))), length(unlist(strsplit(burdenlist[i,"non_case_cohort"],"_"))), n.case, n.noncase))
    fre <- t(fre)
    fre[,3] <- fre[,3] - fre[,1]
    fre[,4] <- fre[,4] - fre[,2]
    fre <- cbind(fre,(fre[,1]/fre[,3])/(fre[,2]/fre[,4]))
    ps <- sapply(1:dim(fre)[1], function(i) fisher.test(matrix(fre[i,1:4],2,2))$p.value )
    fre <- cbind(fre,ps)
    colnames(fre) <- c(paste("#case",group,sep=""),paste("#noncase",group,sep=""),"#noshotcase","#noshotnoncase","#foldcase-non","#p_case_noncase")
    burdenlist <- cbind(burdenlist,fre)
    qwt(burdenlist,file=gsub(".burden.txt","_burden_anno_Fre.txt",vburdenfile),flag=2)
    
    ## add pseudo controls information
    burdenlist <- read.delim(gsub(".burden.txt","_burden_anno_Fre.txt",vburdenfile),check.names=FALSE)
    load(Cohortfile)
    conts <- read.csv(pseudoCont)
    pcontlist <- onelist[onelist[,"Subject_ID"] %in% conts[,3],]
    
    pconvars <- paste(pcontlist[,1],pcontlist[,2],pcontlist[,4],pcontlist[,5],sep="_")
    pcset <- sapply(1:dim(burdenlist)[1],function(i){
        paste(unique(pcontlist[pconvars == burdenlist[i,2],"Subject_ID"]),sep="",collapse="_")	
    })
    
    n.pcont <- length(unique(pcontlist[,"Subject_ID"]))
    fre <- sapply(1:dim(burdenlist)[1],function(i) length(unlist(strsplit(pcset[i],"_"))))
    fre <- cbind(fre,(burdenlist[,paste("#case",group,sep="")]/n.case)/(fre/n.pcont))
    ps <- sapply(1:dim(fre)[1], function(i) fisher.test(matrix(c(burdenlist[i,paste("#case",group,sep="")],fre[i,1],n.case-burdenlist[i,paste("#case",group,sep="")],n.pcont-fre[i,1]),2,2))$p.value )
    fre <- cbind(fre,ps)
    fre <- cbind(pcset,fre)
    colnames(fre) <- c("pseducontrols","#pseducontrol","#foldcase-pseducontrols","#p_case_pseducontrols")
    burdenlist <- cbind(burdenlist,fre)
    qwt(burdenlist,file=gsub(".burden.txt","_burden_anno_Fre_Pseducont.txt",vburdenfile),flag=2)
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

VariantSta(contstaf,contstaf,indexcase,conts,paste(outputpath,"variantSta_overlap/",sep=""))

}

### families without cases in Jewish and Hispanic
FamNoCases <- function(){

source("src.R")
pheno <- phenoinfo()
fams <- setdiff(pheno[pheno[,"BreastCancer"]=="No",1],pheno[pheno[,"BreastCancer"]=="Yes",1])
pheno1 <- pheno[pheno[,1] %in% fams, ]
write.csv(pheno1,file="../data/ALL_cohort_withoutcases_families.csv",row.names=FALSE,quote=FALSE)

AJcasefile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/tablesandlists/All_AJ_samples.list"
AJs <- unlist(read.table(AJcasefile))
AJpheno <- pheno[pheno[,3] %in% AJs,]
pheno <- AJpheno
fams <- setdiff(pheno[pheno[,"BreastCancer"]=="No",1],pheno[pheno[,"BreastCancer"]=="Yes",1])
pheno <- phenoinfo()
pheno2 <- pheno[pheno[,1] %in% fams, ]
write.csv(pheno2,file="../data/Jewish716_withoutcases_families.csv",row.names=FALSE,quote=FALSE)

pheno <- phenoinfo()
pheno <- pheno[pheno[,"HISPFAM"]=="H" & pheno[,"AJFAM"]!="J", ]
fams <- setdiff(pheno[pheno[,"BreastCancer"]=="No",1],pheno[pheno[,"BreastCancer"]=="Yes",1])
pheno <- phenoinfo()
pheno3 <- pheno[pheno[,1] %in% fams, ]
write.csv(pheno3,file="../data/Hispanic548_withoutcases_families.csv",row.names=FALSE,quote=FALSE)
 
}
