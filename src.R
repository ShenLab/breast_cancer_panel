### === src.R
dupSamf <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/CGC_YALE_Rege_repeatedSamples.txt"
pContf <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/WES_CaseControl_PossibleControls_OtherCancer.csv" ### excldued related cancers not as controls
pseudoCont <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/controls_Qiang_1_5_2016.csv"
variantpath <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"

## label unknown samples by admixture model
labelUnknown <- function(){


}

## samples and phenotype informations
pheno_samples <- function(){
    outliers <- read.delim(outlierfile,header=FALSE)[,2]
    dups <- unlist(read.table(dupSamf))
    
    length(intersect(outliers,dups))
    
    files <- list.files(path=variantpath,pattern=".tsv$")
    subjects <- gsub(".AllVariants.tsv","",files)
    pheno <- read.csv(phenofile)
    pheno[pheno[,3]=="220897, 220904",3] <- "220897"
    pheno[pheno[,3]=="222357, 222966",3] <- "222357" 
    
    length(subjects)
    dim(pheno)[1]
    
    pheno <- pheno[pheno[,3] %in% subjects,]
    dim(pheno)[1]
    sum(outliers %in% pheno[,3])
    
    pop <- paste(pheno[,4],pheno[,5],sep="")
    pop[pop=="JH"] <- "J"
    table(pop)
     
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
    
    canfil <- read.csv(pContf)
    tmpid <- canfil[canfil[,"Control.Status"]=="N","Subject_ID"]    
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% tmpid),]
    
    ## excluding outliers
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
    ###dupIDs <- c("222357","222966") ## we keep 222357
    qwt(HIcases,file="../data/HispanicCases549.txt")
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
    
    AJs <- unlist(read.table(AJBRfile))
    AJpheno <- pheno[pheno[,3] %in% AJs,]
    FamilyDis(AJpheno,"../data/phenotype/JewishFamily.pdf","../data/phenotype/Jewishcase_perFamily.pdf")
    
    HIs <- unlist(read.table(HIBRfile))
    Hpheno <- pheno[pheno[,3] %in% HIs,]
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
        source("sourcefiles.R")
        source("misc.R")
        source("src.R")
        popG <- getPopGroup(AJBRfile,HIBRfile,phenofile,pseudoCont)
        
        ### Jewish single variant test
        caselistf <- "../data/Rdata/AJcaselist715_12_17"
        contlistf <- "../data/Rdata/AJcontlist_12_17"
        contlistf2 <- "../data/Rdata/HIcontlist_1_5"
        AJslistf <- "../data/Rdata/AJcaselist715_12_17"
        indexf <- "../data/AJindexcases265.txt"
        vburdenfile <- "../resultf/AJ_variant_level.burden.txt"
        NumOfSubjf <- "../resultf/NumofSubjects_AJ.txt"
        group="AJ715"
        singleVariantTest(caselistf,contlistf,indexf,Cohortfile,AJslistf,popG,contlistf2,vburdenfile,group,NumOfSubjf)
        
        ### Hispanic single variant test
        caselistf <- "../data/Rdata/HIcaselist_1_29"
        contlistf <- "../data/Rdata/HIcontlist_1_5"
        contlistf2 <- "../data/Rdata/AJcontlist_12_17"
        AJslistf <- "../data/Rdata/AJcaselist715_12_17"
        indexf <- "../data/HIindexcases138.txt"
        vburdenfile <- "../resultf/HISP_variant_level.burden.txt"
        NumOfSubjf <- "../resultf/NumofSubjects_HI.txt"
        group="HI550"
        singleVariantTest(caselistf,contlistf,indexf,Cohortfile,AJslistf,popG,contlistf2,vburdenfile,group,NumOfSubjf)
}

singleVariantTest  <- function(caselistf,contlistf,indexf,Cohortfile,AJslistf,popG,contlistf2,vburdenfile,group="AJ715",NumOfSubjf){
    ###======= burden test 
    mis <- "nonsynonymousSNV"
    Ecut=0.01
    indexcases <- unlist(read.table(indexf))
    exSamples <- excluded_samples()
    
    load(caselistf)
    onelist <- onelist[!(onelist[,"Subject_ID"] %in% exSamples), ]
    indexlist <- onelist[onelist[,"Subject_ID"] %in% indexcases, ]
    rm(onelist)
    
    load(contlistf)
    contlist <- onelist
    rm(onelist)
    
    indexlist <- variant_filtering(indexlist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf="",alleleFrefile="",popcut=0.05)
    indexlist <- indexlist[indexlist[,"filtered"], ]
    contlist <- variant_filtering(contlist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf="",alleleFrefile="",popcut=0.05)
    contlist <- contlist[contlist[,"filtered"], ]
    
    ## single variant burden test in index cases and controls
    variantTable <- burden_test(indexlist,contlist,flag=3,sig=FALSE)
    qwt(variantTable,file=vburdenfile,flag=2) 

    ###===============add annotation information=============
    indexvars <- paste(indexlist[,1],indexlist[,2],indexlist[,4],indexlist[,5],sep="_") 
    burdenf <- read.delim(vburdenfile,check.names=FALSE)
    burdenlist <- cbind(burdenf,indexlist[match(burdenf[,"Variant"],indexvars),])
    burdenlist <- burdenlist[,!(colnames(burdenlist) %in% c("Subject_ID"))]
    
    ### variant lists annotation
    load(AJslistf)
    onelist <- onelist[!(onelist[,"Subject_ID"] %in% exSamples), ]
    AJslist <- onelist
    rm(onelist)
    
    load(Cohortfile)
    onelist <- onelist[!(onelist[,"Subject_ID"] %in% exSamples), ]
    AJ.case <- length(intersect(onelist[,"Subject_ID"],popG[[1]]))
    AJ.pcont <- length(intersect(onelist[,"Subject_ID"],popG[[2]])) ## this is pseudo-controls
    HI.case <- length(intersect(onelist[,"Subject_ID"],popG[[5]]))
    HI.pcont <- length(intersect(onelist[,"Subject_ID"],popG[[6]])) ## this is pseudo-controls
    
    Anns <- c()
    Fres <- c()
    # popG: AJs: index cases, pseudo-controls, cases, non-cases; HIs: index cases, pseudo-controls, cases, non-cases;
    for(i in 1:8){
        if(i<=4){ onepop <- AJslist[AJslist[,"Subject_ID"] %in% popG[[i]], ];
        }else{ onepop <- onelist[onelist[,"Subject_ID"] %in% popG[[i]], ];}
        Freone <- onelistCount(onepop,burdenf[,"Variant"])
        Anns <- cbind(Anns,Freone[,1])
        Fres <- cbind(Fres,Freone[,2])
    }
    Frecont1 <- onelistCount(contlist,burdenf[,"Variant"])
    load(contlistf2)
    contlist2 <- onelist
    rm(onelist)
    Frecont2 <- onelistCount(contlist2,burdenf[,"Variant"])
    if(group=="AJ715"){
        Anns <- cbind(Anns,Frecont1[,1],Frecont2[,1])
        Fres <- cbind(Fres,Frecont1[,2],Frecont2[,2])
        AJ.cont <- length(unique(contlist[,"Subject_ID"]))
        HI.cont <- length(unique(contlist2[,"Subject_ID"]))
    }else if(group=="HI550"){
        Anns <- cbind(Anns,Frecont2[,1],Frecont1[,1])
        Fres <- cbind(Fres,Frecont2[,2],Frecont1[,2])
        AJ.cont <- length(unique(contlist2[,"Subject_ID"]))
        HI.cont <- length(unique(contlist[,"Subject_ID"]))
    }
    colnames(Fres) <- c("N.index.AJ","N.pseudoCont.AJ","N.case.AJ","N.non_case.AJ","N.index.HI","N.pseudoCont.HI","N.case.HI","N.non_case.HI","N.cont.AJ","N.cont.HI")
    colnames(Anns) <- c("index.AJ","pseudoCont.AJ","case.AJ","non_case.AJ","index.HI","pseudoCont.HI","case.HI","non_case.HI","cont.AJ","cont.HI")
    mode(Fres) <- "numeric"
    ### get odds and p values 
    tmp1 <- oneOddsPvalue(Fres[,c("N.index.AJ","N.pseudoCont.AJ")],AJ.case,AJ.pcont)
    tmp2 <- oneOddsPvalue(Fres[,c("N.index.AJ","N.cont.AJ")],AJ.case,AJ.cont)
    tmp3 <- oneOddsPvalue(Fres[,c("N.index.HI","N.pseudoCont.HI")],HI.case,HI.pcont)
    tmp4 <- oneOddsPvalue(Fres[,c("N.index.HI","N.cont.HI")],HI.case,HI.cont)
    OddPs <- cbind(tmp1,tmp2,tmp3,tmp4)
    colnames(OddPs) <- c("Odds_AJ_pseudo","p_AJ_pseudo","Odds_AJ_cont","p_AJ_cont","Odds_HI_pseudo","p_HI_pseudo","Odds_HI_cont","p_HI_cont")
        
    burdenlist <- cbind(burdenlist,Fres,OddPs,Anns)
    qwt(burdenlist,file=gsub(".burden.txt","_burden_Pseducont.txt",vburdenfile),flag=2)
    NumOfSubj <- paste(c("Num of AJ index cases:","Num of AJ pseudo-controls:","Num of AJ controls:","Num of HI cases:","Num of HI pseudo-controls:","Num of HI controls:"), c(AJ.case, AJ.pcont, AJ.cont, HI.case, HI.pcont, HI.cont),sep=" ")
    qwt(NumOfSubj, file=NumOfSubjf)
    
}

getPopGroup <- function(AJBRfile,HIBRfile,phenofile,pseudoCont){
    # popG: AJs: index cases, pseudo-controls, cases, non-cases; HIs: index cases, pseudo-controls, cases, non-cases;
    popG <- list()
    
    AJs <- unlist(read.table(AJBRfile))
    HIs <- unlist(read.table(HIBRfile))
    indexcases <- getindexcase(phenofile)
    pseudoconts <- read.csv(pseudoCont)[,3]
    pheno <- read.csv(phenofile)
    pheno[pheno[,3]=="220897, 220904",3] <- "220897"
    pheno[pheno[,3]=="222357, 222966",3] <- "222357" 

    cases <- pheno[pheno[,"BreastCancer"]=="Yes",3]
    conts <- pheno[pheno[,"BreastCancer"]=="No",3]
    
    popG[[1]] <- intersect(AJs,indexcases)
    popG[[2]] <- intersect(AJs,pseudoconts)
    popG[[3]] <- intersect(AJs,setdiff(cases,indexcases))
    popG[[4]] <- intersect(AJs,setdiff(conts,pseudoconts))
    popG[[5]] <- intersect(HIs,indexcases)
    popG[[6]] <- intersect(HIs,pseudoconts)
    popG[[7]] <- intersect(HIs,setdiff(cases,indexcases))
    popG[[8]] <- intersect(HIs,setdiff(conts,pseudoconts))    
    
    
    popG
}

onelistCount <- function(onelist,varlist){
    onevars <- paste(onelist[,1],onelist[,2],onelist[,4],onelist[,5],sep="_") 
    Freone <- sapply(1:length(varlist),function(i){
        tmp <- unique(onelist[onevars == varlist[i],"Subject_ID"])
        c(paste(tmp,sep="",collapse="_"),length(tmp))	
    })
    Freone <- t(Freone)
    Freone
}

oneOddsPvalue <- function(fre,n.case,n.cont){
    fre <- cbind(fre,n.case-fre[,1],n.cont-fre[,2])
    odds <- (fre[,1]/fre[,3])/(fre[,2]/fre[,4])
    ps <- sapply(1:dim(fre)[1], function(i) fisher.test(matrix(fre[i,1:4],2,2))$p.value )

    cbind(odds,ps)
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
    
    outliers <- read.delim(outlierfile,header=FALSE)[,2]
    
    onelist <- onelist[!(onelist[,"Subject_ID"] %in% outliers), ]
    onelist <- variant_filtering(onelist,mis,Ecut=0.01,segd=0.95,pp2=TRUE,hotf="",alleleFrefile="",popcut=0.05)
    onelist <- onelist[onelist[,"filtered"], ]
    
    qwt(Lpheno,file="../resultf/LargeFamilyPhenotype.txt",flag=2)
    qwt(onelist,file="../resultf/LargeFamilyVariants.txt",flag=2)
    
}

filtered_indexcases  <- function(){

        source("misc.R")
        source("sourcefiles.R")
        exSamples <- excluded_samples()
        
        load("../data/Rdata/AJcaselist715_12_17")
        indexcases <- getindexcase(phenofile)
        subjs <- unique(onelist[,"Subject_ID"])
        aa <- intersect(subjs,indexcases)
        aa <- setdiff(aa,exSamples)
        qwt(aa,file="../data/AJindexcases265.txt")
        
        load("../data/Rdata/HIcaselist_1_29")
        subjs <- unique(onelist[,"Subject_ID"])
        aa <- intersect(subjs,indexcases)
        aa <- setdiff(aa,exSamples)
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

source("sourcefiles.R")
source("misc.R")
pheno <- phenoinfo()
fams <- setdiff(pheno[pheno[,"BreastCancer"]=="No",1],pheno[pheno[,"BreastCancer"]=="Yes",1])
pheno1 <- pheno[pheno[,1] %in% fams, ]
write.csv(pheno1,file="../data/ALL_cohort_withoutcases_families.csv",row.names=FALSE,quote=FALSE)

AJs <- unlist(read.table(AJBRfile))
pheno <- pheno[pheno[,3] %in% AJs,]
fams <- setdiff(pheno[pheno[,"BreastCancer"]=="No",1],pheno[pheno[,"BreastCancer"]=="Yes",1])
pheno <- phenoinfo()
pheno2 <- pheno[pheno[,1] %in% fams, ]
write.csv(pheno2,file="../data/Jewish715_withoutcases_families.csv",row.names=FALSE,quote=FALSE)

HIs <- unlist(read.table(HIBRfile))
pheno <- pheno[pheno[,3] %in% HIs, ]
fams <- setdiff(pheno[pheno[,"BreastCancer"]=="No",1],pheno[pheno[,"BreastCancer"]=="Yes",1])
pheno <- phenoinfo()
pheno3 <- pheno[pheno[,1] %in% fams, ]
write.csv(pheno3,file="../data/Hispanic549_withoutcases_families.csv",row.names=FALSE,quote=FALSE)
 
}

###= families with more than two cases
FamL2Cases <- function(){
    source("src.R")
    source("misc.R")
    source("sourcefiles.R")
    pheno <- phenoinfo()
    AJs <- unlist(read.table(AJBRfile))
    AJpheno <- pheno[pheno[,3] %in% AJs,]
    HIs <- unlist(read.table(HIBRfile))
    Hpheno <- pheno[pheno[,3] %in% HIs,]
    
    pheno <- pheno[pheno[,"BreastCancer"]=="Yes",]
    AJpheno <- AJpheno[AJpheno[,"BreastCancer"]=="Yes",]
    Hpheno <- Hpheno[Hpheno[,"BreastCancer"]=="Yes",]
    
    a <- table(pheno[,1])
    aj <- table(AJpheno[,1])
    hi <- table(Hpheno[,1])
    
    fams76 <- names(a)[a>=2]
    famshi <- names(hi)[hi>=2]
    famsaj <- names(aj)[aj>=2]
    
    setdiff(fams76,c(famshi,famsaj))
    
    qwt(fams76,file="../data/phenotype/FamiliesL2Cases.txt")
    qwt(famshi,file="../data/phenotype/FamiliesL2Cases_HI.txt")
    qwt(famsaj,file="../data/phenotype/FamiliesL2Cases_AJ.txt")
    #200165 unknown
    #300369 relabeled as not AJ
}

## double check single gene
FamGenes <- function(){
    source("src.R")
    source("misc.R")
    source("sourcefiles.R")
    pheno <- phenoinfo()
    Samples <- read.delim(Samplelistfile)
    load(Cohortfile)
    
    genes <- matrix(0,4,3)
    genes[1,] <- c("MSH3","5","79950727")
    genes[2,] <- c("PARP4","13","25077802")
    genes[3,] <- c("POT1","7","124482897")
    genes[4,] <- c("SETD2","3","47162897")
    
    for(i in 1:4){
        SingleGene(genes[i,],onelist,pheno,Samples)
    }
}

SingleGene <- function(gene,onelist,pheno,Samples){
    
    subs <- onelist[,"Gene"]==gene[1] & onelist[,"Chromosome"]==gene[2] & onelist[,"Position"]==gene[3]
    tmp <- onelist[subs,]
    samid <- unique(tmp[,"Subject_ID"])
    fams <- unique(pheno[pheno[,3] %in% samid, 1])
    
    PheCarriers <- pheno[match(samid,pheno[,3]),c("FAMILYID","INDIVID","Subject_ID","BreastCancer","LiveAge","Sex")]
    PheCarriers <- cbind(PheCarriers,Samples[match(samid,Samples[,2]),3])
    colnames(PheCarriers) <- c("FAMILYID","INDIVID","Subject_ID","BreastCancer","LiveAge","Sex","Ethnicity")
    PheCarriers <- PheCarriers[order(PheCarriers[,"Ethnicity"],PheCarriers[,"BreastCancer"],PheCarriers[,"Sex"]),]
    qwt(PheCarriers,file=paste("../single_check/",gene[1],"inCarriers.txt",sep=""),flag=2)
    
    #ages <- sapply(1:dim(pheno)[1], function(i) 114 -  as.numeric(unlist(strsplit(pheno[i,"BIRTHDT"],"/"))[3]) )
    
    oneg <- c()
    for(i in 1:length(fams)){
        tmp <- rep(0,14)
        tmp[1] <- fams[i]
        onePhe <- pheno[pheno[,1]==fams[i],]
        samset <- list()
        samset[[1]] <- intersect(onePhe[onePhe[,"BreastCancer"]=="Yes",3],samid)
        samset[[2]] <- setdiff(onePhe[onePhe[,"BreastCancer"]=="Yes",3],samid)
        samset[[3]] <- intersect(onePhe[onePhe[,"BreastCancer"]=="No",3],samid)
        samset[[4]] <- setdiff(onePhe[onePhe[,"BreastCancer"]=="No",3],samid)
        for(j in 1:4){
            tmp[j+1] <- length(samset[[j]])   
            tmp[j+5] <- paste(samset[[j]],sep="",collapse="_")
            tmp[j+9] <- paste(onePhe[onePhe[,3] %in% samset[[j]],"LiveAge"],sep="",collapse="_")
        }
        
        tmp[14] <- paste(unique(Samples[Samples[,2] %in% onePhe[,3],3]),sep="",collapse="_")
        oneg <- rbind(oneg,tmp)
    }
    
    tmp <- rep(0,14)
    tmp[1] <- "All"
    for(i in 2:5){
        tmp[i] <- sum(as.numeric(oneg[,i]))
    }
    oneg <- rbind(oneg,tmp)
    
    tmp <- rep(0,14)
    tmp[1] <- "AJs"
    for(i in 2:5){
        tmp[i] <- sum(as.numeric(oneg[oneg[,14]=="AJ",i]))
    }
    oneg <- rbind(oneg,tmp)
    
    tmp <- rep(0,14)
    tmp[1] <- "Dominicans"
    for(i in 2:5){
        tmp[i] <- sum(as.numeric(oneg[oneg[,14]=="Dominican",i]))
    }
    oneg <- rbind(oneg,tmp)
    colnames(oneg) <- c("FamilyID","#BR+Var+","#BR+Var-","#BR-Var+","#BR-Var-","BR+Var+","BR+Var-","BR-Var+","BR-Var-","BR+Var+age","BR+Var-age","BR-Var+age","BR-Var-age","Ethnicity")
    oneg <- oneg[order(oneg[,"Ethnicity"]),]
    qwt(oneg,file=paste("../single_check/",gene[1],"inFamily.txt",sep=""),flag=2)
    
}

