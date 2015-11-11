phenoinfo <- function(){
    phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv"
    variantpath <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"
    
    files <- list.files(path=variantpath,pattern=".tsv$")
    subjects <- gsub(".AllVariants.tsv","",files)
    pheno <- read.csv(phenofile)
    pheno <- pheno[pheno[,3] %in% subjects,]
    
    pheno
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
    mp <- barplot(t(famdis[,c(4,3,2,1)]),ylim=c(0,max(rowSums(famdis))+5), space=0.4, col=c("green","cyan","blue","red") ,cex.axis=1.6,xlab="Number of family members",ylab="Number of families",cex.lab=2,legend = c("No Breast Cancer(Male)","No Breast Cancer(Female)","Breast Cancer(Male)","Breast Cancer(Female)"),main="Family Distribution",cex.main=2 )
    axis(1, at=mp, labels=tmp1,cex.axis=2)
    text(x=mp,y=rowSums(famdis),labels=testr,pos=3,cex=0.8)
    text(x=mp[5],y=max(rowSums(famdis))/2+50,labels="(#subject, #case(F), #case(M), #non-case(F), #non-case(M))",cex=1)
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
    
    Hpheno <- pheno[pheno[,"HISPFAM"] %in% "H",]
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
