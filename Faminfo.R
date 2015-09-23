## family information and population information
Faminfo <- function(){

    pheno <- read.csv("WES BCFR phenotypic data.csv")
    
    ## excluded outlier samples and undetermined and mismatched sex samples
    outliers <- read.delim("Potential_Outliers.tsv")[,1]
    sexcheck <- read.delim("CUMC_Regeneron.mismatches.sexcheck.tsv")[,1]
    outlfam <- read.delim("Potential_Problem_families.tsv",sep=" ")[,1]
    
    outliers <- union(outliers,sexcheck)
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% outliers),]
    pheno <- pheno[!(pheno[,"FAMILYID"] %in% outlfam),]    
    ### BRCA families
    br <- c(223109,223041,223275,260333,222968) 
    
    unif <- pheno[pheno[,3] %in% br,1]
    allm <- pheno[pheno[,1] %in% unif,3]
    
    
    load("alllist_9_8")
    subs <- (alllist[,"Subject_ID"] %in% allm) & (alllist[,"Gene"] %in% c("BRCA1","BRCA2"))
    tmplist <- alllist[subs,]
    
    phetmp <- pheno[match(tmplist[,"Subject_ID"],pheno[,3]),]
    
    allV <- cbind(phetmp,tmplist)
    qwt(allV,file="BRCAFam.txt",flag=2)

}

Famdis <- function(){

    pheno <- read.csv("WES BCFR phenotypic data.csv")
    ## excluded outlier samples and undetermined and mismatched sex samples
    outliers <- read.delim("Potential_Outliers.tsv")[,1]
    sexcheck <- read.delim("CUMC_Regeneron.mismatches.sexcheck.tsv")[,1]
    outlfam <- read.delim("Potential_Problem_families.tsv",sep=" ")[,1]
    
    outliers <- union(outliers,sexcheck)
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% outliers),]
    pheno <- pheno[!(pheno[,"FAMILYID"] %in% outlfam),]  
    
    
    ### case and control
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
    
    pdf(file="Family.pdf",width=12,height=10)
    par(mai=c(2,2,1,1))
    mp <- barplot(t(famdis[,c(4,3,2,1)]),ylim=c(0,max(rowSums(famdis))+5), space=0.4, col=c("green","cyan","blue","red") ,cex.axis=1.6,xlab="Number of family members",ylab="Number of families",cex.lab=2,legend = c("No(Male)","No(Female)","Breast Cancer(Male)","Breast Cancer(Female)") )
    axis(1, at=mp, labels=tmp1,cex.axis=2)
    text(x=mp,y=rowSums(famdis),labels=testr,pos=3,cex=0.8)
    text(x=mp[5],y=max(rowSums(famdis))/2+30,labels="(#subject, #case(F), #case(M), #non-case(F), #non-case(M))",cex=1)
    dev.off()   
    
}

write_ped <- function(){

    ## phenotype information
    pheno <- read.csv("WES BCFR phenotypic data.csv")
    ## ped file based on vcf file (vcf file from Ashley based on common SNPs)
    pedf <- "family/aa.ped"
    ## family information inferred from Primus
    firf <- "family/FirstDegreeTable_INDIVID.tsv"
    secf <- "family/SecondDegreeTable_INDIVID.tsv"
    firin <- read.delim(firf,sep="\t")
    sedin <- read.delim(secf,sep="\t")
    firin <- firin[!grepl("Missing",firin[,"IID"]),]
    sedin <- sedin[!grepl("Missing",sedin[,"IID"]),]
    
    firped <- "family/FirstDegree.ped"
    secped <- "family/SecondDegree.ped"
    con1 <- file(firped,"w");
    con2 <- file(secped,"w");
    
    con <- file(pedf,'r');
    k <- 1
    line = readLines(con,n=1);
    while(length(line)!=0){
        tmp <- unlist(strsplit(line,"\t"))
        tmp2 <- tmp
        if(tmp[2] %in% pheno[,3]){
            tmp[1] <- pheno[pheno[,3]==tmp[2],1]
            tmp[5] <- pheno[pheno[,3]==tmp[2],"Sex"]
            tmp[6] <- pheno[pheno[,3]==tmp[2],"BreastCancer"]
            tmp2[1] <- tmp[1]
            tmp2[5] <- tmp[5]
            tmp2[6] <- tmp[6]
            
            sub <- which(firin[,"IID"]==pheno[pheno[,3]==tmp[2],2])
            if(length(sub) > 0){
                sub <- sub[1]
            tmp[3] <- ifelse(grepl("Missing",firin[sub,"PID"]),0,firin[sub,"PID"])
            tmp[4] <- ifelse(grepl("Missing",firin[sub,"MID"]),0,firin[sub,"MID"])
            if(tmp[3]!=0){tmp[3] <- pheno[pheno[,2]==tmp[3],3];}
            }
            write(tmp,con1,append=TRUE,sep="\t")
            
            sub <- which(sedin[,"IID"]==pheno[pheno[,3]==tmp[2],2])
            if(length(sub) > 0){
                sub <- sub[1]
            tmp2[3] <- ifelse(grepl("Missing",sedin[sub,"PID"]),0,sedin[sub,"PID"])
            tmp2[4] <- ifelse(grepl("Missing",sedin[sub,"MID"]),0,sedin[sub,"MID"])
            if(tmp2[3]!=0){tmp2[3] <- pheno[pheno[,2]==tmp2[3],3];}
            }
            write(tmp2,con2,append=TRUE,sep="\t")
        }
        k <- k+1
        print(k)
        line = readLines(con,n=1);
    }
    close(con)
    close(con1)
    close(con2)
    
}
