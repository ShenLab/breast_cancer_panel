batch_one <- function(){
        source("InheritedModels.R")
        source("misc.R")
	setwd("/home/local/ARCS/qh2159/breast_cancer/variants/families")
        ##### all families inherited models
#         L2CasesFams <- unlist(read.table("/home/local/ARCS/qh2159/breast_cancer/variants/families/FamiliesL2Cases.txt"))
#         aa <- read.delim("/home/local/ARCS/qh2159/breast_cancer/variants/families/Prioritized43families.txt")[,1]
#         Lars <- unique(read.delim("/home/local/ARCS/qh2159/breast_cancer/variants/families/LargeFamilyPhenotype.txt")[,1])
#         Lars <- setdiff(Lars,"200563")
#         svgfiles <- unique(paste(c(L2CasesFams,aa,Lars),".svg",sep=""))
#         wfile="Inheritance.Pattern.Families.v2.txt"
#         inheritedModels(svgfiles,wfile)
         
        
        #### get AD model frequency
	ncasef <- CaseinFam("Inheritance.Pattern.Families.v2.txt")
	PV <- control_freq(0.05)
	substr <- ".ADfiltered.tsv"
        Pfiles <- list.files("/home/local/ARCS/qh2159/breast_cancer/variants/families/Families78V2",pattern=".ADfiltered.tsv$",full.names = TRUE)
        ### delete the families without any reference sample
        tmp <- read.delim("Inheritance.Pattern.Families.Subjects.v2.txt",header=FALSE)
        tmp <- paste(tmp[tmp[,3]=="" & tmp[,2]=="AD",1],substr,sep="")
        ADfiles <- Pfiles[!(basename(Pfiles) %in% tmp)]
        wstr <- "FreqRef.v3.txt"
        Frequency_inherited(ADfiles,wstr,ncasef,substr,PV)
        
        wstr <- "FreqNonRef.v3.txt"
        Frequency_inherited(Pfiles,wstr,ncasef,substr,PV)
        
}

inheritedModels <- function(svgfiles,wfile){
        
        files <- list.files("/home/local/ARCS/qh2159/breast_cancer/variants/pedigree",pattern=".svg",full.names = TRUE,recursive = TRUE)
        tmpaa <- intersect(basename(files),svgfiles)
        files <- files[match(tmpaa,basename(files))]
        ### step 3: generate the AD and AR inherited model for each family
        pheno <- phenoin()
        pedis <- read.csv("/home/local/ARCS/qh2159/breast_cancer/variants/pedigree/ALL_pedigree.csv")
        
        fams <- gsub(".svg","",basename(files))
        n.fam <- length(fams)
        ADs <- matrix("",n.fam,8)
        ADs[,1] <- fams
        ADs[,2] <- "AD"
        ARs <- matrix("",n.fam,8)
        ARs[,1] <- fams
        ARs[,2] <- "AR"
        
        for(i in 1:length(fams)){
                onephe <- pheno[pheno[,1] %in% fams[i], ]
                oneped <- pedis[pedis[,1] %in% fams[i], ]
                
                ### AD model
                ## AD reference samples; heterozygous samples; NotFiltered Samples
                ## rule 1: Breast Cancer cases are heterozygous
                ## rule 2: Reference samples: No Breast Cancer and other Cancers; if there is any male, male older than 50 is ref; if there is any female, female older than 50 is ref.
                ## rule 3: samples not in rule 1 and rule 3 are not filtered samples
                tmp0 <- onephe[onephe[,"BreastCancer"]=="Yes", 2]
                ADs[i,4] <- paste(tmp0,sep="",collapse = ",")
                #tmp1tmp <- oneped[oneped[,3] %in% tmp0 | oneped[,4] %in% tmp0,2]
                tmp1 <- onephe[onephe[,"BreastCancer"]=="No" & onephe[,"CANCODE_1"]=="" & onephe[,"Sex"]=="Male" & onephe[,"UNIage"] >= 50, 2]
                #tmp <- onephe[onephe[,"BreastCancer"]=="No" & onephe[,"CANCODE_1"]=="" & onephe[,"Sex"]=="Male" & !(onephe[,2] %in% tmp1tmp), 2]
                tmp2 <- onephe[onephe[,"BreastCancer"]=="No" & onephe[,"CANCODE_1"]=="" & onephe[,"Sex"]=="Female" & onephe[,"UNIage"] >= 50, 2]
                ADs[i,3] <- paste(c(tmp1,tmp2),sep="",collapse = ",")
                ADs[i,8] <- paste(setdiff(onephe[,2],c(tmp0,tmp1,tmp2)),sep="",collapse = ",")
                
                ### AR model
                ## AR heterozygous samples; Alternate samples; NotAlternate samples; NotFiltered samples
                ## rule 1: Breast Cancer cases are alternate samples
                ## rule 2: Heterozygous samples: alternate samples' parents without Breast Cancers;  alternate samples' female children older than 50, alternate samples' male children, without any other type of cancers.
                ## rule 3: NotReference Sample IDs: samples with other type of cancers;
                ## rule 4: NotAlternate samples: alternate samples's siblings without any type of cancers
                ## rule 5: NotFiltered samples: samples not in rule 1,2 and 3,4
                tmp0 <- onephe[onephe[,"BreastCancer"]=="Yes", 2]
                ARs[i,5] <- paste(tmp0,sep="",collapse = ",")
                tmp1 <- intersect(onephe[onephe[,"BreastCancer"]=="No", 2], c(oneped[oneped[,2] %in% tmp0,3],oneped[oneped[,2] %in% tmp0,4]))
                tmp2 <- intersect(onephe[onephe[,"BreastCancer"]=="No" & onephe[,"CANCODE_1"]=="" & onephe[,"Sex"]=="Male", 2], oneped[oneped[,3] %in% tmp0 | oneped[,4] %in% tmp0,2])
                tmp3 <- intersect(onephe[onephe[,"BreastCancer"]=="No" & onephe[,"CANCODE_1"]=="" & onephe[,"Sex"]=="Female" & onephe[,"UNIage"] >= 50, 2], oneped[oneped[,3] %in% tmp0 | oneped[,4] %in% tmp0,2])
                ARs[i,4] <- paste(c(tmp1,tmp2,tmp3),sep="",collapse = ",")
                tmp4tmp <- c(oneped[oneped[,2] %in% tmp0,3], oneped[oneped[,2] %in% tmp0,4])
                #tmp4 <- intersect(onephe[onephe[,"BreastCancer"]=="No" & onephe[,"CANCODE_1"]!="", 2], oneped[oneped[,3] %in% tmp4tmp | oneped[,4] %in% tmp4tmp, 2])
                tmp4 <- onephe[onephe[,"BreastCancer"]=="No" & onephe[,"CANCODE_1"]!="", 2]
                ARs[i,6] <- paste(tmp4,sep="",collapse = ",")
                tmp5 <- intersect(onephe[onephe[,"BreastCancer"]=="No" & onephe[,"CANCODE_1"]=="", 2], oneped[oneped[,3] %in% tmp4tmp | oneped[,4] %in% tmp4tmp, 2])
                ARs[i,7] <- paste(tmp5,sep="",collapse = ",")
                tmp6 <- setdiff(onephe[,2],c(tmp0,tmp1,tmp2,tmp3,tmp4,tmp5))
                ARs[i,8] <- paste(tmp6,sep="",collapse = ",")
        }
        
        # write the AD and AR model for each family
        InModels <- matrix("",2*n.fam,8)
        InModels[seq(1,2*n.fam,2), ] <- ADs
        InModels[seq(2,2*n.fam,2), ] <- ARs
        qwt(InModels,file=wfile)
        
        ### step 4: mapping individual ids to subject ids
        b <- InModels
        for(i in 1:dim(b)[1]){
                for(j in 3:8){
                        if(b[i,j]!=""){
                                tmp <- unlist(strsplit(b[i,j],","))
                                b[i,j] <- paste(pheno[match(tmp,pheno[,2]),3],sep="",collapse=",")
                        }
                }
        }
        qwt(b,file=gsub("Families.","Families.Subjects.",wfile))
        
        ### run Ashley Inherited model script

}

Frequency_inherited <- function(ADfiles,wstr,ncasef,substr=".AD.tsv",PV=""){
        ### step 5: most frequency genes and variants in families
        aa <- ADfiles
        oneT <- c()
        for(i in 1:length(aa)){
                tmp <- read.delim(aa[i])
                if(dim(tmp)[1]>0){
                        TMP0 <- tmp[,grepl("X.",names(tmp))]
                        ssub <- grepl("GT",names(TMP0))
                        stmp <- sum(ssub)
                        tmpn <- gsub("\\D","",names(TMP0)[ssub])
                        nCarrier <- sapply(1:dim(TMP0)[1],function(kk) stmp-sum(grepl("0/0",TMP0[kk,ssub])) ) ## note: not 0/0 
                        Carriers <- sapply(1:dim(TMP0)[1],function(kk) paste(tmpn[!grepl("0/0",TMP0[kk,ssub])],sep="",collapse="_")  )
                        tmp <- tmp[,!grepl("X.",names(tmp))]
                        FAMID <- gsub(substr,"",basename(aa[i]))
                        tmp <- cbind(FAMID,nCarrier,Carriers,tmp)
                        oneT <- rbind(oneT,tmp)
                }
        }
        oneT[is.na(oneT[,"AlleleFrequency.ExAC"]),"AlleleFrequency.ExAC"] <- 0
        oneT[oneT[,"AlleleFrequency.ExAC"]==".","AlleleFrequency.ExAC"] <- 0
        oneT <- oneT[as.numeric(oneT[,"AlleleFrequency.ExAC"]) <= 0.01, ]
        vars <- paste(oneT[,"Chromosome"],oneT[,"Position"],sep="_")
        
        if(any(PV!="")){
                #source("/home/local/ARCS/qh2159/breast_cancer/Panel/breast_cancer_panel/sourcefiles.R")
                #AJs <- unlist(read.table(AJBRfile))
                #HIs <- unlist(read.table(AJBRfile))
                pheno <- phenoin()
                phe1 <- paste(pheno[,4],pheno[,5],sep="")
                phe1[phe1=="JH"] <- "J"
                oneT <- oneT[(oneT[,"FAMID"] %in% pheno[phe1=="J",1] & !(vars %in% PV[[1]])) | (oneT[,"FAMID"] %in% pheno[phe1=="H",1] & !(vars %in% PV[[2]])) | (oneT[,"FAMID"] %in% pheno[phe1=="",1]),  ]
        }
        
        ### order by variants
        varC <- sort(table(vars),decreasing = TRUE)
        varN <- names(varC)
        varOr <- c()
        for(i in 1:length(varC)){
                subs <- vars==varN[i]
                nCase <- ncasef[match(oneT[subs,"FAMID"],ncasef[,1]),2]
                nFam_Case_Carrier <- paste(varC[i],sum(nCase),sum(oneT[subs,"nCarrier"]),sep="_")
                oner <- cbind(nFam_Case_Carrier,nCase,oneT[subs, ])
                varOr <- rbind(varOr,oner)
        }
        qwt(varOr,file=paste("Variants",wstr,sep=""),flag=2)

        ### order by genes
        geneC <- sort(table(oneT[,"Gene"]),decreasing = TRUE)
        geneN <- names(geneC)
        geneOr <- c()
        for(i in 1:length(geneC)){
                onetmp <- oneT[oneT[,"Gene"]==geneN[i], ]
                #Freq <- geneC[i]
                nFam <- length(unique(onetmp[,"FAMID"]))
                nVar <- length(unique(onetmp[,"Position"]))
                tmpgs <- c()
                for(kk in 1:dim(onetmp)[1]){
                        tmpgs <- union(tmpgs,unlist(strsplit(onetmp[kk,"Carriers"],"_")))
                }
                nCase <- length(tmpgs)
                oner <- cbind(nFam,nVar,nCase,oneT[oneT[,"Gene"]==geneN[i], ])
                geneOr <- rbind(geneOr,oner)
        }
        
        qwt(geneOr,file=paste("Genes",wstr,sep=""),flag=2)
        
}

phenoin <- function(){
        pheno <- read.csv("/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/WES BCFR phenotypic data.csv")
        pheno[pheno[,3]=="220897, 220904",3] <- "220897"
        pheno[pheno[,3]=="222357, 222966",3] <- "222357" 
        subsb <- sapply(1:dim(pheno)[1],function(i) 114-as.numeric(unlist(strsplit(pheno[i,"BIRTHDT"],"/"))[3]) )
        UNIage <- sapply(1:dim(pheno)[1],function(i) if(is.na(pheno[i,"LiveAge"])){subsb[i];}else{as.numeric(pheno[i,"LiveAge"]);})
        pheno <- cbind(pheno,UNIage)
        pheno
}

CaseinFam <- function(wfile,m="AD"){
        inh <- read.delim(gsub("Families.","Families.Subjects.",wfile),header=FALSE)   
        inh <- inh[inh[,2]==m, ]
        ncasef <- matrix(0,dim(inh)[1],2)
        ncasef[,1] <- inh[,1]
        for(i in 1:dim(inh)[1]){
              ncasef[i,2] <- length(unlist(strsplit(inh[i,4],",")))  
        }
        ncasef
}

Variants_IGV <- function(){
        ### Only carriers
        aa <- list.files("./variants",pattern=".AD.tsv$",full.names = TRUE)
        wstr <- "IGV.v1.txt"
        oneT <- c()
        for(i in 1:length(aa)){
                tmp <- read.delim(aa[i],check.names = FALSE)
                if(dim(tmp)[1]>0){
                        subs <- grepl(".GT",colnames(tmp))
                        sams <- colnames(tmp)[subs]
                        sams <- gsub("\\D","",sams)
                        oner <- do.call(rbind,lapply(1:dim(tmp)[1], function(k) cbind(tmp[k,"Chromosome"],tmp[k,"Position"],sams[!grepl("0/0",tmp[k,subs])])))
                        oneT <- rbind(oneT,oner)
                }
        }
        qwt(oneT,file=paste("ADvariants",wstr,sep=""))
}

Variants_IGV_100 <- function(){
        ### Only variants occurred in at least two families with hihly frequency in 43 families with at least one reference samples
        aa <- read.delim("VariantsFreq.v2.txt")
        vars <- paste(aa[,"Chromosome"],aa[,"Position"],sep="_")
        vars <- unique(vars[duplicated(vars)])
        
        alligvs <- read.table("ADvariantsIGV.v1.txt")
        igvs <- paste(alligvs[,1],alligvs[,2],sep="_")
        subs <- igvs %in% vars
        oner <- alligvs[subs, ]
       
        qwt(oner,file="ADvariants_L2fams.txt")
}

control_freq <- function(Pcut=0.05){
        source("/home/local/ARCS/qh2159/breast_cancer/Panel/breast_cancer_panel/sourcefiles.R")
        n.cont <- c(557,341)
        PV <- lapply(1:2, function(i){
                load(contlistfs[i])
                vars <- paste(onelist[,"Chromosome"],onelist[,"Position"],sep="_")
                tmp <- table(vars)
                names(tmp)[tmp >= n.cont[i]*Pcut]
                })
        PV
}

filtered_LargeFam <- function(){
        
        source("InheritedModels.R")
        source("misc.R")
        mis <- "nonsynonymousSNV"
        PV <- control_freq(0.05)
        pheno <- phenoin()
        phe1 <- paste(pheno[,4],pheno[,5],sep="")
        phe1[phe1=="JH"] <- "J"
        path="/home/local/ARCS/qh2159/breast_cancer/variants/families/Families78V2"
        
        samf  <- list.files(path=path,pattern=".tsv$",full.names=TRUE)
        for(i in 1:length(samf)){
                oner <- read.delim(samf[i],check.names=FALSE)
                if(dim(oner)[1]>0){
                        oner <- variant_filtering(oner,mis,Ecut=0.01,segd=0.95,pp2=TRUE,hotf="",alleleFrefile="",popcut=0.05)
                        oner <- oner[oner[,"ExACfreq"] & oner[,"VCFPASS"] & oner[,"noneSegmentalDup"],]
                        ### further filtered by others
                        ## meta-svm D or CADD >= 15 or polyphen as D
                        subs <- rep(TRUE,dim(oner)[1])
                        subs[oner[,"VariantClass"] %in% mis] <- FALSE
                        oner <- oner[subs | (oner[,"MetaSVM"]=="D" | (oner[,"PP2prediction"]=="D" & oner[,"CADDscore"] >= 15)), ]
                        ### filtered by population control
                        famid <- gsub("\\D","",basename(samf[i]))
                        if(famid %in% pheno[phe1=="J",1]) oner <- oner[!(paste(oner[,"Chromosome"],oner[,"Position"],sep="_") %in% PV[[1]]),  ]
                        if(famid %in% pheno[phe1=="H",1]) oner <- oner[!(paste(oner[,"Chromosome"],oner[,"Position"],sep="_") %in% PV[[2]]),  ]
                        qwt(oner,file=gsub(".tsv","filtered.tsv",samf[i]),flag=2)
                }
        }
}
