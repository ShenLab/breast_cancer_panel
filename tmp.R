batch_one <- function(){
        source("tmp.R")
        source("~/.Rprofile")
        ##### all families inherited models
        L2CasesFams <- unlist(read.table("FamiliesL2Cases.txt"))
        aa <- read.delim("Prioritized43families.txt")[,1]
        Lars <- unique(read.delim("../LargeFamily_7/LargeFamilyPhenotype.txt")[,1])
        Lars <- setdiff(Lars,"200563")
        svgfiles <- unique(paste(c(L2CasesFams,aa,Lars),".svg",sep=""))
        wfile="Inheritance.Pattern.Families.v2.txt"
        inheritedModels(svgfiles,wfile)
        ncasef <- CaseinFam("Inheritance.Pattern.Families.v2.txt")
        
        ##### run batches
#         L2CasesFams <- unlist(read.table("FamiliesL2Cases.txt"))
#         ### batch one have run
#         aa <- read.delim("Prioritized43families.txt")
#         #svgfiles <- paste(aa[,1],".svg",sep="")
#         #wfile="Inheritance.Pattern.Prior.Families.v1.txt"
#         #inheritedModels(svgfiles,wfile)
#         ### batch two 
#         aa <- setdiff(L2CasesFams,aa[,1])
#         svgfiles <- paste(aa,".svg",sep="")
#         wfile="Inheritance.Pattern.L2.Families.v1.txt"
#         inheritedModels(svgfiles,wfile)

        #### get AD model frequency
        #Lfiles <- list.files("../LargeFamily_7/variants",pattern=".AD.tsv$",full.names = TRUE)
        Pfiles <- list.files("./variants",pattern=".AD.tsv$",full.names = TRUE)
        ### delete the families without any reference sample
        tmp <- read.delim("Inheritance.Pattern.Families.Subjects.v1.txt",header=FALSE)
        tmp <- paste(tmp[tmp[,3]=="" & tmp[,2]=="AD",1],".AD.tsv",sep="")
        ADfiles <- Pfiles[!(basename(Pfiles) %in% tmp)]
        wstr <- "Freq.v2.txt"
        Frequency_inherited(ADfiles,wstr,ncasef)
        
        ADfiles <- list.files("./variants",pattern=".AD.tsv$",full.names = TRUE)
        wstr <- "Freq.v1.txt"
        Frequency_inherited(ADfiles,wstr,ncasef)
        
}

inheritedModels <- function(svgfiles,wfile){
        ### step 1: copy svg files for families
        svgbatch1 <- list.files("../../Analysis_Wendy/150504 Family Pedigree/",pattern=".svg$",full.names = TRUE)
        svgbatch2 <- list.files("../Last43MissedFamilies/MissingFamilies/",pattern=".svg$",full.names = TRUE)
        svgbatch3 <- list.files("../LargeFamily_7/pedigrees/",pattern=".svg$",full.names = TRUE)
        famid1 <- basename(svgbatch1)
        famid2 <- basename(svgbatch2)
        famid3 <- basename(svgbatch3)
        fs <- intersect(svgfiles,c(famid1,famid2,famid3))
#         file.copy(svgbatch1[famid1 %in% fs],"./pedigree/")
#         file.copy(svgbatch2[famid2 %in% fs],"./pedigree/")
#         file.copy(svgbatch3[famid3 %in% fs],"./pedigree/")
        
#         print(length(intersect(svgfiles,c(famid1,famid2,famid3))))
#         print(length(intersect(svgfiles,famid3)))
        
        ### step 2: delete file DIV to make it open with inkscape
         files <- paste("./pedigree/",fs,sep="") #list.files("./pedigree",pattern=".svg$",full.names = TRUE)
#         for(i in 1:length(files)){
#                 con <- file(files[i],"r")
#                 tmp <- readLines(con,warn=FALSE)
#                 tmp <- gsub("<DIV id=\"svgContainer\">","",tmp)
#                 tmp <- gsub("</DIV>","",tmp)
#                 #writeLines(tmp,con=file(files[i],"w"))
#                 close(con)
#                 qwt(tmp,files[i])
#         }
        
        ### step 3: generate the AD and AR inherited model for each family
        pheno <- phenoinfo()
        pedis <- read.csv("../../ALL_pedigree.csv")
        
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

Frequency_inherited <- function(ADfiles,wstr,ncasef){
        ### step 5: most frequency genes and variants in families
        aa <- ADfiles
        oneT <- c()
        for(i in 1:length(aa)){
                tmp <- read.delim(aa[i])
                if(dim(tmp)[1]>0){
                        TMP0 <- tmp[,grepl("X.",names(tmp))]
                        ssub <- grepl("GT",names(TMP0))
                        stmp <- sum(ssub)
                        nCarrier <- sapply(1:dim(TMP0)[1],function(kk) stmp-sum(grepl("0/0",TMP0[kk,ssub])) ) ## note: not 0/0 are considered
                        tmp <- tmp[,!grepl("X.",names(tmp))]
                        FAMID <- gsub(".AD.tsv","",basename(aa[i]))
                        tmp <- cbind(FAMID,nCarrier,tmp)
                        oneT <- rbind(oneT,tmp)
                }
        }
        oneT[is.na(oneT[,"AlleleFrequency.ExAC"]),"AlleleFrequency.ExAC"] <- 0
        oneT[oneT[,"AlleleFrequency.ExAC"]==".","AlleleFrequency.ExAC"] <- 0
        oneT <- oneT[as.numeric(oneT[,"AlleleFrequency.ExAC"]) <= 0.01, ]
        
        ### order by genes
#         geneC <- sort(table(oneT[,"Gene"]),decreasing = TRUE)
#         geneN <- names(geneC)
#         geneOr <- c()
#         for(i in 1:length(geneC)){
#                 Freq <- geneC[i]
#                 oner <- cbind(Freq,oneT[oneT[,"Gene"]==geneN[i], ])
#                 geneOr <- rbind(geneOr,oner)
#         }
#         qwt(geneOr,file="FrequencyGenes.txt",flag=2)
        
        ### order by variants
        vars <- paste(oneT[,"Chromosome"],oneT[,"Position"],sep="_")
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
        
#         ## the number of variant in each gene
#         aa <- read.delim("FrequencyGenes.txt")
#         genes <- unique(aa[,"Gene"])
#         Ng <- rep(0,length(genes))
#         for(i in 1:length(genes)){
#                 Ng[i] <- length(unique(aa[aa[,"Gene"]==genes[i],"Position"]))        
#         }
#         qwt(cbind(genes,Ng),file="Gene_N_variant.txt")
        
}

phenoinfo <- function(){
        pheno <- read.csv("../LargeFamily_7/WES BCFR phenotypic data.csv")
        pheno[pheno[,3]=="220897, 220904",3] <- "220897"
        pheno[pheno[,3]=="222357, 222966",3] <- "222357" 
        subsb <- sapply(1:dim(pheno)[1],function(i) 114-as.numeric(unlist(strsplit(pheno[i,"BIRTHDT"],"/"))[3]) )
        UNIage <- sapply(1:dim(pheno)[1],function(i) if(is.na(pheno[i,"LiveAge"])){subsb[i];}else{as.numeric(pheno[i,"LiveAge"]);})
        pheno <- cbind(pheno,UNIage)
        pheno
}

CaseinFam <- function(wfiles,m="AD"){
        inh <- read.delim(gsub("v1.txt","Subjects.v1.txt",wfile),header=FALSE)   
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
