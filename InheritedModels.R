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
        
	frecols <- c("N.index.AJ","N.pseudoCont.AJ","N.case.AJ","N.non_case.AJ","N.cont.AJ","N.index.HI","N.pseudoCont.HI","N.case.HI","N.non_case.HI","N.cont.HI")
	samcols <- c("index.AJ","pseudoCont.AJ","case.AJ","non_case.AJ","cont.AJ","index.HI","pseudoCont.HI","case.HI","non_case.HI","cont.HI")
	oddcols <- c("Odds_AJ_pseudo","p_AJ_pseudo","Odds_AJ_cont","p_AJ_cont","Odds_HI_pseudo","p_HI_pseudo","Odds_HI_cont","p_HI_cont")
	cols <- c(frecols,oddcols,samcols)
	HIvars <- read.delim("/home/local/ARCS/qh2159/breast_cancer/Panel/resultf/HISP_variant_level_burden_Pseducont.txt")
	AJvars <- read.delim("/home/local/ARCS/qh2159/breast_cancer/Panel/resultf/AJ_variant_level_burden_Pseducont.txt")
	
	ncasef <- CaseinFam("Inheritance.Pattern.Families.v2.txt")
	PV <- control_freq(0.05)
	substr <- ".ADfiltered.tsv"
        Pfiles <- list.files("/home/local/ARCS/qh2159/breast_cancer/variants/families/Families78V2",pattern=".ADfiltered.tsv$",full.names = TRUE)
        tmp <- read.delim("Inheritance.Pattern.Families.Subjects.v2.txt",header=FALSE)
        tmp <- paste(tmp[tmp[,3]=="" & tmp[,2]=="AD",1],substr,sep="") ### delete the families without any reference sample
        ADfiles <- Pfiles[!(basename(Pfiles) %in% tmp)]
        ADfiles1 <- Pfiles[(basename(Pfiles) %in% tmp)]
        
        #### get AD model frequency
        wstr <- "Families_Reference_42.v3.txt"
        Frequency_inherited(ADfiles,wstr,ncasef,substr,PV)
        wstr <- "All_Families_77.v3.txt"
        Frequency_inherited(Pfiles,wstr,ncasef,substr,PV)
        wstr <- "Families_without_Reference_35.v3.txt"
        Frequency_inherited(ADfiles1,wstr,ncasef,substr,PV)
        
        #### get AD model Statistics
        wstr1 <- "Families_ReferenceSta_42.v3.txt"
        Statistic_inherited(ADfiles,wstr1,ncasef,substr,PV,HIvars,AJvars,cols)
        wstr1 <- "All_FamiliesSta_77.v3.txt"
        Statistic_inherited(Pfiles,wstr1,ncasef,substr,PV,HIvars,AJvars,cols)
        wstr1 <- "Families_without_ReferenceSta_35.v3.txt"
        Statistic_inherited(ADfiles1,wstr1,ncasef,substr,PV,HIvars,AJvars,cols)
        
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
        
        pheno <- phenoin()
        if(any(PV!="")){
                #source("/home/local/ARCS/qh2159/breast_cancer/Panel/breast_cancer_panel/sourcefiles.R")
                #AJs <- unlist(read.table(AJBRfile))
                #HIs <- unlist(read.table(AJBRfile))
                phe1 <- paste(pheno[,4],pheno[,5],sep="")
                phe1[phe1=="JH"] <- "J"
                oneT <- oneT[(oneT[,"FAMID"] %in% pheno[phe1=="J",1] & !(vars %in% PV[[1]])) | (oneT[,"FAMID"] %in% pheno[phe1=="H",1] & !(vars %in% PV[[2]])) | (oneT[,"FAMID"] %in% pheno[phe1=="",1]),  ]
        }
        
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
                nCase <- length(intersect(tmpgs,pheno[pheno[,"BreastCancer"]=="Yes",3]))
                oner <- cbind(nFam,nVar,nCase,oneT[oneT[,"Gene"]==geneN[i], ])
                geneOr <- rbind(geneOr,oner)
        }
        
        qwt(geneOr,file=paste("Genes",wstr,sep=""),flag=2)
        
}

Statistic_inherited <- function(ADfiles,wstr,ncasef,substr=".AD.tsv",PV="",HIvars,AJvars,cols){

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
        
        kcut <- 0.0625
        library(kinship)
        pedis <- read.csv("/home/local/ARCS/qh2159/breast_cancer/variants/pedigree/ALL_pedigree.csv")
        K <- kinship(pedis[,2],pedis[,3],pedis[,4])
        pheno <- phenoin()
        pheno <- pheno[pheno[,2] %in% colnames(K), ]
        subs <- match(pheno[,2],colnames(K))
        Kp <- K[subs,subs]
        colnames(Kp) <- pheno[,3]
        rownames(Kp) <- pheno[,3]
        Kp[lower.tri(Kp,diag=TRUE)] <- 0
        #pop <- paste(pheno[,4],pheno[,5],sep="")
        #pop[pop=="JH" | pop=="HJ"] <- "J"
        #mutR <- read.csv("/home/local/ARCS/qh2159/breast_cancer/Panel/data/PCGCSCI.csv")
        
        ### order by variants
        vars <- paste(oneT[,"Chromosome"],oneT[,"Position"],oneT[,"REF"],oneT[,"ALT"],sep="_")
        varC <- sort(table(vars),decreasing = TRUE)
        varN <- names(varC)
        varOr <- c()
        for(i in 1:length(varC)){
                subs <- vars==varN[i]
                nCase <- ncasef[match(oneT[subs,"FAMID"],ncasef[,1]),2]
                nFam_Case_Carrier <- paste(varC[i],sum(nCase),sum(oneT[subs,"nCarrier"]),sep="_")
                
                tmp <- oneT[subs, ]
                fams <- unique(tmp[,"FAMID"])
                fsam <- pheno[pheno[,1] %in% fams,3]
                cases <- pheno[pheno[,1] %in% fams & pheno[,"BreastCancer"]=="Yes",3]
                non_cases <- pheno[pheno[,1] %in% fams & pheno[,"BreastCancer"]=="No" & pheno[,"UNIage"] > 50, 3]
                carriers <- c()
                for(ii in 1:dim(tmp)[1]) carriers <- union(carriers, unlist(strsplit(tmp[ii,"Carriers"],"_")) )
                non_carriers <- setdiff(fsam,carriers)
                
                aS <- intersect(cases,carriers)
                a <- ifelse(length(aS) - sum(Kp[aS,aS] > kcut) >= 0, length(aS) - sum(Kp[aS,aS] > kcut), 1)
                bS <- intersect(non_cases,carriers)
                b <- ifelse(length(bS) - sum(Kp[bS,bS] > kcut) >= 0, length(bS) - sum(Kp[bS,bS] > kcut), 1)
                cS <- intersect(non_carriers,cases)
                c <- ifelse(length(cS) - sum(Kp[cS,cS] > kcut) >= 0, length(cS) - sum(Kp[cS,cS] > kcut), 1)
                dS <- intersect(non_carriers,non_cases)
                d <- ifelse(length(dS) - sum(Kp[dS,dS] > kcut) >= 0, length(dS) - sum(Kp[dS,dS] > kcut), 1)
                pvalue <- fisher.test(matrix(c(a,c,b,d),2,2))$p.value
        
                popN <- c(HIvars[HIvars[,"Variant"]==varN[i],cols], AJvars[AJvars[,"Variant"]==varN[i],cols])
                popN <- matrix(popN,dim(tmp)[1],length(cols)*2,byrow=TRUE)
                
                oner <- cbind(pvalue,nFam_Case_Carrier,nCase,tmp,popN)
                varOr <- rbind(varOr,oner)
                
        }
        
        varOr <- as.matrix(varOr[order(as.numeric(varOr[,"pvalue"])), ])
        qwt(varOr,file=paste("Variants",wstr,sep=""),flag=2)
        
        ### order by genes
        geneC <- sort(table(oneT[,"Gene"]),decreasing = TRUE)
        geneN <- names(geneC)
        geneOr <- c()
        for(i in 1:length(geneC)){
                onetmp <- oneT[oneT[,"Gene"]==geneN[i], ]
                nFam <- length(unique(onetmp[,"FAMID"]))
                nVar <- length(unique(onetmp[,"Position"]))
                
                onevar <- paste(onetmp[,"Chromosome"],onetmp[,"Position"],onetmp[,"REF"],onetmp[,"ALT"],sep="_")
                carriers <- c()
                popN <- c()
                for(kk in 1:dim(onetmp)[1]){
                        carriers <- union(carriers,unlist(strsplit(onetmp[kk,"Carriers"],"_")))
                        popN <- rbind(popN,c(HIvars[HIvars[,"Variant"]==onevar[kk],cols], AJvars[AJvars[,"Variant"]==onevar[kk],cols]))
                }
                
                fams <- unique(onetmp[,"FAMID"])
                fsam <- pheno[pheno[,1] %in% fams,3]
                cases <- pheno[pheno[,1] %in% fams & pheno[,"BreastCancer"]=="Yes",3]
                non_cases <- pheno[pheno[,1] %in% fams & pheno[,"BreastCancer"]=="No" & pheno[,"UNIage"] > 50, 3]
                n.case <- ifelse(length(cases) - sum(Kp[cases,cases] > kcut) >= 0, length(cases) - sum(Kp[cases,cases] > kcut), 1)
                n.cont <- ifelse(length(non_cases) - sum(Kp[non_cases,non_cases] > kcut) >= 0, length(non_cases) - sum(Kp[non_cases,non_cases] > kcut), 1)
                aS <- intersect(cases,carriers)
                a <- ifelse(length(aS) - sum(Kp[aS,aS] > kcut) >= 0, length(aS) - sum(Kp[aS,aS] > kcut), 1)
                bS <- intersect(non_cases,carriers)
                b <- ifelse(length(bS) - sum(Kp[bS,bS] > kcut) >= 0, length(bS) - sum(Kp[bS,bS] > kcut), 1)
                pvalue <- ifelse( (a+b)>0, binom.test(a,a+b,n.case/(n.case+n.cont))$p.value,1)

                oner <- cbind(pvalue,nFam,nVar,length(aS),onetmp,popN)
                geneOr <- rbind(geneOr,oner)
        }
        geneOr <- as.matrix(geneOr[order(as.numeric(geneOr[,"pvalue"])), ])
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

caselistFams <- function(){

        caselist <- onlyFilterCaselist(Ecut=0.01,hotf=4,swi=1)
        caselist2 <- onlyFilterCaselist(Ecut=0.01,hotf=4,swi=2)   
        save(caselist,file="/home/local/ARCS/qh2159/breast_cancer/variants/trios/caselist")
        save(caselist2,file="/home/local/ARCS/qh2159/breast_cancer/variants/trios/caselist2")
        caselistA <- rbind(caselist,caselist2)       
        caselistA 
}

onlyFilterCaselist <- function(Ecut=0.01,hotf=4,swi=1){
        source("misc.R")
        source("sourcefiles.R")
        
        if(hotf==1){ hotspotfile <- hotHMM; ## HongjianPred hotspots file
        }else if(hotf==2){ hotspotfile <- hotCOSMIC; ## COSMIC hotspots file
        }else if(hotf==3){ hotspotfile="";
        }else if(hotf==4){ hotspotfile <- hotHMM;}
        ## switch to the right files for burden test
        alleleFrefile <- alleleFrefiles[swi] 
        caselistf <- caselistfs[swi]
        mis <- "nonsynonymousSNV"
        load(caselistf)
        caselist <- onelist
        rm(onelist)
        exSamples <- excluded_samples() ## exclude subjects
        caselist <- caselist[!(caselist[,"Subject_ID"] %in% exSamples), ]
        caselist <- variant_filtering(caselist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf=hotspotfile,alleleFrefile,popcut=0.05)
        if(hotf==4){ caselist <- caselist[caselist[,"ExACfreq"] & caselist[,"VCFPASS"] & caselist[,"noneSegmentalDup"] & caselist[,"meta-SVM_PP2"] & caselist[,"alleleFre"], ]; }else{ caselist <- caselist[caselist[,"filtered"], ];}
        
        caselist
}

contlistFams <- function(){
        
        contlist <- onlyFilterContlist(Ecut=0.01,hotf=4,swi=1)
        contlist2 <- onlyFilterContlist(Ecut=0.01,hotf=4,swi=2)
        save(contlist,file="/home/local/ARCS/qh2159/breast_cancer/variants/trios/contlist")
        save(contlist2,file="/home/local/ARCS/qh2159/breast_cancer/variants/trios/contlist2")

}

onlyFilterContlist <- function(Ecut=0.01,hotf=4,swi=1){
        source("misc.R")
        source("sourcefiles.R")
        
        if(hotf==1){ hotspotfile <- hotHMM; ## HongjianPred hotspots file
        }else if(hotf==2){ hotspotfile <- hotCOSMIC; ## COSMIC hotspots file
        }else if(hotf==3){ hotspotfile="";
        }else if(hotf==4){ hotspotfile <- hotHMM;}
        ## switch to the right files for burden test
        alleleFrefile <- alleleFrefiles[swi] 
        contlistf <- contlistfs[swi]
        mis <- "nonsynonymousSNV"
        load(contlistf)
        contlist <- onelist
        rm(onelist)

        contlist <- variant_filtering(contlist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf=hotspotfile,alleleFrefile,popcut=0.05)
        if(hotf==4){ contlist <- contlist[contlist[,"ExACfreq"] & contlist[,"VCFPASS"] & contlist[,"noneSegmentalDup"] & contlist[,"meta-SVM_PP2"] & contlist[,"alleleFre"], ]; }else{ contlist <- contlist[contlist[,"filtered"], ];}
        
        contlist
}

writePhenoFams <- function(flag=2){
        source("InheritedModels.R")
        # Family ID, H,J, individual ID, subject ID, father ID, mother ID, case or not
        L2CasesFams <- unlist(read.table("/home/local/ARCS/qh2159/breast_cancer/variants/families/FamiliesL2Cases.txt"))
        aa <- read.delim("/home/local/ARCS/qh2159/breast_cancer/variants/families/Prioritized43families.txt")[,1]
        Lars <- unique(read.delim("/home/local/ARCS/qh2159/breast_cancer/variants/families/LargeFamilyPhenotype.txt")[,1])
        fams <- unique(c(L2CasesFams,aa,Lars))
        
        pheno <- phenoin()
        pop <- paste(pheno[,4],pheno[,5],sep="")
        pop[pop=="JH" | pop=="HJ"] <- "J"
        subs <- pheno[,1] %in% fams
        pedis <- read.csv("/home/local/ARCS/qh2159/breast_cancer/variants/pedigree/ALL_pedigree.csv")
        
        phenoFams <- cbind(pheno[subs,1],pop[subs],pheno[subs,c("INDIVID","Subject_ID","BreastCancer","Sex","UNIage")])
        phenoFams <- cbind(phenoFams,pedis[match(phenoFams[,3],pedis[,2]),3:4])
        colnames(phenoFams) <- c("FamilyID","Ethnic","IndividualID","SubjectID","Status","Sex","UNIage","FatherID","MotherID")
        phenoFams[,"FatherID"] <- pheno[match(phenoFams[,"FatherID"],pheno[,2]),3]
        phenoFams[,"MotherID"] <- pheno[match(phenoFams[,"MotherID"],pheno[,2]),3]
        
        if(flag==2) qwt(phenoFams,file="/home/local/ARCS/qh2159/breast_cancer/variants/families/PhenotypeInfo.txt",flag=2)
        phenoFams
}

writeGenoFams <- function(){
        phenoFams <- writePhenoFams(1)
        ## excluding yonger than 50 and no cancer
        phenoFams <- phenoFams[!(phenoFams[,"Status"]=="No" & phenoFams[,"UNIage"] <= 50),  ]
        sams <- unique(phenoFams[,"SubjectID"])
        
        LOF <- c(".","none","stopgain","stoploss")
        DMIS <- "nonsynonymousSNV"
        indLOF <- c("frameshiftdeletion","frameshiftinsertion")
        indMIS <- c("nonframeshiftdeletion","nonframeshiftinsertion")
        syn <- "synonymousSNV"
        vTyp <- list(LOF,DMIS,indLOF,indMIS,syn)
        vstr <- c("LOF","MIS","indelLOF","indelMIS","SYN")
        
#         load("/home/local/ARCS/qh2159/breast_cancer/variants/trios/caselist")
#         for(i in 1:length(vTyp)){
#                 tmp <- GenoOneCaseL(caselist,sams,paste("AJ",vstr[i],sep=""),vTyp[[i]])       
#         }
#         
#         load("/home/local/ARCS/qh2159/breast_cancer/variants/trios/caselist2")
#         for(i in 1:length(vTyp)){
#                 tmp <- GenoOneCaseL(caselist2,sams,paste("HI",vstr[i],sep=""),vTyp[[i]])       
#         }
        
        caselistA <- caselistFams()
        for(i in 1:length(vTyp)){
                tmp <- GenoOneCaseL(caselistA,sams,paste("All",vstr[i],sep=""),vTyp[[i]])       
        }
}

GenoOneCaseL <- function(caselist,sams,wstr,Vtype){
        
        caselist <- caselist[caselist[,"Subject_ID"] %in% sams, ]
        caselist <- caselist[caselist[,"VariantClass"] %in% Vtype, ]
        
        vars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
        uniV <- unique(vars)
        genes <- caselist[match(uniV,vars),"Gene"]
        uniG <- unique(genes)
        uniID <- unique(caselist[,"Subject_ID"])
        
        genotypeV <- matrix(0,length(uniV),length(uniID)+2)
        genotypeV[,1] <- genes
        genotypeV[,2] <- uniV
        rownames(genotypeV) <- uniV
        colnames(genotypeV) <- c("Gene","SNP",uniID)
        for(i in 3:dim(genotypeV)[2]){
                oner <- caselist[caselist[,"Subject_ID"] == uniID[i-2], ]
                oneV <- paste(oner[,1],oner[,2],oner[,4],oner[,5],sep="_")
                subs1 <- grepl("0/0",oner[,"GT"]) | grepl("\\./\\.",oner[,"GT"])
                
                genotypeV[oneV[subs1],i] <- 0
                subs2 <- grepl("0/1",oner[,"GT"])
                genotypeV[oneV[subs2],i] <- 1
                genotypeV[oneV[!subs1 & !subs2],i] <- 2
        }
        if(wstr!="") qwt(genotypeV,file=paste("/home/local/ARCS/qh2159/breast_cancer/variants/families/Genotype",wstr,"Info.txt",sep=""),flag=2)
        genotypeV
}

RunFSKAT_QQplot <- function(){
        library(kinship)
        library(CompQuadForm)
        source("InheritedModels.R")
        source("/home/local/ARCS/qh2159/breast_cancer/variants/families/F_SKAT/F-SKAT/glmmPQL.s")
        source("/home/local/ARCS/qh2159/breast_cancer/variants/families/F_SKAT/F-SKAT/FSKAT_HQ.R")
        
        # Subject IDs are character
        vstr <- c("LOF","MIS","indelLOF","indelMIS","SYN")
        for(i in 1:length(vstr)){
                #oneRunFSKAT(paste("AJ",vstr[i],sep=""))
                #oneRunFSKAT(paste("HI",vstr[i],sep=""))
                oneRunFSKAT(paste("All",vstr[i],sep=""))
        }
        
}

oneRunFSKAT <- function(wstr){
        source("misc.R")
        
        y <- read.delim("/home/local/ARCS/qh2159/breast_cancer/variants/families/PhenotypeInfo.txt")
        gene <- read.delim(paste("/home/local/ARCS/qh2159/breast_cancer/variants/families/Genotype",wstr,"Info.txt",sep=""),check.names = FALSE)
        y <- y[y[,"SubjectID"] %in% colnames(gene), ]
        subs <- match(y[,"SubjectID"],colnames(gene))
        lab <- y[,"Status"]
        lab[lab=="Yes"] <- 1
        lab[lab=="No"] <- 0
        covs <- y[,c("Sex","UNIage")]
        covs[covs[,1]=="Male",1] <- 0
        covs[covs[,1]=="Female",1] <- 1
        covs[,1] <- as.numeric(covs[,1])
        covs[,2] <- as.numeric(covs[,2])
        gene <- gene[rowSums(gene[,3:dim(gene)[2]])>0, ]
        
        pvalue1 <- FSKAT(phenotype=as.numeric(lab), genotypes=as.data.frame(gene[,c(1,2,subs)]), id=as.character(y[,"SubjectID"]), fa=as.character(y[,"FatherID"]), mo=as.character(y[,"MotherID"]), family="binomial", covariates=NULL, weights=NULL)
        pva1 <- as.numeric(pvalue1[,2])
        pva1[pva1>=1] <- 1
        pdf(file=paste("/home/local/ARCS/qh2159/breast_cancer/variants/families/F_SKAT/F-SKAT/",wstr,"QQ.pdf",sep=""),width=12,height=10)
        par(mai=c(2,1,1,1))
        PQQ(pva1,main="",labels=pvalue1[,1],nlab=20)
        dev.off()
        qwt(pvalue1,file=paste("/home/local/ARCS/qh2159/breast_cancer/variants/families/F_SKAT/F-SKAT/",wstr,"pVas.txt",sep=""))
        
        pvalue1 <- FSKAT(phenotype=as.numeric(lab), genotypes=as.data.frame(gene[,c(1,2,subs)]), id=as.character(y[,"SubjectID"]), fa=as.character(y[,"FatherID"]), mo=as.character(y[,"MotherID"]), family="binomial", covariates=data.frame("UNIage"=covs[,2]), weights=NULL)
        pva1 <- as.numeric(pvalue1[,2])
        pva1[pva1>=1] <- 1
        pdf(file=paste("/home/local/ARCS/qh2159/breast_cancer/variants/families/F_SKAT/F-SKAT/",wstr,"covsQQ.pdf",sep=""),width=12,height=10)
        par(mai=c(2,1,1,1))
        PQQ(pva1,main="",labels=pvalue1[,1],nlab=20)
        dev.off()
        qwt(pvalue1,file=paste("/home/local/ARCS/qh2159/breast_cancer/variants/families/F_SKAT/F-SKAT/",wstr,"covspVas.txt",sep=""))
}

Statistic_Sigene <- function(){

        genes <- c("NEIL1","ADRA2B","ANK3","FMN2","POU4F2","PDIA2","C20orf141","DPY19L4")
        tests <- c("indelLOF","indelMIS","indelMIS","indelMIS","indelMIS","indelMIS","LOF","LOF")
        
        frecols <- c("N.index.AJ","N.pseudoCont.AJ","N.case.AJ","N.non_case.AJ","N.cont.AJ","N.index.HI","N.pseudoCont.HI","N.case.HI","N.non_case.HI","N.cont.HI")
        samcols <- c("index.AJ","pseudoCont.AJ","case.AJ","non_case.AJ","cont.AJ","index.HI","pseudoCont.HI","case.HI","non_case.HI","cont.HI")
        oddcols <- c("Odds_AJ_pseudo","p_AJ_pseudo","Odds_AJ_cont","p_AJ_cont","Odds_HI_pseudo","p_HI_pseudo","Odds_HI_cont","p_HI_cont")
        cols <- c(frecols,oddcols,samcols)
        HIvars <- read.delim("/home/local/ARCS/qh2159/breast_cancer/Panel/resultf/HISP_variant_level_burden_Pseducont.txt")
        AJvars <- read.delim("/home/local/ARCS/qh2159/breast_cancer/Panel/resultf/AJ_variant_level_burden_Pseducont.txt")
        ncol <- length(cols)
        
        ## variant number in populations and ## IGV plot for each variant
        phenoFams <- writePhenoFams(1)
        phenoFams <- phenoFams[!(phenoFams[,"Status"]=="No" & phenoFams[,"UNIage"] <= 50),  ]
        BRp <- phenoFams[phenoFams[,"Status"]=="Yes", "SubjectID"]
        BRn <- phenoFams[phenoFams[,"Status"]=="No", "SubjectID"]
        
        PopSta <- c()
        IGVS <- c()
        for(i in 1:length(genes)){
                geno <- read.delim(paste("/home/local/ARCS/qh2159/breast_cancer/variants/families/Genotype",paste("All",tests[i],sep=""),"Info.txt",sep=""),check.names = FALSE)     
                n <- dim(geno)[2]
                sams <- colnames(geno)[3:n]
                oner <- geno[geno[,1]==genes[i], ]
                oneSta <- c()
                oneIGV <- c()
                for(j in 1:dim(oner)[1]){
                        carriers <- sams[oner[j,3:n] > 0]
                        non_cars <- sams[oner[j,3:n] == 0]
                        BRpVarp <- intersect(carriers,BRp)
                        BRpVarn <- intersect(non_cars,BRp)
                        BRnVarp <- intersect(carriers,BRn)
                        BRnVarn <- intersect(non_cars,BRn)
                        popN <- c(HIvars[HIvars[,"Variant"]==oner[j,2],cols], AJvars[AJvars[,"Variant"]==oner[j,2],cols])
                        if(all(HIvars[,"Variant"]!=oner[j,2])) popN[1:ncol] <- ""
                        if(all(AJvars[,"Variant"]!=oner[j,2])) popN[(ncol+1):(2*ncol)] <- ""
                        
                        tmp <- c(oner[j,1:2],length(BRpVarp),length(BRpVarn),length(BRnVarp),length(BRnVarn),paste(BRpVarp,sep="",collapse = "_"), paste(BRnVarp,sep="",collapse = "_"),popN)
                        oneSta <- rbind(oneSta,tmp)
                        
                        tmp <- unlist(strsplit(oner[j,2],"_"))[1:2]
                        oneIGV <- rbind(oneIGV,cbind(tmp[1],tmp[2],carriers))
                }
                
                PopSta <- rbind(PopSta,oneSta)
                IGVS <- rbind(IGVS,oneIGV)
        }
        colnames(PopSta) <- c("Gene","SNP","#BR+Var+","#BR+Var-","#BR-Var+","#BR-Var-","BR+Var+","BR-Var+",cols,cols)
        qwt(PopSta,file="/home/local/ARCS/qh2159/breast_cancer/variants/families/F_SKAT/SignGenesPopSta.txt",flag=2)
        qwt(IGVS,file="/home/local/ARCS/qh2159/breast_cancer/variants/families/F_SKAT/SignGenesIGV.txt")
}
