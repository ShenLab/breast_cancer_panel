rvTDTtest <- function(Vtype,wstr){
        source("/home/local/ARCS/qh2159/breast_cancer/Panel/breast_cancer_panel/InheritedModels.R")
        source("/home/local/ARCS/qh2159/breast_cancer/Panel/breast_cancer_panel/sourcefiles.R")
        pheno <- phenoin()
        pedis <- read.csv("/home/local/ARCS/qh2159/breast_cancer/variants/pedigree/ALL_pedigree.csv")
        subs <- (pedis[,2] %in% pheno[,2]) & (pedis[,3] %in% pheno[,2]) & (pedis[,4] %in% pheno[,2])
        trioped <- pedis[subs, ]
        trioped <- trioped[trioped[,2] %in% pheno[pheno[,"BreastCancer"]=="Yes",2], ]
        #trioped <- trioped[!(trioped[,3] %in% trioped[,2]) & !(trioped[,4] %in% trioped[,2]), ]
       
        kids <- trioped[,2]
        mots <- trioped[,4]
        fats <- trioped[,3]
        kids <- pheno[match(kids,pheno[,2]), 3]
        mots <- pheno[match(mots,pheno[,2]), 3]
        fats <- pheno[match(fats,pheno[,2]), 3]
        sams <- unique(c(kids,mots,fats))
        
        #### get trio variant information
        load("caselist")
        load("caselist2")
        caselist <- rbind(caselist,caselist2)
        caselist <- caselist[caselist[,"Subject_ID"] %in% sams, ]
        caselist <- caselist[caselist[,"VariantClass"] %in% Vtype, ]
        
        vars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
        uniV <- unique(vars)
        genes <- caselist[match(uniV,vars),"Gene"]
        uniG <- unique(genes)
        
        trioV <- matrix(0,3*length(kids),length(uniV))
        colnames(trioV) <- uniV
        rownames(trioV) <- c(kids,mots,fats)
        rows <- c(kids,mots,fats)
        for(i in 1:dim(trioV)[1]){
                oner <- caselist[caselist[,"Subject_ID"] == rows[i], ]
                oneV <- paste(oner[,1],oner[,2],oner[,4],oner[,5],sep="_")
                subs1 <- grepl("0/0",oner[,"GT"]) | grepl("\\./\\.",oner[,"GT"])
                trioV[i,oneV[subs1]] <- 0
                subs2 <- grepl("0/1",oner[,"GT"])
                trioV[i,oneV[subs2]] <- 1
                trioV[i,oneV[!subs1 & !subs2]] <- 2
        }
        
        #### get control variant information
        load("contlist")
        load("contlist2")
        contlist <- rbind(contlist,contlist2)
        contV <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
        subs <- contV %in% uniV
        contlist <- contlist[subs, ]
        contV <- contV[subs]
        
        n.cont <- length(unique(contlist[,"Subject_ID"]))
        n.dep <- 50
        a1 <- contV[grepl("0/1",contlist[,"GT"])]
        ta1 <- table(a1)
        na1 <- intersect(names(ta1),uniV)
        a2 <- contV[!grepl("0/1",contlist[,"GT"])]
        ta2 <- table(a2)
        na2 <- intersect(names(ta2),uniV)
        deps <- sapply(1:dim(contlist)[1], function(kk) sum(as.numeric(unlist(strsplit(contlist[kk,"AD"],",")))))
        
        cons <- matrix(0,length(uniV),4)
        rownames(cons) <- uniV
        cons[,4] <- n.dep
        cons[na2,1] <- ta2[na2]
        cons[na1,2] <- ta1[na1]
        cons[,3] <- n.cont - cons[,2] - cons[,1]
        na3 <- intersect(uniV,contV)
        tmp <- sapply(1:length(na3), function(kk) mean(deps[contV==na3[kk]]))
        cons[na3,4] <- tmp 
        
        #### running rvTDT
        library(rvTDT)
        #aa <- rvTDT(trioV,cons,1) ## time consuming
        ### single gene rvTDT test
        rvTr <- sapply(uniG, function(gene){ if(sum(genes==gene)>=2) c(gene,unlist(rvTDT(trioV[ ,genes==gene],cons[genes==gene, ],1))) else c(gene,rep(1,9)); })
        rvST <- t(rvTr)
        rvST <- rvST[order(as.numeric(rvST[,7])),]
        colnames(rvST) <- c("Gene","nSNP","nfamily","nSNPc","p_lc_1","p_lc_maf","p_lc_pc","p_k_1","p_k_maf","p_k_pc")
        qwt(rvST,file=paste("rvTDTtestGene",wstr,".txt",sep=""),flag=2)
    
}

AJ_rvTDT <- function(Vtype,wstr){
        ### AJ trios only with controls
        source("/home/local/ARCS/qh2159/breast_cancer/Panel/breast_cancer_panel/InheritedModels.R")
        source("/home/local/ARCS/qh2159/breast_cancer/Panel/breast_cancer_panel/sourcefiles.R")
        pheno <- phenoin()
        AJs <- unlist(read.table(AJBRfile))
        pheno <- pheno[pheno[,3] %in% AJs, ]
        pedis <- read.csv("/home/local/ARCS/qh2159/breast_cancer/variants/pedigree/ALL_pedigree.csv")
        subs <- (pedis[,2] %in% pheno[,2]) & (pedis[,3] %in% pheno[,2]) & (pedis[,4] %in% pheno[,2])
        trioped <- pedis[subs, ]
        ### unique trios and probands are breast cancer cases
        trioped <- trioped[trioped[,2] %in% pheno[pheno[,"BreastCancer"]=="Yes",2], ]
        #trioped <- trioped[!(trioped[,3] %in% trioped[,2]) & !(trioped[,4] %in% trioped[,2]), ]
        #qwt(trioped,file="AJcaseTrios.txt")
        kids <- trioped[,2]
        mots <- trioped[,4]
        fats <- trioped[,3]
        kids <- pheno[match(kids,pheno[,2]), 3]
        mots <- pheno[match(mots,pheno[,2]), 3]
        fats <- pheno[match(fats,pheno[,2]), 3]
        sams <- unique(c(kids,mots,fats))
        
        #### get trio variant information
        load("caselist")
        caselist <- caselist[caselist[,"Subject_ID"] %in% sams, ]
        caselist <- caselist[caselist[,"VariantClass"] %in% Vtype, ]
        
        vars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
        uniV <- unique(vars)
        genes <- caselist[match(uniV,vars),"Gene"]
        uniG <- unique(genes)
        
        trioV <- matrix(0,3*length(kids),length(uniV))
        colnames(trioV) <- uniV
        rownames(trioV) <- c(kids,mots,fats)
        rows <- c(kids,mots,fats)
        for(i in 1:dim(trioV)[1]){
                oner <- caselist[caselist[,"Subject_ID"] == rows[i], ]
                oneV <- paste(oner[,1],oner[,2],oner[,4],oner[,5],sep="_")
                subs1 <- grepl("0/0",oner[,"GT"]) | grepl("\\./\\.",oner[,"GT"])
                trioV[i,oneV[subs1]] <- 0
                subs2 <- grepl("0/1",oner[,"GT"])
                trioV[i,oneV[subs2]] <- 1
                trioV[i,oneV[!subs1 & !subs2]] <- 2
        }
        
        #### get control variant information
        load("contlist")
        contV <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
        subs <- contV %in% uniV
        contlist <- contlist[subs, ]
        contV <- contV[subs]
        
        n.cont <- length(unique(contlist[,"Subject_ID"]))
        n.dep <- 50
        a1 <- contV[grepl("0/1",contlist[,"GT"])]
        ta1 <- table(a1)
        na1 <- intersect(names(ta1),uniV)
        a2 <- contV[!grepl("0/1",contlist[,"GT"])]
        ta2 <- table(a2)
        na2 <- intersect(names(ta2),uniV)
        deps <- sapply(1:dim(contlist)[1], function(kk) sum(as.numeric(unlist(strsplit(contlist[kk,"AD"],",")))))
        
        cons <- matrix(0,length(uniV),4)
        rownames(cons) <- uniV
        cons[,4] <- n.dep
        cons[na2,1] <- ta2[na2]
        cons[na1,2] <- ta1[na1]
        cons[,3] <- n.cont - cons[,2] - cons[,1]
        na3 <- intersect(uniV,contV)
        tmp <- sapply(1:length(na3), function(kk) mean(deps[contV==na3[kk]]))
        cons[na3,4] <- tmp 

        #### -----------running rvTDT
        library(rvTDT)
        #aa <- rvTDT(trioV,cons,1)
        ### single gene rvTDT test
        rvTr <- sapply(uniG, function(gene){ if(sum(genes==gene)>=2) c(gene,unlist(rvTDT(trioV[ ,genes==gene],cons[genes==gene, ],1))) else c(gene,rep(1,9)); })
        #subs <- as.numeric(rvTr[7,]) < 0.05
        rvST <- t(rvTr)
        rvST <- rvST[order(as.numeric(rvST[,7])),]
        colnames(rvST) <- c("Gene","nSNP","nfamily","nSNPc","p_lc_1","p_lc_maf","p_lc_pc","p_k_1","p_k_maf","p_k_pc")
        qwt(rvST,file=paste("rvTDTtestGene_AJ",wstr,".txt",sep=""),flag=2)

}

rundemo <- function(){
        source("rvTDTtest.R")
        SINs <- c(".","none","nonsynonymousSNV","stopgain","stoploss")
        LOF <- c(".","none","stopgain","stoploss")
        DMIS <- "nonsynonymousSNV"
        indLOF <- c("frameshiftdeletion","frameshiftinsertion")
        indMIS <- c("nonframeshiftdeletion","nonframeshiftinsertion")
        syn <- "synonymousSNV"
        vTyp <- list(SINs,LOF,DMIS,indLOF,indMIS,syn)
        vstr <- c("DAM","LOF","MIS","indelLOF","indelMIS","SYN")

        for(i in 1:length(vTyp)){
                #rvTDTtest(vTyp[[i]],vstr[i])
                AJ_rvTDT(vTyp[[i]],vstr[i])
        }
        
        
#         nfamily: number of total families in computation
#         
#         nsnptot: the total number of snps that in the input files
#         
#         nsnpcompute: the number of snps that pass the QC
#         
#         p_lc_1: p value of unweighted linear combinated TDT
#         
#         p_lc_maf: p value of linear combinated TDT weighted by MAF
#         (dbeta(1,25,maf))
#         
#         p_lc_pc: p value of linear combinated TDT weighted by population
#         controls
#         
#         p_k_1: p value of unweighted kernel TDT
#         
#         p_k_maf: p value of kernel TDT weighted by MAF (dbeta(1,25,maf))
#         
#         p_k_pc: p value of kernel TDT weighted by population controls
        
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
