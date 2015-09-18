pre <- function(){
    ##  Sinai target gene list
    target <- read.delim("0627081_Covered_3col_b37.withGene.tab",sep="\t")
    genes <- setdiff(unique(target[,4]),"")
    qwt(genes,file="Sinai_targets.txt",sep="")
    
    # CNV information for TCGA samples
    cnv1f <- "../somatic_mutation/CNV_cbio/Breast Invasive Carcinoma TCGA Nature_825.txt"
    cnv2f <- "../somatic_mutation/CNV_cbio/Breast Invasive Carcinoma_1104.txt"
    cnv1 <- read.delim(cnv1f,check.names=FALSE)
    cnv2 <- read.delim(cnv2f,check.names=FALSE)
    cnv <- rbind(cnv1,cnv2)
    ## cutoff for CNV fraction
    rc <- 0.5
    samples <- unique(cnv[cnv[,2] <= rc,1])
    
    ## somatic mutations in TCGA excluding the CNV
    filename <- "../somatic_mutation/Somatic_Mutations_TCGA/genome.wustl.edu__IlluminaGA_curated_DNA_sequencing_level2.maf"
    tcgaV <- read.delim(filename,check.names=FALSE)
    tcgaSample <- sapply(1:dim(tcgaV)[1], function(i) substr(tcgaV[i,"Tumor_Sample_Barcode"],1,15))
    subs <- tcgaSample %in% samples
    tcgaV <- tcgaV[subs,]
    muta <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Site")
    tcgaV <- tcgaV[tcgaV[,"Variant_Classification"] %in% muta,]
    lof <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonsense_Mutation","Nonstop_Mutation","Splice_Site")
    mis <- "Missense_Mutation"
    
    genes <- unique(tcgaV[,1])
    geneT <- matrix(0,length(genes),3,dimnames=list(genes,c("Gene","dn.LoF","dn.mis")))
    geneT[,1] <-  genes
    tmp <- tcgaV[tcgaV[,"Variant_Classification"] %in% lof,]
    geneT[match(names(table(tmp[,1])),geneT[,1]),2] <- table(tmp[,1])
    tmp <- tcgaV[tcgaV[,"Variant_Classification"] %in% mis,]
    geneT[match(names(table(tmp[,1])),geneT[,1]),3] <- table(tmp[,1])
    
    ## add the mutation rate
    rate <- read.csv("GeneMutaRate.csv")
    ig <- intersect(rate[,1],geneT[,1])
    rate <- cbind(rate,0,0)
    rate[match(ig,rate[,1]),6] <- geneT[match(ig,geneT[,1]),2]
    rate[match(ig,rate[,1]),7] <- geneT[match(ig,geneT[,1]),3]
    colnames(rate)[6:7] <- c("dn.LoF","dn.mis")
    write.csv(rate,file="TCGAmut.csv",row.names=FALSE)
    ### poisson test
    n <- length(unique(tcgaV[,"Tumor_Sample_Barcode"])) ##803
    ratep <- Poisson_test_hq("TCGAmut.csv",n)
    write.csv(ratep,file="TCGAmutp.csv",row.names=FALSE)
    
    ## /ifs/scratch/c2b2/ys_lab/yshen/WENDY/BreastCancer/Regeneron/Filtering_for_Qiang
    control <- read.csv("controls.csv")
    pheno <- read.csv("WES BCFR phenotypic data.csv")
    
    ## excluded outlier samples and undetermined and mismatched sex samples
    outliers <- read.delim("Potential_Outliers.tsv")[,1]
    sexcheck <- read.csv("CUMC_Regeneron.sexcheck.csv")
    checked <- sexcheck[sexcheck[,"Match"]==TRUE,1]
    chsamples <- setdiff(checked,outliers)
    
    pheno <- pheno[pheno[,"Subject_ID"] %in% chsamples,]
    control <- control[control[,"Subject_ID"] %in% chsamples,]
    
    cases <- pheno[pheno[,"BreastCancer"]=="Yes" & pheno[,"Sex"]=="Female",]
    subject_ID <- gsub(" ","",cases[,"Subject_ID"])
    control_ID <- gsub(" ","",control[,"Subject_ID"])
    
    casef <- paste(subject_ID,".tsv",sep="")
    contf <- paste(control_ID,".tsv",sep="")
    path="/ifs/scratch/c2b2/ys_lab/yshen/WENDY/BreastCancer/Regeneron/Filtering_for_Qiang/"
    files <-  list.files(path=path,pattern=".tsv$")
    ##"222966.tsv" %in%  files  # TRUE
    ##"222357.tsv" %in%  files # TRUE
    casef <- gsub("222357,222966.tsv","222357.tsv",casef)
    ##"220897.tsv" %in% files #TRUE
    ##"220904.tsv" %in% files # FALSE
    casef <- gsub("220897,220904.tsv","220897.tsv",casef)
    casef <- intersect(casef,files)
    contf <- intersect(contf,files)
    ##casef <- setdiff(casef,files)
    
    caselist <- c()
    for(i in 1:length(casef)){
        tmp <- paste(path,casef[i],sep="")
        oner <- read.delim(tmp)
        oner <- cbind(oner,gsub(".tsv","",casef[i]))
        colnames(oner)[c(23,24,29,30)] <- c("GT","AD","Subject_INFO","SubID")
        caselist <- rbind(caselist,oner)
    }
    caselist[caselist[,"AlleleFrequency.ExAC"]==".","AlleleFrequency.ExAC"] <- 0
    ## max(as.numeric(caselist[as.numeric(caselist[,"AlleleFrequency.ExAC"]) >0,"AlleleFrequency.ExAC"]))
    ##ExAC < 1% 
    save(caselist,file="caselist")
    
    length(unique(caselist[,"Gene"]))
    length(unique(caselist[,"SubID"]))
    
    contlist <- c()
    for(i in 1:length(contf)){
        tmp <- paste(path,contf[i],sep="")
        oner <- read.delim(tmp)
        oner <- cbind(oner,gsub(".tsv","",contf[i]))
        colnames(oner)[c(23,24,29,30)] <- c("GT","AD","Subject_INFO","SubID")
        contlist <- rbind(contlist,oner)
    }
    contlist[contlist[,"AlleleFrequency.ExAC"]==".","AlleleFrequency.ExAC"] <- 0
    save(contlist,file="contlist")
    
    ### filtered: control
    #a1 <- paste(caselist[,1],caselist[,2],sep="_")
    #a2 <- paste(contlist[,1],contlist[,2],sep="_")
    #subs <- which(!(a1 %in% a2))
    #caselist <- caselist[subs,]
    
    ### ExAC 0.1%
    subs <- as.numeric(caselist[,"AlleleFrequency.ExAC"]) < 0.01
    caselist <- caselist[subs,]
    
    ### excluding unknown and none mutations
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss")                 
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    caselist <- caselist[caselist[,"VariantClass"] %in% c(lof,mis),]
    
    ### split damaging missense and misense: MetaSVM
    cosm <- read.delim("CosmicMutantExportCensus.tsv",check.names=FALSE)
    cosms <- cosm[,"Mutation genome position"]
    cosms <- sapply(1:length(cosms),function(i) unlist(strsplit(cosms[i],"-"))[1])
    idcos=FALSE
    
    casems <- paste(caselist[,1],caselist[,2],sep=":")
    subs <- rep("MIS",dim(caselist)[1])
    if(idcos){
        subs[caselist[,"MetaSVM"]=="D" & casems %in% cosms] <- "dMIS"
    }else{
        subs[caselist[,"MetaSVM"]=="D"] <- "dMIS"
    }
    subs[caselist[,"VariantClass"] %in% lof] <- "LOF"
    caselist <- cbind(caselist,subs)
    
    n.case <- length(unique(caselist[,"SubID"]))
    ### poisson test
    genes <- unique(caselist[,"Gene"])
    geneT <- matrix(0,length(genes),4,dimnames=list(genes,c("Gene","dn.LoF","dn.mis3","dn.mis")))
    geneT[,1] <-  genes
    tmp <- caselist[caselist[,"VariantClass"] %in% lof,]
    geneT[match(names(table(tmp[,"Gene"])),geneT[,1]),2] <- table(tmp[,"Gene"])
    tmp <- caselist[caselist[,"VariantClass"] %in% mis & subs=="dMIS",]
    geneT[match(names(table(tmp[,"Gene"])),geneT[,1]),3] <- table(tmp[,"Gene"])
    tmp <- caselist[caselist[,"VariantClass"] %in% mis & subs=="MIS",]
    geneT[match(names(table(tmp[,"Gene"])),geneT[,1]),4] <- table(tmp[,"Gene"])
    write.csv(geneT,file="case_num.csv",row.names=FALSE)
    
    
    idmis=FALSE
    if(idmis){
        ratep <- add_rate(geneT[,1:3],"Remuta.csv",n.case,idmis=TRUE)
        ig <- intersect(ratep[,1],geneT[,1])
        ratep[match(ig,ratep[,1]),"dn.mis"] <- geneT[match(ig,geneT[,1]),4] 
        write.csv(ratep,file="Remutall.csv",row.names=FALSE)
    }else{
        geneT[,3] <- as.numeric(geneT[,3]) + as.numeric(geneT[,4])
        ratep <- add_rate(geneT[,1:3],"Remutamis.csv",n.case,idmis)
    }
    
    ###========================================================
    ### control analysis
    ### ExAC 0.1%
    subs <- as.numeric(contlist[,"AlleleFrequency.ExAC"]) < 0.01
    contlist <- contlist[subs,]
    
    ### excluding unknown and none mutations
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss")                 
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    contlist <- contlist[contlist[,"VariantClass"] %in% c(lof,mis),]
    
    ### split damaging missense and misense: MetaSVM
    conms <- paste(contlist[,1],contlist[,2],sep=":")
    
    subs <- rep("MIS",dim(contlist)[1])
    if(idcos){
        subs[contlist[,"MetaSVM"]=="D" & conms %in% cosms] <- "dMIS"
    }else{
        subs[contlist[,"MetaSVM"]=="D"] <- "dMIS"
    }
    subs[contlist[,"VariantClass"] %in% lof] <- "LOF"
    contlist <- cbind(contlist,subs)
    
    n.case <- length(unique(contlist[,"SubID"]))
    ### poisson test
    genes <- unique(contlist[,"Gene"])
    geneT <- matrix(0,length(genes),4,dimnames=list(genes,c("Gene","dn.LoF","dn.mis3","dn.mis")))
    geneT[,1] <-  genes
    tmp <- contlist[contlist[,"VariantClass"] %in% lof,]
    geneT[match(names(table(tmp[,"Gene"])),geneT[,1]),2] <- table(tmp[,"Gene"])
    tmp <- contlist[contlist[,"VariantClass"] %in% mis & subs=="dMIS",]
    geneT[match(names(table(tmp[,"Gene"])),geneT[,1]),3] <- table(tmp[,"Gene"])
    tmp <- contlist[contlist[,"VariantClass"] %in% mis & subs=="MIS",]
    geneT[match(names(table(tmp[,"Gene"])),geneT[,1]),4] <- table(tmp[,"Gene"])
    write.csv(geneT,file="control_num.csv",row.names=FALSE)
    
    ### control mutation fraction estimate
    dlof <- sum(contlist[,"VariantClass"] %in% lof)
    dmis <- sum(contlist[,"VariantClass"] %in% mis & contlist[,"MetaSVM"]=="D")
    mu.frac <- 1:2
    mu.frac[1] <- dlof/dim(contlist)[1]
    mu.frac[2] <- dmis/dim(contlist)[1]
}

add_genelist <- function(){
    tadar <- read.csv("TADAdenovo_breastcase.csv")
    tadar[,"Sinai_target"] <- FALSE
    tadar[,"COSMIC"] <- FALSE
    
    sit <- unlist(read.table("Sinai_targets.txt"))
    cot <- unlist(read.table("COSMIC_censusgenes.txt"))
    tadar[match(intersect(sit,tadar[,1]),tadar[,1]),"Sinai_target"] <- TRUE
    tadar[match(intersect(cot,tadar[,1]),tadar[,1]),"COSMIC"] <- TRUE
    
    tcga <- read.csv("TCGAmutp.csv")
    ig <- intersect(tadar[,1],tcga[,1])
    tadar[,paste("TCGA",c("dn.LoF","dn.mis","p_LOF","p_mis","p_both","min_p"),sep="_")] <- "."
    tadar[match(ig,tadar[,1]),paste("TCGA",c("dn.LoF","dn.mis","p_LOF","p_mis","p_both","min_p"),sep="_")] <- tcga[match(ig,tcga[,1]),c("dn.LoF","dn.mis","p_LOF","p_mis","p_both","min_p")]
    
    write.csv(tadar,file="TADAresult_anno.csv",row.names=FALSE)
}

add_rate <- function(geneT,filename,n.case,idmis){
    ## add the mutation rate
    rate <- read.csv("GeneMutaRate.csv")
    ig <- intersect(rate[,1],geneT[,1])
    rate <- cbind(rate,0,0)
    rate[match(ig,rate[,1]),6] <- geneT[match(ig,geneT[,1]),2]
    rate[match(ig,rate[,1]),7] <- geneT[match(ig,geneT[,1]),3]
    if(idmis){
        colnames(rate)[6:7] <- c("dn.LoF","dn.mis3")
    }else{ colnames(rate)[6:7] <- c("dn.LoF","dn.mis");}
    write.csv(rate,file=filename,row.names=FALSE)
    ### poisson test
    ratep <- Poisson_test_hq(filename,n.case,idmis)
    write.csv(ratep,file=gsub(".csv","p.csv",filename),row.names=FALSE)
    
    ratep
}

Poisson_test_hq <- function(filename,ntrio,idmis=FALSE){
    tmp <- read.csv(filename)
    pV <- poisson_test(tmp,ntrio,idmis)
    tmp <- cbind(tmp,pV)
    tmp[,"min_p"] <- apply(pV,1,min)
    tmp <- tmp[order(tmp[,"min_p"]),]
 
    tmp
}

poisson_test <- function(data,N,idmis){
    
    n <- dim(data)[1]
    if(idmis){
        nmis="dn.mis3"
        rmis="dmis"
    }else{
        nmis="dn.mis"
        rmis="mis"
    }
    logpM <- sapply(1:n, function(i){
        a1 <- poisson.test(x=data[i,"dn.LoF"], T = N *2* data[i,"LOF"], alternative = "greater", conf.level = 0.95)$p.value
        a2 <- poisson.test(x=data[i,nmis], T = N *2* data[i,rmis], alternative = "greater", conf.level = 0.95)$p.value
        a3 <- poisson.test(x=data[i,"dn.LoF"] +  data[i,nmis], T = N *2* (data[i,"LOF"]+ data[i,rmis]), alternative = "greater", conf.level = 0.95)$p.value 
        c(a1,a2,a3)
    })
    logpM <- t(logpM)
    colnames(logpM) <- c("p_LOF","p_mis","p_both")
    
    logpM
}

twogenelist <- function(){
    glist1 <- unlist(read.table("../genelist/Genelist1.txt"))
    glist2 <- unlist(read.table("../genelist/Genelist2.txt"))
    
    ### case number; control number; burden test; poisson test; TADA test,  
    Tglist1 <- add_table(glist1)
    Tglist2 <- add_table(glist2)
    write.csv(Tglist1,file="Genelist1_analysis.csv",row.names=FALSE)
    write.csv(Tglist2,file="Genelist2_analysis.csv",row.names=FALSE)

    glist <- read.csv("case_num.csv")[,1]
    Tglist <- add_table(glist)
    write.csv(Tglist,file="All_analysis.csv",row.names=FALSE)
    
}

add_table <- function(glist,ncase=538,ncont=72){
    
    caseall <- read.csv("case_num.csv")
    contall <- read.csv("control_num.csv")
    
    glist <- intersect(caseall[,1],glist)
    
    Ts <- caseall[match(glist,caseall[,1]),]
    Ts <- cbind(Ts,contall[match(glist,contall[,1]),2:4])
    for(i in 5:7){
        Ts[is.na(Ts[,i]),i] <- 0
    }
    Ts <- cbind(Ts,1,1,1,1)
    
    colnames(Ts)[2:11] <- c("LoF","dmis","mis",paste(c("LoF","dmis","mis"),"control",sep="_"),"odds.lof","odds.dmis","odds.both","odds.mis")
    
    n <- dim(Ts)[2]
    for(i in 1:dim(Ts)[1]){
        #Ts[i,n-2] <- fisher.test(matrix(c(Ts[i,"LoF"],ncase,Ts[i,"LoF_control"],ncont),2,2))$p.value
        #Ts[i,n-1] <- fisher.test(matrix(c(Ts[i,"mis"],ncase,Ts[i,"mis_control"],ncont),2,2))$p.value
        #Ts[i,n] <- fisher.test(matrix(c(sum(Ts[i,c("LoF","mis3","mis")]) ,ncase,sum(Ts[i,c("LoF_control","mis3_control","mis_control")]),ncont),2,2))$p.value
        
        Ts[i,n-3] <- (Ts[i,"LoF"]/ncase)/(Ts[i,"LoF_control"]/ncont)
        Ts[i,n-2] <- (Ts[i,"dmis"]/ncase)/(Ts[i,"dmis_control"]/ncont)
        Ts[i,n-1] <- ((Ts[i,"dmis"]+Ts[i,"LoF"])/ncase)/((Ts[i,"dmis_control"]+Ts[i,"LoF_control"])/ncont)
        Ts[i,n]   <- (Ts[i,"mis"]/ncase)/(Ts[i,"mis_control"]/ncont)
    }
    
    Ts <- Ts[order(-as.numeric(Ts[,n-1])),]
    
    Ts

}

batch_genes <- function(){
    genes <- c("MLH1","RAD51C","TOX3","BMPR1A")
    load("caselist")
    info <- c()
    for(i in 1:length(genes)){
        tmp <- one_gene(genes[i],caselist)
        info <- rbind(info,tmp)
    }
    write.csv(info,file=paste("Genes","check.csv",sep=""),row.names=FALSE)
}

one_gene <- function(gene,caselist){
    #gene="BMPR1A"
    n <- dim(caselist)[2]
    subs <- caselist[,"Gene"]==gene
    info <- caselist[subs,c(n,1:8,10,13)]
    
    pheno <- read.csv("WES BCFR phenotypic data.csv")
    subs1 <- match(info[,1],pheno[,3])
    info <- cbind(info,pheno[subs1,])
    print(length(unique(pheno[subs1,1])))
    
    info
}

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
caseSKAT <- function(){
    source("pre.R")
    cases <- caseonefam()
    subject_ID <- gsub(" ","",cases[,"Subject_ID"])
    
    casef <- paste(subject_ID,".tsv",sep="")
    path="/ifs/scratch/c2b2/ys_lab/yshen/WENDY/BreastCancer/Regeneron/Filtering_for_Qiang/"
    files <-  list.files(path=path,pattern=".tsv$")
    ##"222966.tsv" %in%  files  # TRUE
    ##"222357.tsv" %in%  files # TRUE
    casef <- gsub("222357,222966.tsv","222357.tsv",casef)
    #casef <- gsub("280002.tsv",files[1328],casef)
    casef <- intersect(casef,files) ## three missed: "222488" "222417" "280002"
    
    caselist <- c()
    for(i in 1:length(casef)){
        tmp <- paste(path,casef[i],sep="")
        oner <- read.delim(tmp)
        oner <- cbind(oner,gsub(".tsv","",casef[i]))
        colnames(oner)[c(23,24,29,30)] <- c("GT","AD","Subject_INFO","SubID")
        caselist <- rbind(caselist,oner)
    }
    save(caselist,file="caselist_7_28")
    
  
    ###==============================================================================
    load("caselist_7_28")
    ### filtered 
    ### excluding unknown and none mutations
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss")                 
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    caselist <- caselist[caselist[,"VariantClass"] %in% c(lof,mis),]
    
    ### split damaging missense and misense: MetaSVM
    subs <- rep("MIS",dim(caselist)[1])
    subs[caselist[,"MetaSVM"]=="D"] <- "dMIS"
    subs[caselist[,"VariantClass"] %in% lof] <- "LOF"
    caselist <- cbind(caselist,subs)
    
    ### ExAC 0.1% :  used exac.nfe and exac.afr (minimum)
    freN <- c("ExAC.nfe.freq","ExAC.afr.freq")
    subs <- sapply(1:dim(caselist)[1], function(i){
        tmp <- unlist(strsplit(caselist[i,"INFO"],";"))
        tmp2 <- unlist(strsplit(tmp,"="))
        min(tmp2[match(freN,tmp2)+1]) < 0.001
    })
    subs[is.na(subs)] <- TRUE
    caselist <- caselist[subs,]
    
    ## caselist <- caselist[caselist[,"VariantClass"] %in% lof | (caselist[,"VariantClass"] %in% mis & caselist[,"MetaSVM"]=="D"),]
    
    ## MAF (minor-allele frequency): 
    freN <- c("ExAC.nfe.freq","ExAC.afr.freq")
    fres <- sapply(1:dim(caselist)[1], function(i){
        tmp <- unlist(strsplit(caselist[i,"INFO"],";"))
        tmp2 <- unlist(strsplit(tmp,"="))
        min(tmp2[match(freN,tmp2)+1])
    })
    fres[is.na(fres)] <- 0
    names(fres) <- paste(caselist[,1],caselist[,2],sep="_")
    vars <- unique(paste(caselist[,1],caselist[,2],sep="_"))
    fres <- fres[match(vars,names(fres))]
    save(fres,file="MAF")
    
    ### genotype matrix in SKAT
    caseids <- unique(caselist[,"SubID"])
    n.case <- length(caseids)
    vars <- unique(paste(caselist[,1],caselist[,2],sep="_"))
    n.var <- length(vars)
    Z <- matrix(0,n.case,n.var,dimnames=list(caseids,vars))
    for(i in 1:n.case){
        tmp <- caselist[caselist[,"SubID"]==caseids[i],]
        svar <- paste(tmp[,1],tmp[,2],sep="_")
        geo <- rep(2,length(svar))
        geo[tmp[,"GT"]=="0/0"] <- 0
        geo[tmp[,"GT"]=="0/1"] <- 1
        Z[i,match(svar,vars)] <- geo
    }
    save(Z,file="genotype")
    
}

## control SKAT analysis
controlSKAT <- function(){
    source("pre.R")
    control <- controlpheno()

    control_ID <- control[,"Subject_ID"]
    contf <- paste(control_ID,".tsv",sep="")
    #path="/ifs/scratch/c2b2/ys_lab/yshen/WENDY/BreastCancer/Regeneron/Filtering_for_Qiang/"
    path="/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_for_Qiang/"
    files <-  list.files(path=path,pattern=".tsv$")
    contf <- intersect(contf,files)
    
    contlist <- c()
    for(i in 1:length(contf)){
        tmp <- paste(path,contf[i],sep="")
        oner <- read.delim(tmp)
        oner <- cbind(oner,gsub(".tsv","",contf[i]))
        colnames(oner)[c(23,24,29,30)] <- c("GT","AD","Subject_INFO","SubID")
        contlist <- rbind(contlist,oner)
    }
    save(contlist,file="contlist_8_6")
    

    load("contlist_7_28")
    ### filtered 
    ### excluding unknown and none mutations
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss")                 
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    contlist <- contlist[contlist[,"VariantClass"] %in% c(lof,mis),]
    
    ### split damaging missense and misense: MetaSVM
    subs <- rep("MIS",dim(contlist)[1])
    subs[contlist[,"MetaSVM"]=="D"] <- "dMIS"
    subs[contlist[,"VariantClass"] %in% lof] <- "LOF"
    contlist <- cbind(contlist,subs)
    
    ### ExAC 0.1% :  used exac.nfe and exac.afr (minimum)
    freN <- c("ExAC.nfe.freq","ExAC.afr.freq")
    subs <- sapply(1:dim(contlist)[1], function(i){
        tmp <- unlist(strsplit(contlist[i,"INFO"],";"))
        tmp2 <- unlist(strsplit(tmp,"="))
        min(tmp2[match(freN,tmp2)+1]) < 0.001
    })
    subs[is.na(subs)] <- TRUE
    contlist <- contlist[subs,]
    
    ## MAF (minor-allele frequency): 
    freN <- c("ExAC.nfe.freq","ExAC.afr.freq")
    fres <- sapply(1:dim(contlist)[1], function(i){
        tmp <- unlist(strsplit(contlist[i,"INFO"],";"))
        tmp2 <- unlist(strsplit(tmp,"="))
        min(tmp2[match(freN,tmp2)+1])
    })
    fres[is.na(fres)] <- 0
    names(fres) <- paste(contlist[,1],contlist[,2],sep="_")
    vars <- unique(paste(contlist[,1],contlist[,2],sep="_"))
    fres <- fres[match(vars,names(fres))]
    fres1 <- fres
    save(fres1,file="MAF_cont")
    
    ### genotype matrix in SKAT
    contids <- unique(contlist[,"SubID"])
    n.cont <- length(contids)
    vars <- unique(paste(contlist[,1],contlist[,2],sep="_"))
    n.var <- length(vars)
    Z1 <- matrix(0,n.cont,n.var,dimnames=list(contids,vars))
    for(i in 1:n.cont){
        tmp <- contlist[contlist[,"SubID"]==contids[i],]
        svar <- paste(tmp[,1],tmp[,2],sep="_")
        geo <- rep(2,length(svar))
        geo[tmp[,"GT"]=="0/0"] <- 0
        geo[tmp[,"GT"]=="0/1"] <- 1
        Z1[i,match(svar,vars)] <- geo
    }
    save(Z1,file="genotype_cont")
    
}

## single variants SKAT analysis
runSKAT <- function(){
    library(SKAT)
    
    pheno <- read.csv("WES BCFR phenotypic data.csv")
    load("genotype")
    load("genotype_cont")
    load("MAF")
    load("MAF_cont")
    
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
    
    X <- matrix(unlist(pheno[match(sams,pheno[,"Subject_ID"]),c("Sex","LiveAge")]),ncol=2)
    X[113,] <- unlist(pheno[512,c("Sex","LiveAge")])
    X[X[,1]=="Female",1] <- 0
    X[X[,1]=="Male",1] <- 1
    X <- as.numeric(X)
    X <- matrix(X,ncol=2)
    y <- rep(0,length(sams))
    y[match(rownames(Z1),sams)] <- 1
    
    obj <- SKAT_Null_Model(y ~ X, out_type="D")
    ## weighted
    #p0 <- SKAT(G, obj, kernel = "linear.weighted", weights=wts)
    #save(p0,file="SKATr0")
    #pV <- SKAT(G, obj)
    #save(pV,file="SKATr")
    
    ## weighted rare variants
#     pV <- rep(1,dim(G)[2])
#     for(i in 1:dim(G)[2]){
#         pV[i] <- SKAT(as.matrix(G[,i],ncol=1), obj, kernel = "linear.weighted", weights=wts[i])$p.value
#     }
#     
#     save(pV,file="singleSKATr")
    
    ## weighted combined and rare variants
    p0 <- SKAT_CommonRare(G, obj)$p.value
    save(p0,file="SKATr_RC")

    pV <- rep(1,dim(G)[2])
    for(i in 1:dim(G)[2]){
        pV[i] <- SKAT_CommonRare(as.matrix(G[,i],ncol=1), obj)$p.value
    }
    
    save(pV,file="singleSKATr_RC")
    
    
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
    load("caselist_7_28")
    casevar <- paste(caselist[,1],caselist[,2],sep="_")
    
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

    load("caselist_8_18")
    ### excluding unknown and none mutations
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss")                 
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    caselist <- caselist[caselist[,"VariantClass"] %in% c(lof,mis),]
    ### split damaging missense and misense: MetaSVM
    subs <- rep("MIS",dim(caselist)[1])
    subs[caselist[,"MetaSVM"]=="D"] <- "dMIS"
    subs[caselist[,"VariantClass"] %in% lof] <- "LOF"
    caselist <- cbind(caselist,subs)
    ncase <- length(unique(caselist[,"Subject_ID"]))
    
    load("contlist_8_6")
    ### filtered 
    contlist <- contlist[contlist[,"VariantClass"] %in% c(lof,mis),]
    ### split damaging missense and misense: MetaSVM
    subs <- rep("MIS",dim(contlist)[1])
    subs[contlist[,"MetaSVM"]=="D"] <- "dMIS"
    subs[contlist[,"VariantClass"] %in% lof] <- "LOF"
    contlist <- cbind(contlist,subs)
    
    ### ExAC 1% : 
    subs <- contlist[,"AlleleFrequency.ExAC"] < 0.01
    contlist <- contlist[subs,]
    
    ncont <- length(unique(contlist[,"SubID"]))
    
    a <- rep(0,5)
    a[1] <- dim(caselist)[1]
    a[2] <- sum(caselist[,"subs"]=="LOF")
    a[3] <- sum(caselist[,"subs"]=="dMIS")
    a[4] <- sum(caselist[,"subs"]=="MIS")
    a[5] <- a[2] + a[3]
    
    b <- rep(0,5)
    b[1] <- dim(contlist)[1]
    b[2] <- sum(contlist[,"subs"]=="LOF")
    b[3] <- sum(contlist[,"subs"]=="dMIS")
    b[4] <- sum(contlist[,"subs"]=="MIS")
    b[5] <- b[2] + b[3]
    
    NM <- matrix(0,5,5)
    for(i in 1:5){
        NM[i,1] <- a[i]
        NM[i,2] <- b[i]
        NM[i,3] <- (a[i]/ncase)/(b[i]/ncont)
        NM[i,4] <- binom.test(a[i],a[i]+b[i],ncase/(ncase+ncont))$p.value
        NM[i,5] <- fisher.test(matrix(c(a[i],b[i],ncase,ncont),2,2))$p.value
    }
    
    ##===========================================================================
    ## gene level test
    source("pre.R")
    geneburden(caselist,contlist)
    
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

## gene level burden test
geneburden <- function(caselist,contlist){
    ncase <- length(unique(caselist[,"Subject_ID"]))
    ncont <- length(unique(contlist[,"SubID"]))
    
    ngene <- length(unique(caselist[,"Gene"]))
    genes <- unique(caselist[,"Gene"])
    
    cols <- c("Gene","Variants","n.case","LOF_dMIS","n.case_1","case_LOF","case_dmis","case_mis","control_LOF","control_dmis","control_mis","p_LOF","p_dmis","p_mis","p_min")
    Sta <- matrix(0,ngene,15,dimnames=list(genes,cols))
    Sta[,1] <- genes
    Sta[,2] <- table(caselist[,"Gene"])[genes]    
    Sta[,3] <- sapply(1:length(genes),function(i) length(unique(caselist[caselist[,"Gene"] %in% genes[i],"Subject_ID"])))
    Sta[,4] <- table(caselist[caselist[,"subs"] %in% c("LOF","dMIS"),"Gene"])[genes]
    tmp <- caselist[caselist[,"subs"] %in% c("LOF","dMIS"),]
    Sta[,5] <- sapply(1:length(genes),function(i) length(unique(tmp[tmp[,"Gene"] %in% genes[i],"Subject_ID"])))
    
    vt <- c("LOF","dMIS","MIS")
    for(i in 1:3){
        tmp <- caselist[caselist[,"subs"]==vt[i],]
        Sta[,6+i-1] <- table(tmp[,"Gene"])[genes]
        tmp <- contlist[contlist[,"subs"]==vt[i],]
        Sta[,9+i-1] <- table(tmp[,"Gene"])[genes]
    }
    
    p <- ncase/(ncase+ncont)
    Sta[is.na(Sta)] <- 0
    tmp <- sapply(1:dim(Sta)[1],function(i){
        tmpn <- as.numeric(Sta[i,6:11])
#         a1 <- binom.test(tmpn[1],tmpn[1]+tmpn[4],p=p)$p.value
#         a2 <- binom.test(tmpn[2],tmpn[2]+tmpn[5],p=p)$p.value
#         a3 <- binom.test(tmpn[3],tmpn[3]+tmpn[6],p=p)$p.value
        a1 <- fisher.test(matrix(c(tmpn[1],tmpn[4],ncase,ncont),2,2))$p.value
        a2 <- fisher.test(matrix(c(tmpn[2],tmpn[5],ncase,ncont),2,2))$p.value
        a3 <- fisher.test(matrix(c(tmpn[3],tmpn[6],ncase,ncont),2,2))$p.value
        c(a1,a2,a3)
        })
    
    Sta[,12:14] <- t(tmp)
    Sta[,15] <- apply(t(tmp),1,min)
    Sta <- Sta[order(-as.numeric(Sta[,4])),]
    write.csv(Sta,file="GeneBurden_1.csv",row.names=FALSE)

}

## Panel genelists 
Panelg <- function(){
    burt <- read.csv("GeneBurden.csv")
    burg <- burt[1:100,1]
    
    filenames <- c("../genelist/Genelist1.txt","../genelist/Genelist2.txt","../genelist/Genelist3.txt","../genelist/CHEA_TFs.txt","../genelist/ENCODE_TFs.txt","../genelist/TRANSFAC_TFs.txt")
    for(i in 1:length(filenames)){
        g1 <- unlist(read.table(filenames[i]))
        print(length(intersect(g1,burg)))
    }
    
    g3 <- unlist(read.table(filenames[3]))
    qwt(intersect(g3,burg),file="Panel3_Bur.txt")
    
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

