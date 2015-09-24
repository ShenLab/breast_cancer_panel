### Step 1
case_cont <- function(){
    source("Faminfo.R")
    #load("caselist_8_26")
    load("caselist_9_15")
    caselist <- onlycase
    load("contlist_8_6")
    hotf <- "hotspots/hotf_cos_2"
    
    caselist <- remove_out(caselist)
    contlist <- remove_out(contlist)
    
    sig=FALSE
    caselist <- filter_variant(caselist,sig)
    write.csv(caselist[caselist[,"ExACfreq"],],file="IndexCasesvariants.csv",row.names=FALSE) # only write the rare variant based on ExAc frequency
    contlist <- filter_variant(contlist,sig)
    write.csv(contlist[contlist[,"ExACfreq"],],file="Controlvariants.csv",row.names=FALSE) # ontly write the rare variant based on ExAc Frequency
    singleana(caselist,contlist,sig,hotf)
    
    #load("caselist_8_26")
    load("caselist_9_15")
    caselist <- onlycase
    load("contlist_8_6")
    sig=TRUE
    caselist <- filter_variant(caselist,sig)
    write.csv(caselist[caselist[,"ExACfreq"],],file="IndexCasesvariants_single.csv",row.names=FALSE)
    contlist <- filter_variant(contlist,sig)
    write.csv(contlist[contlist[,"ExACfreq"],],file="Controlvariants_single.csv",row.names=FALSE)
    singleana(caselist,contlist,sig,hotf)
    
}

remove_out <- function(onelist){
    ## BRCA1/2 pathogenic mutations cases: br <- c(223109,223041,223275,260333,222968)
    br <- c(223109,223041,223275,260333,222968)
    if(grepl("SubID",colnames(onelist))){
        onelist <- onelist[!(onelist[,"SubID"] %in% br),]
    }else{
        onelist <- onelist[!(onelist[,"Subject_ID"] %in% br),]
    }
    onelist
}

singleana <- function(caselist,contlist,sig,hotf){
    
    ### filtered variants to analysis
    caselist <- caselist[caselist[,"Variantfiltering"],]
    contlist <- contlist[contlist[,"Variantfiltering"],]
        
    ### populations
    bc.pop <- read.delim("WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    ###table(bc.pop[,4])
    
    caseind <- nchar(caselist[,"REF"]) != nchar(caselist[,"ALT"])
    contind <- nchar(contlist[,"REF"]) != nchar(contlist[,"ALT"])
        
    ## LOF mutations analysis 
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss","none")                 
    case1 <- caselist[caselist[,"VariantClass"] %in% lof,]    
    cont1 <- contlist[contlist[,"VariantClass"] %in% lof,]
    strf <- paste("LOF",sig,sep="_")
    aa = geneburden_LOF(case1,cont1,bc.pop,strf)

    ## missense mutations in cosmic hotspots analysis
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    case2 <- caselist[caselist[,"VariantClass"] %in% mis & !caseind,]
    cont2 <- contlist[contlist[,"VariantClass"] %in% mis & !contind,]
    #bb = geneburden_MIS(case2,cont2,bc.pop,hotf)
    strf <- paste("MIS",sig,sep="_")
    bb = variantburden_MIS(case2,cont2,bc.pop,hotf,hotk=1,strf)

    ## indels analysis
    case4 <- caselist[caselist[,"VariantClass"] %in% mis & caseind,]
    cont4 <- contlist[contlist[,"VariantClass"] %in% mis & contind,]
    strf <- paste("indel",sig,sep="_")
    cc = variantburden_MIS(case4,cont4,bc.pop,hotf="",hotk=2,strf)

    ## any LOF and missense mutations in hotspots
    case3 <- rbind(case1,bb$caselist,case4)
    cont3 <- rbind(cont1,bb$contlist,cont4)
    
    save(case3,file=paste("case3_single",sig,sep=""))
    save(cont3,file=paste("cont3_single",sig,sep=""))
    
    ## missense mutations analysis
    #hotf <- ""
    #cc = geneburden_MIS(case2,cont2,bc.pop,hotf,hotk=2)

} 

filter_variant <- function(onelist,sig=FALSE){
    
    filters <- c("Variantfiltering","ExACfreq","Popfreq","VCFPASS","noneSegmentalDup","meta-SVM_PP2","GTEXexp","singleton")
    filS <- matrix(FALSE,dim(onelist)[1],length(filters))
    colnames(filS) <- filters
    onelist <- cbind(onelist,filS)
    
    # EXAC cut <- 0.001 ## rare variants
    #Ecut=0.001;topf <- 1;
    Ecut=0.01;topf <- 0.5;
    
    onelist[onelist[,"AlleleFrequency.ExAC"]< Ecut,"ExACfreq"] <- TRUE
    
    # population frequency < 5%
    onelist[,"Popfreq"] <- TRUE
    Pcut <- 0.05
    Pn <- 1160 * 0.05 ## 1160 population
    load("varT_9_10")
    vs <- paste(onelist[,1],onelist[,2],onelist[,4],onelist[,5],sep="_")
    subs1 <- rep(FALSE,length(vs));
    iv <- intersect(vs,names(varT))
    subs1[match(iv,vs)] <- varT[match(iv,names(varT))] >= Pn
    onelist[subs1,"Popfreq"] <- FALSE
    
    ### filtered more details
    subs1 <- onelist[,"FILTER"] == "PASS"
    #subs2 <- onelist[,"SegmentalDuplication"] == "none"
    subs2 <- sapply(1:dim(onelist)[1], function(i) {
        if(onelist[i,"SegmentalDuplication"] == "none"){ TRUE;
        }else{ tmp <- unlist(strsplit(onelist[i,"SegmentalDuplication"],","))[1]
               as.numeric(unlist(strsplit(tmp,":"))[2]) < 0.95
        }
        })
    onelist[subs1,"VCFPASS"] <- TRUE
    onelist[subs2,"noneSegmentalDup"] <- TRUE
    
    ### missense predicted by meta-SVM and PP2
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    onelist[,"meta-SVM_PP2"] <- TRUE
    subs2 <- onelist[,"VariantClass"] %in% mis
    subs1 <- onelist[,"PP2prediction"]=="D" | onelist[,"MetaSVM"]=="D"
    ###subs1 <- onelist[,"MetaSVM"]=="D"
    subs3 <- nchar(onelist[,"REF"]) != nchar(onelist[,"ALT"]) ### indels
    subs2 <- subs2 & !subs3
    onelist[subs2 & !subs1,"meta-SVM_PP2"] <- FALSE
    
    ### expressed genes in GTEX or ranked genes by TCGA (or 29 inherited genes)
    ##expg <- read.table("../geneexpression/GTEx/expg_ranked.txt")[,1]
    expg <- read.table("data/TCGA57g.txt")[,1]
    expg <- expg[1:(length(expg)*topf)]
    onelist[onelist[,"Gene"] %in% expg,"GTEXexp"] <- TRUE
    #a <- read.csv("../genelist/Tumor_supressor/TS.csv")[,1]
    #length(intersect(a,expg)) ## 18114
    
    onelist[,"singleton"] <- TRUE
    if(sig){
        ### singleton variants only !!!!  
        vs <- paste(onelist[,1],onelist[,2],onelist[,4],onelist[,5],sep="_") ####!!!!
        vsdup <- vs[duplicated(vs)]
        subs1 <- vs %in% vsdup
        #subs2 <- onelist[,"VariantClass"] %in% lof
        onelist[subs1, "singleton"] <- FALSE
    }
    
    ### no dbsnp
    ### onelist <- onelist[onelist[,"ID"]==".",]
    
    onelist[,"Variantfiltering"] <- onelist[,"ExACfreq"] & onelist[,"Popfreq"] & onelist[,"VCFPASS"] & onelist[,"noneSegmentalDup"] & onelist[,"meta-SVM_PP2"] & onelist[,"GTEXexp"] & onelist[,"singleton"]
    
    #### write log file for each version
    con <- file(paste("variant_filtering_single_",sig,".log",sep=""),"w");
    writeLines(paste("Filtering log: ",date(),sep=""),con)
    writeLines(paste("ExAC frequency: ",Ecut," (ExAc",Ecut,")",sep=""),con)
    writeLines(paste("Population frequency: ",Pcut," (Pop",Pcut,")",sep=""),con)
    #writeLines(paste("LOF and missense mutations (LOF: frameshiftdeletion,frameshiftinsertion,stopgain,stoploss,splicing; Mis: nonframeshiftdeletion,nonframeshiftinsertion,nonsynonymousSNV) "," (LOF_MIS)",sep=""),con)
    writeLines(paste("VCF FILTER: PASS"," (VCFPASS)",sep=""),con)
    writeLines(paste("VCF SegmentalDuplication: < 0.95"," (noneSegmentalDup)",sep=""),con)
    writeLines(paste("Meta-SVM or PP2 prediction: D for missense mutations"," (meta-SVM_PP2)",sep=""),con)
    writeLines(paste("GTEX expressed genes "," (GTExexp)",sep=""),con)
    if(sig){
        writeLines(paste("Singleton variants "," (singleton)",sep=""),con)
    }
    close(con)
    
    #### write updated log file for each version
#     con <- file("variant_filtering_updated.log","w");
#     writeLines(paste("Updated filtering log: ",date(),sep=""),con)
#     writeLines(paste("add LOF with splicing; ","(LOF_MIS)",sep=""),con)
#     writeLines(paste("Stricted Meta-SVM and PP2 predictions as: (D .),(.,D),(D,D),(.,.) for missense mutations excluding predicted as (.,T),(T,.)"," (meta-SVM_PP2)",sep=""),con)
#     writeLines(paste("add multiallele LOF"," (singleton)",sep=""),con)
#     writeLines(paste("burden test for each missense mutation"," (missense mutation)",sep=""),con)
#     close(con)
    
    onelist
}

after_check <- function(){
    
    sig <- FALSE
    if(sig){
        load("case3_singleTRUE")
        load("cont3_singleTRUE")
    }else{
        load("case3_singleFALSE")
        load("cont3_singleFALSE")
    }
    #fils <- read.csv("filtered.csv")
    fils <- read.csv("Famcheck.csv")
    
    vfil <- paste(fils[,1],fils[,2],fils[,3],fils[,6],fils[,7],sep="_")
    vcase <- paste(case3[,"Gene"],case3[,1],case3[,2],case3[,"REF"],case3[,"ALT"],sep="_")
    vcont <- paste(cont3[,"Gene"],cont3[,1],cont3[,2],cont3[,"REF"],cont3[,"ALT"],sep="_")
    case3 <- case3[vcase %in% vfil,]
    #cont3 <- cont3[vcont %in% vfil,]
    ##==============================================================================
    caselist <- case3
    contlist <- cont3
    ### populations
    bc.pop <- read.delim("WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    ###table(bc.pop[,4])
    
    caseind <- nchar(caselist[,"REF"]) != nchar(caselist[,"ALT"])
    contind <- nchar(contlist[,"REF"]) != nchar(contlist[,"ALT"])
    
    ## LOF mutations analysis 
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss","none")                 
    case1 <- caselist[caselist[,"VariantClass"] %in% lof,]    
    cont1 <- contlist[contlist[,"VariantClass"] %in% lof,]
    strf <- paste("LOFcheck",sig,sep="_")
    aa = geneburden_LOF(case1,cont1,bc.pop,strf)
    
    ## missense mutations in cosmic hotspots analysis
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    case2 <- caselist[caselist[,"VariantClass"] %in% mis & !caseind,]
    cont2 <- contlist[contlist[,"VariantClass"] %in% mis & !contind,]
    hotf <- "hotspots/hotf_cos"
    strf <- paste("MIScheck",sig,sep="_")
    bb = variantburden_MIS(case2,cont2,bc.pop,hotf,hotk=1,strf)
    
    ## indels analysis
    case4 <- caselist[caselist[,"VariantClass"] %in% mis & caseind,]
    cont4 <- contlist[contlist[,"VariantClass"] %in% mis & contind,]
    strf <- paste("indelcheck",sig,sep="_")
    cc = variantburden_MIS(case4,cont4,bc.pop,hotf="",hotk=2,strf)
    
}

## lof and missense mutation burden test
geneburden_LOF <- function(caselist,contlist,bc.pop,strf){
    
    cols <- c("Gene","UniLOF","n.case","case_LOF","control_LOF","odds_ratio","p_value","case_J_LOF","control_J_LOF","J_odds_ratio","J_p_value","case_H_LOF","control_H_LOF","H_odds_ratio","H_p_value")
    result <- geneburden(caselist,contlist,bc.pop,cols)
    Sta <- result$Sta
    n.case <- result$n.case
    n.cont <- result$n.cont
    write.csv(Sta,file=paste("GeneBurden_",strf,".csv",sep=""),row.names=FALSE)
    
    list(n.case=n.case,n.cont=n.cont)
}

geneburden_MIS <- function(caselist,contlist,bc.pop,hotf,hotk=1,strf){
    
    if(hotk==1){
        Fil <- mis_hotspot(hotf,caselist,contlist)
        caselist <- Fil$caselist
        contlist <- Fil$contlist
    }
    
    cols <- c("Gene","UniMIS","n.case","case_MIS","control_MIS","odds_ratio","p_value","case_J_MIS","control_J_MIS","J_odds_ratio","J_p_value","case_H_MIS","control_H_MIS","H_odds_ratio","H_p_value")
    result <- geneburden(caselist,contlist,bc.pop,cols)
    Sta <- result$Sta
    n.case <- result$n.case
    n.cont <- result$n.cont
    
    if(hotk==1){
        write.csv(Sta,file=paste("GeneBurden_hotspot_",strf,".csv",sep=""),row.names=FALSE)
    }else{
        write.csv(Sta,file=paste("GeneBurden_",strf,".csv",sep=""),row.names=FALSE)
    }
    list(n.case=n.case,n.cont=n.cont,caselist=caselist,contlist=contlist)
    
}

variantburden_MIS <- function(caselist,contlist,bc.pop,hotf,hotk=1,strf){
    
    if(hotk==1){
        Fil <- mis_hotspot(hotf,caselist,contlist)
        caselist <- Fil$caselist
        contlist <- Fil$contlist
    }
    
    cols <- c("Variant","Gene","n.case","case_MIS","control_MIS","odds_ratio","p_value","case_J_MIS","control_J_MIS","J_odds_ratio","J_p_value","case_H_MIS","control_H_MIS","H_odds_ratio","H_p_value")
    result <- variantburden(caselist,contlist,bc.pop,cols)
    Sta <- result$Sta
    n.case <- result$n.case
    n.cont <- result$n.cont
    
    if(hotk==1){
        write.csv(Sta,file=paste("VariantBurden_hotspot_",strf,".csv",sep=""),row.names=FALSE)
    }else{
        write.csv(Sta,file=paste("VariantBurden_",strf,".csv",sep=""),row.names=FALSE)
    }
    list(n.case=n.case,n.cont=n.cont,caselist=caselist,contlist=contlist)
    
}

mis_hotspot <- function(hotf,caselist,contlist){

    load(hotf)
    hotg <- names(hots)
    #caselist <- caselist[caselist[,"Gene"] %in% hotg,]
    #contlist <- contlist[contlist[,"Gene"] %in% hotg,]
    
    ### only missense mutations in hotspots and indels
    p.case <- sapply(1:dim(caselist)[1], function(i) {
        tmp <- unlist(strsplit(caselist[i,"AAchange"],":"))[5];
        if(grepl("_",tmp) | grepl("del",tmp) | grepl("ins",tmp)){ TRUE; }else{ gsub("\\D","",tmp);}
    })
    s.case <- sapply(1:length(p.case),function(i){
        if(p.case[i]!=TRUE){
            k <- which(hotg==caselist[i,"Gene"])
            if(length(k)>0){p.case[i] %in% hots[[k]];
            }else{FALSE;}
        }else{TRUE;}
    })
    caselist <- caselist[s.case,]
    
    p.cont <- sapply(1:dim(contlist)[1], function(i) {
        tmp <- unlist(strsplit(contlist[i,"AAchange"],":"))[5];
        if(grepl("_",tmp) | grepl("del",tmp) | grepl("ins",tmp)){ TRUE; }else{ gsub("\\D","",tmp);}
    })
    s.cont <- sapply(1:length(p.cont),function(i){
        if(p.cont[i]!=TRUE){
            k <- which(hotg==contlist[i,"Gene"])
            if(length(k)>0){p.cont[i] %in% hots[[k]];
            }else{FALSE;}
        }else{TRUE;}
    })
    contlist <- contlist[s.cont,]
    print(dim(caselist))
    print(dim(contlist))  
    
    list(caselist=caselist,contlist=contlist)

}

geneburden <- function(caselist,contlist,bc.pop,cols){
    
    ##genes <- union(caselist[,"Gene"],contlist[,"Gene"])
    genes <- unique(caselist[,"Gene"])
    ngene <- length(genes)
    Sta <- matrix(0,ngene,length(cols),dimnames=list(genes,cols))
    
    caseV <- paste(caselist[,1],caselist[,2],sep="_")
    contV <- paste(contlist[,1],contlist[,2],sep="_") 
    
    Sta[,2] <- sapply(1:length(genes),function(i) length(unique(caseV[caselist[,"Gene"] %in% genes[i]])))
    ##Sta[,3] <- sapply(1:length(genes),function(i) length(unique(contV[contlist[,"Gene"] %in% genes[i]])))
    Sta[,3] <- sapply(1:length(genes),function(i) length(unique(caselist[caselist[,"Gene"] %in% genes[i],"Subject_ID"])))
    
    Sta[,4] <- table(caselist[,"Gene"])[genes]
    Sta[,5] <- table(contlist[,"Gene"])[genes]
    
    Jp <-  bc.pop[bc.pop[,4] %in% "J",3]
    Sta[,8] <- table(caselist[caselist[,"Subject_ID"] %in% Jp,"Gene"])[genes]
    Sta[,9] <- table(contlist[contlist[,"SubID"] %in% Jp,"Gene"])[genes]
    
    Hp <-  bc.pop[bc.pop[,4] %in% "H",3]
    Sta[,12] <- table(caselist[caselist[,"Subject_ID"] %in% Hp,"Gene"])[genes]
    Sta[,13] <- table(contlist[contlist[,"SubID"] %in% Hp,"Gene"])[genes]
    
    caseid <- unique(caselist[,"Subject_ID"])
    contid <- unique(contlist[,"SubID"])
    #n.case <- c(length(caseid),length(intersect(caseid,Jp)),length(intersect(caseid,Hp)))
    #n.cont <- c(length(contid),length(intersect(contid,Jp)),length(intersect(contid,Hp)))
    n.case <- c(356,223,99)
    n.cont <- c(114,59,55)
    
    Sta[is.na(Sta)] <- 0
    for(i in 1:3){
        ntmp <- sapply(1:length(genes), function(j){
            tmp <- fisher.test(matrix(c(Sta[j,4*i],Sta[j,4*i+1],n.case[i],n.cont[i]),2,2));
            c(tmp$estimate,tmp$p.value)
        })
        Sta[,c(4*i+2,4*i+3)] <- t(ntmp)
    }
    
    Sta[,1] <- genes
    Sta <- Sta[order(-as.numeric(Sta[,4])),]
    
    list(Sta=Sta,n.case=n.case,n.cont=n.cont)
}

variantburden <- function(caselist,contlist,bc.pop,cols){
    
    allv <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
    contv <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
    vars <- unique(allv)
    
    nvar <- length(vars)
    Sta <- matrix(0,nvar,length(cols),dimnames=list(vars,cols))
    
    Sta[,3] <- sapply(1:length(vars),function(i) length(unique(caselist[allv %in% vars[i],"Subject_ID"])))
    
    Sta[,4] <- table(allv)[vars]
    Sta[,5] <- table(contv)[vars]
    
    Jp <-  bc.pop[bc.pop[,4] %in% "J",3]
    Sta[,8] <- table(allv[caselist[,"Subject_ID"] %in% Jp])[vars]
    Sta[,9] <- table(contv[contlist[,"SubID"] %in% Jp])[vars]
    
    Hp <-  bc.pop[bc.pop[,4] %in% "H",3]
    Sta[,12] <- table(allv[caselist[,"Subject_ID"] %in% Hp])[vars]
    Sta[,13] <- table(contv[contlist[,"SubID"] %in% Hp])[vars]
    
    caseid <- unique(caselist[,"Subject_ID"])
    contid <- unique(contlist[,"SubID"])
    #n.case <- c(length(caseid),length(intersect(caseid,Jp)),length(intersect(caseid,Hp)))
    #n.cont <- c(length(contid),length(intersect(contid,Jp)),length(intersect(contid,Hp)))
    n.case <- c(356,223,99)
    n.cont <- c(114,59,55)
    
    Sta[is.na(Sta)] <- 0
    for(i in 1:3){
        ntmp <- sapply(1:length(vars), function(j){
            tmp <- fisher.test(matrix(c(Sta[j,4*i],Sta[j,4*i+1],n.case[i],n.cont[i]),2,2));
            c(tmp$estimate,tmp$p.value)
        })
        Sta[,c(4*i+2,4*i+3)] <- t(ntmp)
    }
    
    Sta[,1] <- vars
    Sta[,2] <- caselist[match(vars,allv),"Gene"]
    Sta <- Sta[order(-as.numeric(Sta[,4])),]
    
    list(Sta=Sta,n.case=n.case,n.cont=n.cont)

}

### Step 2
### related gene set burden analysis
geneset_burden <- function(){
    ## gene set Panel burden 
    ts <- unlist(read.table("hotspots/TScell_filtered.txt")) ## tumor suppresspor
    drs <- unlist(read.table("hotspots/Driver_filtered.txt")) ## cancer drivers
    dna <- unlist(read.table("../genelist/allrepair_gene.txt")) ## DNA repair genes from KEGG
    Cicc <- read.csv("BC_Candidategene_Summary_09262014 for Ciccia.csv")
    Cis <- Cicc[Cicc[,"Priority"]=="Likely",1]  ## Cicca likely genes
    
    ts <- gsub("Sep-05","SEPT5",ts);ts <- gsub("Sep-06","SEPT6",ts);ts <- gsub("Sep-09","SEPT9",ts);
    drs <- gsub("Sep-05","SEPT5",drs);drs <- gsub("Sep-06","SEPT6",drs);drs <- gsub("Sep-09","SEPT9",drs);
    dna <- gsub("Sep-05","SEPT5",dna);dna <- gsub("Sep-06","SEPT6",dna);dna <- gsub("Sep-09","SEPT9",dna);
    Cis <- gsub("Sep-05","SEPT5",Cis);Cis <- gsub("Sep-06","SEPT6",Cis);Cis <- gsub("Sep-09","SEPT9",Cis);
    
#     allg <- union(ts,union(drs,union(dna,Cis)))
#     allg <- gsub("Sep-05","SEPT5",allg);allg <- gsub("Sep-06","SEPT6",allg);allg <- gsub("Sep-09","SEPT9",allg);  
#     expg <- unlist(read.table("../geneexpression/GTEx/expg.txt"))
#     print(sum(allg %in% expg))
#     print(length(unique(allg)))
#     a <- setdiff(allg,expg)
#     qwt(a,file="un_express.txt")
    
    strnames <- c("LGD","MIS","indel")
    genestr <- c("TS","driver","DNA_repair","Ciccia")
    sig=c(FALSE,TRUE)
    
    genel <- list()
    genel[[1]] <- ts; genel[[2]] <- drs; genel[[3]] <- dna; genel[[4]] <- Cis
    
    for(i in 1:2){
        filenames <- c(paste("GeneBurden_LOF_",sig[i],".csv",sep=""),paste("VariantBurden_hotspot_MIS_",sig[i],".csv",sep=""),paste("VariantBurden_indel_",sig[i],".csv",sep=""))
        lgdl <- c()
        varl <- c()
        for(j in 1:3){
            for(k in 1:4){
                cc <- write_panel(genel[[k]],filenames[j], paste("details/",strnames[j],"_",genestr[k],"_",sig[i],"_burden.csv",sep="")) 
                if(j==1){
                    lgdl <- rbind(lgdl,cc)
                }else{
                    varl <- rbind(varl,cc)
                }
            }
        }
        onef <- paste(lgdl[,1],lgdl[,2],lgdl[,3],lgdl[,4],lgdl[,5],sep="_")
        tmp <- lgdl[match(unique(onef),onef),]
        write.csv(tmp,file=paste("ALLcheckLOF_uni_",sig[i],".csv",sep=""),row.names=FALSE)
        
        onef <- paste(varl[,1],varl[,2],varl[,3],varl[,4],varl[,5],sep="_")
        tmp <- varl[match(unique(onef),onef),]
        write.csv(tmp,file=paste("ALLcheckMIS_indel_uni_",sig[i],".csv",sep=""),row.names=FALSE)
    }

}

filter_gene <- function(){
    ## tumor suppressors 
    ts1 <- unlist(read.table("../genelist/Genelist2.txt"))
    ts2  <- unlist(read.table("../genelist/Tumor_supressor/TS_Cosmic.txt")) ## tumor supressors
    ts3 <- unlist(read.table("hotspots/TScell_filtered.txt"))
    ts <- union(ts1,union(ts2,ts3))
    
    write.table(ts,file="hotspots/TS_collected.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
    
    ### COSMIC LOF percentage >= 15%
    cos <- read.delim("hotspots/CosmicMutantExport.tsv",sep="\t")
    cos1 <- cos[cos[,"Gene.name"] %in% ts,]
    #save(cos1,file="hotspots/TScheck")
    #load("hotspots/TScheck")
    genes <- unique(cos1[,1])
    subnon <- cos1[,"Mutation.Description"] == "Substitution - Nonsense"
    nonfrac <- sapply(1:length(genes),function(i){
        subs1 <- cos1[,"Gene.name"] == genes[i]
        subs2 <- subs1 & subnon
        length(unique(cos1[subs2,"Sample.name"]))/ length(unique(cos1[subs1,"Sample.name"]))
    })
    genes <- genes[nonfrac>=quantile(nonfrac,0.5)]
    
    ### at least one LOF in TCGA somatic mutation
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
    muta <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Nonsense_Mutation","Nonstop_Mutation","Splice_Site")
    tcgaV <- tcgaV[tcgaV[,"Variant_Classification"] %in% muta,]
    
    genes <- genes[genes %in% tcgaV[,1]]
    write.table(genes,file="hotspots/TS_filtered.txt")
    
    tsl <- read.csv("../genelist/Tumor_supressor/TS.csv")
    ts <- tsl[tsl[,"TUSON_q_value_TSG"] <= 0.1,1]
    write.table(ts,file="hotspots/TScell_filtered.txt")
    
    ### cancer driver genes
    ## set 1: cancer census genes
    filename <- "../somatic_mutation/Somatic_Mutations_COSMIC/cancer_gene_census.csv"
    cosV <- read.csv(filename,check.names=FALSE)
    set1 <- unique(cosV[,1])
    ## set 2: mouse candidate genes: CCGD
    filename <- "../genelist/cancer_driver/CCGD_export.csv"
    ccgd <- read.csv(filename)
    set2 <- ccgd[ccgd[,"Relative.Rank"] %in% "A","Human.Symbol"]
    set2 <- unique(set2)
    ## set 3: cbioportal
    set3 <- unlist(read.table("../genelist/cancer_driver/all_MutSig_Q-lt-0.1.genes.txt"))
    
    drs <- union(set1,set3)
    write.table(drs,file="hotspots/Driver_filtered.txt")
    
}

write_panel <- function(genes,burdenf, wf){
    n.case <- c(356,223,99)
    n.cont <- c(114,59,55)
    
    aa <- read.csv(burdenf)
    cc <- aa[aa[,"Gene"] %in% genes,]
    cols <- names(cc)
    #cc <- cc[order(cc[,"p_value"]),]
    tmp <- rep(0,dim(cc)[2])
    if(names(cc)[1]=="Gene"){
        tmp[c(2:5,8:9,12:13)] <- colSums(cc[,c(2:5,8:9,12:13)]);   
    }else{
        tmp[c(3:5,8:9,12:13)] <- colSums(cc[,c(3:5,8:9,12:13)]);
    }
    
    tmp[6] <- (tmp[4]/n.case[1]) / (tmp[5]/n.cont[1])
    tmp[10] <- (tmp[8]/n.case[2]) / (tmp[9]/n.cont[2])
    tmp[14] <- (tmp[12]/n.case[3]) / (tmp[13]/n.cont[3])
    tmp[7] <-  ifelse(tmp[4]+tmp[5] > 0, binom.test(tmp[4],tmp[4]+tmp[5],n.case[1]/(n.case[1]+n.cont[1]))$p.value,1)
    tmp[11] <-  ifelse(tmp[8]+tmp[9] > 0, binom.test(tmp[8],tmp[8]+tmp[9],n.case[2]/(n.case[2]+n.cont[2]))$p.value,1)
    tmp[15] <-  ifelse(tmp[12]+tmp[13] >0, binom.test(tmp[12],tmp[12]+tmp[13],n.case[3]/(n.case[3]+n.cont[3]))$p.value,1)
    tmp[1] <- "Total"      
    names(tmp) <- names(cc)
    cc <- rbind(cc,tmp)
    names(cc) <- cols
    write.csv(cc,file=wf,row.names=FALSE)
    
    cc
}

### Step 3  variant type burden analysis
variant_type_burden0 <- function(){
    
    ts <- unlist(read.table("hotspots/TS_collected.txt")) ## tumor suppresspor
    drs <- unlist(read.table("hotspots/Driver_filtered.txt")) ## cancer drivers
    
    #agenes=ts;
    
    source("Faminfo.R")
    for(sig in c(FALSE,TRUE)){
    load(paste("case3_single",sig,sep=""))
    load(paste("cont3_single",sig,sep=""))
    
    a <- variant_type_burden_sig(case3,cont3)
    a[,1] <- paste(a[,1],"(",a[,5],")",sep="")
    a[,2] <- paste(a[,2],"(",a[,6],")",sep="")
    print(a)
    print(length(union(case3[,"Gene"],cont3[,"Gene"])))
    
    agenes=ts;
    case1 <- case3[case3[,"Gene"] %in% agenes,]
    cont1 <- cont3[cont3[,"Gene"] %in% agenes,]
    a <- variant_type_burden_sig(case1,cont1)
    a[,1] <- paste(a[,1],"(",a[,5],")",sep="")
    a[,2] <- paste(a[,2],"(",a[,6],")",sep="")
    print(a)
    
    agenes=drs;
    case2 <- case3[case3[,"Gene"] %in% agenes,]
    cont2 <- cont3[cont3[,"Gene"] %in% agenes,]
    a <- variant_type_burden_sig(case2,cont2)
    a[,1] <- paste(a[,1],"(",a[,5],")",sep="")
    a[,2] <- paste(a[,2],"(",a[,6],")",sep="")
    print(a)    
    }
    
    #===================synonmous===============================================
    source("Faminfo.R")
    load("casesy_9_15")
    load("contsy_9_15")

    for(sig in c(FALSE,TRUE)){
    casesy <- filter_variant(casesy,sig)
    contsy <- filter_variant(contsy,sig)
    casesy <- casesy[casesy[,"Variantfiltering"],]
    contsy <- contsy[contsy[,"Variantfiltering"],]
    cc <- synonymous_burden(casesy,contsy)
    cc[,1] <- paste(cc[,1],"(",cc[,5],")",sep="")
    cc[,2] <- paste(cc[,2],"(",cc[,6],")",sep="")
    print(cc)
        
    casesy1 <- casesy[casesy[,"Gene"] %in% ts,]
    contsy1 <- contsy[contsy[,"Gene"] %in% ts,]
    cc <- synonymous_burden(casesy1,contsy1)
    cc[,1] <- paste(cc[,1],"(",cc[,5],")",sep="")
    cc[,2] <- paste(cc[,2],"(",cc[,6],")",sep="")
    print(cc)
    
    casesy2 <- casesy[casesy[,"Gene"] %in% drs,]
    contsy2 <- contsy[contsy[,"Gene"] %in% drs,]
    cc <- synonymous_burden(casesy2,contsy2)
    cc[,1] <- paste(cc[,1],"(",cc[,5],")",sep="")
    cc[,2] <- paste(cc[,2],"(",cc[,6],")",sep="")
    print(cc)    
    }
}

synonymous_burden <- function(casesy,contsy){
    n.case <- c(356,223,99)
    n.cont <- c(114,59,55)
    
    ### population
    bc.pop <- read.delim("WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    Jp <-  bc.pop[bc.pop[,4] %in% "J",3]
    Hp <-  bc.pop[bc.pop[,4] %in% "H",3]
    
    cc <- matrix(0,3,6)
    i=1
    list1 <- casesy
    list2 <- contsy
    cc[i,1] <- dim(list1)[1]
    cc[i,2] <- dim(list2)[1]
    cc[i,3] <- (cc[i,1]/n.case[i]) / (cc[i,2]/n.cont[i])
    cc[i,4] <- binom.test(cc[i,1],cc[i,1]+cc[i,2],n.case[i]/(n.case[i]+n.cont[i]))$p.value
    cc[i,5] <- length(unique(list1[,"Gene"]))
    cc[i,6] <- length(unique(list2[,"Gene"]))
    
    i=2
    list1 <- casesy[casesy[,"Subject_ID"] %in% Jp,]
    list2 <- contsy[contsy[,"Subject_ID"] %in% Jp,]
    cc[i,1] <- dim(list1)[1]
    cc[i,2] <- dim(list2)[1]
    cc[i,3] <- (cc[i,1]/n.case[i]) / (cc[i,2]/n.cont[i])
    cc[i,4] <-  binom.test(cc[i,1],cc[i,1]+cc[i,2],n.case[i]/(n.case[i]+n.cont[i]))$p.value
    cc[i,5] <- length(unique(list1[,"Gene"]))
    cc[i,6] <- length(unique(list2[,"Gene"]))
    
    i=3
    list1 <- casesy[casesy[,"Subject_ID"] %in% Hp,]
    list2 <- contsy[contsy[,"Subject_ID"] %in% Hp,]    
    cc[i,1] <- dim(list1)[1]
    cc[i,2] <- dim(list2)[1]
    cc[i,3] <- (cc[i,1]/n.case[i]) / (cc[i,2]/n.cont[i])
    cc[i,4] <-  binom.test(cc[i,1],cc[i,1]+cc[i,2],n.case[i]/(n.case[i]+n.cont[i]))$p.value
    cc[i,5] <- length(unique(list1[,"Gene"]))
    cc[i,6] <- length(unique(list2[,"Gene"]))

    cc
}

variant_type_burden_sig <- function(caselist,contlist){
    
    n.case <- c(356,223,99)
    #n.case <- c(349,219,95)
    n.cont <- c(114,59,55)
    # LOF, MIS, indels, slient
    subcase <- variant_types(caselist)
    subcont <- variant_types(contlist)
    
    cc <- matrix(0,9,6)
    for(i in 1:3){
        cc[i,1] <- sum(subcase[,i])
        cc[i,2] <- sum(subcont[,i])
        cc[i,3] <- (cc[i,1]/n.case[1]) / (cc[i,2]/n.cont[1])
        cc[i,4] <-  ifelse(cc[i,1]+cc[i,2] >0 ,binom.test(cc[i,1],cc[i,1]+cc[i,2],n.case[1]/(n.case[1]+n.cont[1]))$p.value, 1)
        cc[i,5] <- length(unique(caselist[subcase[,i],"Gene"]))
        cc[i,6] <- length(unique(caselist[subcont[,i],"Gene"]))
    }
    
    ### population
    bc.pop <- read.delim("WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    Jp <-  bc.pop[bc.pop[,4] %in% "J",3]
    Hp <-  bc.pop[bc.pop[,4] %in% "H",3]
    
    case1 <- caselist[caselist[,"Subject_ID"] %in% Jp,]
    cont1 <- contlist[contlist[,"SubID"] %in% Jp,]
    case2 <- caselist[caselist[,"Subject_ID"] %in% Hp,]
    cont2 <- contlist[contlist[,"SubID"] %in% Hp,]    
    
    subcase <- variant_types(case1)
    subcont <- variant_types(cont1)
    for(i in 4:6){
        cc[i,1] <- sum(subcase[,i-3])
        cc[i,2] <- sum(subcont[,i-3])
        cc[i,3] <- (cc[i,1]/n.case[2]) / (cc[i,2]/n.cont[2])
        cc[i,4] <-  ifelse(cc[i,1]+cc[i,2] > 0, binom.test(cc[i,1],cc[i,1]+cc[i,2],n.case[2]/(n.case[2]+n.cont[2]))$p.value, 1)
        cc[i,5] <- length(unique(case1[subcase[,i-3],"Gene"]))
        cc[i,6] <- length(unique(cont1[subcont[,i-3],"Gene"]))
    }
    
    subcase <- variant_types(case2)
    subcont <- variant_types(cont2)
    for(i in 7:9){
        cc[i,1] <- sum(subcase[,i-6])
        cc[i,2] <- sum(subcont[,i-6])
        cc[i,3] <- (cc[i,1]/n.case[3]) / (cc[i,2]/n.cont[3])
        cc[i,4] <-  ifelse(cc[i,1]+cc[i,2]>0, binom.test(cc[i,1],cc[i,1]+cc[i,2],n.case[3]/(n.case[3]+n.cont[3]))$p.value, 1)
        cc[i,5] <- length(unique(case2[subcase[,i-6],"Gene"]))
        cc[i,6] <- length(unique(cont2[subcont[,i-6],"Gene"]))
    }
    
    cc
}

variant_types <- function(onelist){
    
    # LOF, MIS, indels,  slient(others unknown)
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss","none")   
    tmpsub <- nchar(onelist[,"REF"]) != nchar(onelist[,"ALT"])
    subM <- matrix(FALSE,dim(onelist)[1],4)
    subM[,1] <- onelist[,"VariantClass"] %in% lof
    subM[,2] <- (onelist[,"VariantClass"] %in% mis) & !tmpsub
    subM[,3] <- (onelist[,"VariantClass"] %in% mis) & tmpsub
    subM[,4] <- !subM[,1] & !subM[,2] & !subM[,3]
        
    subM
}

variant_type_burden <- function(){
    
    n.case <- c(345,216,94)
    n.cont <- c(114,59,55)
    
    ## all class: LOF and missense in hotspots 
    loff <- "GeneBurden_LOF.csv"
    misf <- "VariantBurden_MIS_hotspot.csv"
    cc <- write_variant(loff,misf,n.case,n.cont)
    print(cc)

    #     ## LOF in tumor suppressors
    #     loff <- "LGD_TS_burden.csv"
    #     ## missense in tumor suppressors
    #     misf <- "MIShotspot_TS_burden.csv"
    #     cc <- write_variant(loff,misf,n.case,n.cont)
    #     print(cc)
    #     
    #     ## LOF in cancer dirver
    #     loff <- "LGD_driver_burden.csv"
    #     ## missense in cancer dirver
    #     misf <- "MIShotspot_driver_burden.csv"
    #     cc <- write_variant(loff,misf,n.case,n.cont)
    #     print(cc)
}

write_variant <- function(loff,misf,n.case,n.cont){
    
    aa <- read.csv(loff)
    bb <- read.csv(misf)
    cc <- matrix(0,9,4)
    cc[,1] <- c(sum(aa[,4]),sum(bb[,4]),sum(c(aa[,4],bb[,4])), sum(aa[,8]),sum(bb[,8]),sum(c(aa[,8],bb[,8])), sum(aa[,12]),sum(bb[,12]),sum(c(aa[,12],bb[,12])))
    cc[,2] <- c(sum(aa[,5]),sum(bb[,5]),sum(c(aa[,5],bb[,5])), sum(aa[,9]),sum(bb[,9]),sum(c(aa[,9],bb[,9])), sum(aa[,13]),sum(bb[,13]),sum(c(aa[,13],bb[,13])))
    for(i in 1:9){
        j <- floor((i-1)/3) + 1
        if(cc[i,1]+cc[i,2] > 0){
            cc[i,3] <- (cc[i,1]/n.case[j]) / (cc[i,2]/n.cont[j])
            cc[i,4] <- binom.test(cc[i,1],cc[i,1]+cc[i,2],n.case[j]/(n.case[j]+n.cont[j]))$p.value
        }
    }
    
    cc
}

### Step 4
check_genes <- function(){
    
    #filenames <- c("ALLcheckLOF_uni_TRUE.csv","ALLcheckMIS_indel_uni_TRUE.csv")
    #load("case3_singleTRUE")
    #load("cont3_singleTRUE")
    
    filenames <- c("ALLcheckLOF_uni_FALSE.csv","ALLcheckMIS_indel_uni_FALSE.csv")
    load("case3_singleFALSE")
    load("cont3_singleFALSE")
  
   a <- c()
   for(i in 1:length(filenames)){
       tmp <- read.csv(filenames[i])
       a <- union(a,tmp[,"Gene"])
   }
   print(length(a))
  
   genes <- a
   #gene_specific(a,case3,cont3)
   allv <- case3[case3[,"Gene"] %in% genes,]
   allv <- cbind(allv[,c("Gene","Subject_ID")],allv[,setdiff(names(allv),c("Gene","Subject_ID"))])
   write.table(allv,file="filtered_Qiang/mutationsinfo.txt",row.names=FALSE,quote=FALSE,sep="\t")
   igvf <- allv[,c("Chromosome","Position","Subject_ID")]
   write.table(igvf,file="filtered_Qiang/igvfile.txt",row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t")
    
}

gene_specific <- function(genes,case3,cont3){
    source("Faminfo.R")  
    ### populations
    bc.pop <- read.delim("WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    
    ## any LOF and missense mutations in hotspots
    caseid <- unique(case3[,"Subject_ID"])
    path="/ifs/scratch/c2b2/ys_lab/yshen/WENDY/BreastCancer/Regeneron/CaseControl_Filtering/"
    
    #     ## each case with LOF and MIS mutations information 
    #     pang <- read.csv("Panelg_burden.csv")[,1]
    #     for(i in 1:length(caseid)){
    #         famid <- bc.pop[bc.pop[,3]==caseid[i],1]
    #         onefile <- paste(path,"Fam_",famid,".tsv",sep="")
    #         tmp <- read.delim(onefile)
    #         tmp <- tmp[tmp[,"Gene"] %in% pang,]
    #         if(dim(tmp)[1] >0 ){
    #             vs <- paste(tmp[,1],tmp[,2],caseid[i],sep="_")
    #             subs <- vs %in% paste(case3[,1],case3[,2],case3[,"Subject_ID"],sep="_")
    #             tmp <- tmp[subs,]
    #             if(dim(tmp)[1] >0){
    #                 write.csv(tmp,file=paste("filtered_Qiang/",caseid[i],".csv",sep=""),row.names=FALSE)
    #             }
    #         }
    #     }
    
    ### gene based analysis:
    cols <- c("Gene","Chromosome","Position","ID","REF","ALT","VariantFunction","VariantClass","AAchange","AlleleFrequency.ExAC","AlleleFrequency.1KG","AlleleFrequency.ESP","MetaSVM","SIFTprediction","PP2prediction","MAprediction","MTprediction","GERP..","CADDscore","SegmentalDuplication","PredictionSummary","VariantCallQuality","AlternateAlleles","MetaSVMScore","FILTER","INFO")
    ###genes <- unique(read.csv("Panelg_burden.csv")[,1])
    mutalist <- c()
    
    for(i in 1:length(genes)){
        onegid <- unique(case3[case3[,"Gene"] %in% genes[i],"Subject_ID"])
        famid <- bc.pop[bc.pop[,3] %in% onegid,1]
        onelist <- c()
        f1=1
        for(j in 1:length(famid)){
            onefile <- paste(path,"Fam_",famid[j],".tsv",sep="")
            tmp <- read.delim(onefile)
            tmp <- tmp[(tmp[,"Gene"] %in% genes[i]),]
            tmp <- tmp[tmp[,"Position"] %in% case3[case3[,"Gene"]==genes[i],"Position"],]
            
            if(dim(tmp)[1]>0){
                vs <- paste(tmp[,1],tmp[,2],sep="_")
                onecase <- intersect(caseid,bc.pop[bc.pop[,1] %in% famid[j],3])
                col3 <- c(paste("X",onecase,".GT",sep=""),paste("X",onecase,".AD",sep=""),paste("X",onecase,sep=""))
                
                if(f1==1){
                    onelist <- cbind(tmp[,cols],tmp[,col3])
                    f1=f1+1
                }else{
                    vs1 <- paste(onelist[,2],onelist[,3],sep="_")
                    svs1 <- intersect(vs,vs1)
                    svs2 <- setdiff(vs,vs1)
                    if(length(svs1) > 0 ){
                        onelist <- cbind(onelist,"","","")
                        colnames(onelist)[(dim(onelist)[2]-2):dim(onelist)[2]] <- col3
                        onelist[match(svs1,vs1),(dim(onelist)[2]-2):dim(onelist)[2]] <- tmp[match(svs1,vs),col3]
                    }
                    if(length(svs2) > 0 ){
                        n <- length(svs2)
                        ##onelist <- rbind(onelist,rep("",n))
                        for(kj in 1:n) onelist <- rbind(onelist,"")
                        if(!(col3[1] %in% colnames(onelist))){
                            onelist <- cbind(onelist,"","","")
                            colnames(onelist)[(dim(onelist)[2]-2):dim(onelist)[2]] <- col3
                        }
                        onelist[(dim(onelist)[1]-n+1):dim(onelist)[1],c(cols,col3)] <- tmp[match(svs2,vs),c(cols,col3)]
                    }
                } #end if f1==1   
            }
        } # end for j
        write.csv(onelist,file=paste("filtered_Qiang/",genes[i],".csv",sep=""),row.names=FALSE)
        
        for(k in 1:dim(onelist)[1]){
            b <- onelist[k,27:dim(onelist)[2]]
            tmp <- which(b!="")
            sub <- tmp[seq(1,length(tmp),3)]
            sam <- gsub("X","",colnames(onelist)[sub+26])
            sam <- gsub("\\.GT","",sam)
            tmp <- c("","","","")
            for(ki in 1:length(sub)){
                tmp[1] <- paste(tmp[1],sam[ki],":",onelist[k,sub[ki]+26],";",sep="")
                tmp[2] <- paste(tmp[2],sam[ki],":",onelist[k,sub[ki]+27],";",sep="")
                tmp[3] <- paste(tmp[3],sam[ki],":",onelist[k,sub[ki]+28],";",sep="")
                tmp[4] <- paste(tmp[4],sam[ki],",",sep="")
            }
            tmp[4] <- gsub(",$","",tmp[4])
            
            tmp <- c(onelist[k,1:26],tmp)
            names(tmp) <- c(colnames(onelist)[1:26],c("GT","AD","GT:AD:DP:GQ:PL","Subject_IDs"))
            mutalist <- rbind(mutalist,tmp)        
        }
    }
    
    write.table(mutalist,file="filtered_Qiang/genes_mutations.txt",row.names=FALSE,quote=FALSE,sep="\t")
    
    ###genes <- genes[1:10]
    muta10 <- mutalist[mutalist[,"Gene"] %in% genes,c("Chromosome","Position","Subject_IDs")]
    muta30 <- mutalist[mutalist[,"Gene"] %in% genes,]
    write.table(muta30,file="variantsinfo.txt",row.names=FALSE,quote=FALSE,sep="\t")
    #write.table(muta10,file="filtered_Qiang/genes_selected.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t") ## one variant one file
    muta20 <- c()
    for(i in 1:dim(muta10)[1]){
        tmp <- unlist(strsplit(muta10[i,3][[1]],","))
        for(j in 1:length(tmp)){
            tmp1 <- c(muta10[i,1:2],tmp[j])
            muta20 <- rbind(muta20,tmp1)
        }
    }
    write.table(muta20,file="filtered_Qiang/genes_selected.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")  ## one variant in one bam one snapshot
    
}

## hotspots and sample batches
hotspots <- function(){
    hotf <- "cluster_pos3.txt"
    hotT <- read.delim(hotf,sep=";",header=FALSE)
    ###!!!!!!!!!!!!!!!!!
    hotT <- hotT[!duplicated(hotT[,2]),]
    hots <- list()
    for(i in 1:dim(hotT)[1]){
        gene <- hotT[i,2]
        pots <- union(unlist(strsplit(hotT[i,3],"\t")),union(unlist(strsplit(hotT[i,4],"\t")),unlist(strsplit(hotT[i,5],"\t"))))
        pots <- setdiff(pots,"")
        hots[[gene]] <- pots
    }
    save(hots,file="hotspot")
}

hotspots_cosmic <- function(){
  
    ### cbioportal 
#     library("cgdsr")
#     mycgds = CGDS("http://www.cbioportal.org/public-portal/")
#     stul <- getCancerStudies(mycgds)
#     mycl = getCaseLists(mycgds,stul)
#     cbio <- getMutationData(mycgds,caseList=mycl,genes=genes)
    
    ### cosmic 
    cos <- read.delim("hotspots/CosmicMutantExport.tsv",sep="\t")
    mis <- c("Substitution - Missense")
    cos <- cos[cos[,"Mutation.Description"] %in% mis,]
    cutn <- 3
    mgp3 <- names(table(cos[,"Mutation.genome.position"]))[table(cos[,"Mutation.genome.position"])>=cutn]
    cos <- cos[cos[,"Mutation.genome.position"] %in% mgp3,]
    save(cos,file="hotspots/cosmic_missense_3")

    cosS <- cos[,c("Accession.Number","Gene.name","Mutation.AA","Mutation.genome.position")]
    tmp <- gsub("\\D","",cos[,"Mutation.AA"])
    cosS <- cbind(cosS,tmp)
    save(cosS,file="hotspots/cosmic_hotspots")
    
    genes <- unique(cosS[,"Gene.name"])
    hots <- list()
    for(i in 1:length(genes)){
        gene <- genes[i]
        pots <- cosS[cosS[,"Gene.name"]==gene,"tmp"]
        pots <- setdiff(pots,"")
        hots[[gene]] <- pots
    }
    save(hots,file="hotspots/hotf_cos")
    
    
    ## relax hotspot with >=2 mutations in COSMIC ==================
    cos <- read.delim("hotspots/CosmicMutantExport.tsv",sep="\t")
    mis <- c("Substitution - Missense")
    cos <- cos[cos[,"Mutation.Description"] %in% mis,]
    cutn <- 2
    mgp2 <- names(table(cos[,"Mutation.genome.position"]))[table(cos[,"Mutation.genome.position"])>=cutn]
    cos <- cos[cos[,"Mutation.genome.position"] %in% mgp2,]
    save(cos,file="hotspots/cosmic_missense_2")
    
    cosS <- cos[,c("Accession.Number","Gene.name","Mutation.AA","Mutation.genome.position")]
    tmp <- gsub("\\D","",cos[,"Mutation.AA"])
    cosS <- cbind(cosS,tmp)
    save(cosS,file="hotspots/cosmic_hotspots_2")
    
    genes <- unique(cosS[,"Gene.name"])
    hots <- list()
    for(i in 1:length(genes)){
        gene <- genes[i]
        pots <- cosS[cosS[,"Gene.name"]==gene,"tmp"]
        pots <- setdiff(pots,"")
        hots[[gene]] <- pots
    }
    save(hots,file="hotspots/hotf_cos_2")   

}

### single gene double check
WRN_check <- function(){
    load("caselist_8_26")
    onelist <- caselist
    Ecut=0.001
    onelist <- onelist[onelist[,"AlleleFrequency.ExAC"]< Ecut,]
    a <- onelist[onelist[,"Gene"]=="WRN",]
    a <- a[a[,"Position"] %in% c(31007890,30945376,30941288),]
    a <- a[order(a[,"Position"]),]
    qwcsv(a,file="WRNcheck.csv")
    pheno <- read.csv("WES BCFR phenotypic data.csv")
    b <- pheno[match( a[,"Subject_ID"],pheno[,"Subject_ID"]),]
    qwcsv(b,"WRNFam.csv")

}

BRCA_check <-  function(){
    a2 <- read.csv("BRCA2_check.csv")
    a1 <- read.csv("BRCA1_check.csv")
    dups <- unlist(read.table("../WES_Sinai/OverlappedSamples_dups.txt"))
    
    Reps1 <- rep(FALSE,dim(a1)[1])
    Reps2 <- rep(FALSE,dim(a2)[2])
    
    for(i in 1:dim(a1)[1]){
        Reps1[i] <- any(grepl(a1[i,"Subject_ID"],dups))
    }
    
    for(i in 1:dim(a2)[1]){
        Reps2[i] <- any(grepl(a2[i,"Subject_ID"],dups))
    }    

}

### slient mutation number check
vcf_check <- function(){
    
    load("caselist_9_15")
    caselist <- onlycase
    bc.pop <- read.delim("WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    Jp <-  bc.pop[bc.pop[,4] %in% "J",3]
    Hp <-  bc.pop[bc.pop[,4] %in% "H",3]
    unids <- setdiff(unique(caselist[,"Subject_ID"]),c(Jp,Hp))
    gvcf <- list.files(path="/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/gvcfs",pattern=".gvcf.gz$")
    tmpn <- rep("",length(unids))
    for(i in 1:length(unids)){
        tmpn[i] <- gvcf[grepl(unids[i],gvcf)]
    }
    ## all samples from Yale samples

    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample_M_240041_005_005
    cols <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","Sample")
    # 240041 
    a <- read.delim("./single_check/240041_1.vcf",comment.char="#",sep="\t",header=FALSE) # 980 unique slient variants
    # double check vcf files
    colnames(a) <- cols
    
    
}

## sample checked 
sample_check <- function(){

    load("SampleBatches")
    bc.pop <- read.delim("WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    Jp <-  bc.pop[bc.pop[,4] %in% "J",3]
    Hp <-  bc.pop[bc.pop[,4] %in% "H",3]
    
    load("caselist_9_15")
    caseid <- unique(onlycase[,"Subject_ID"])
    
    YaleSam <- bats[bats[,3]=="Yale",2]
    caseYa <- intersect(caseid,YaleSam)
    
    caseY1 <- intersect(c(Jp,Hp),caseYa)
    caseY2 <- intersect(setdiff(caseid,c(Jp,Hp)),caseYa)
    others <- setdiff(caseid,caseYa)
    
    sig=FALSE
    load("casesy_9_15")
    load("contsy_9_15")
    
    a <- table(casesy[,"Subject_ID"])[caseY1]
    b <- table(casesy[,"Subject_ID"])[caseY2]
    c <- table(casesy[,"Subject_ID"])[others]
    
    median(a)
    median(b)
    median(c)
    
    casesy <- filter_variant(casesy,sig)
    contsy <- filter_variant(contsy,sig)
    casesy <- casesy[casesy[,"Variantfiltering"],]
    contsy <- contsy[contsy[,"Variantfiltering"],]
    
    a <- table(casesy[,"Subject_ID"])[caseY1]
    b <- table(casesy[,"Subject_ID"])[caseY2]
    c <- table(casesy[,"Subject_ID"])[others]
    
    median(a)
    median(b)
    median(c)

}

sampleBatch <- function(){
    
    bams <- unlist(read.table("bam.txt"))
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
    qwcsv(bats,file="Sample_batches.csv")
    save(bats,file="SampleBatches")
    
    ### batches for case and control
    load("cont3")
    load("case3")  
    
    contid <- unique(cont3[,"SubID"])
    table(bats[match(contid,bats[,2]),3])   
    
    caseid <- unique(case3[,"Subject_ID"])
    table(bats[match(caseid,bats[,2]),3])
    ## batches for each gene or mutation
    oneg="PABPC3"
    caseid <- unique(case3[case3[,"Gene"]==oneg,"Subject_ID"])
    table(bats[match(caseid,bats[,2]),3])
    contid <- unique(cont3[cont3[,"Gene"]==oneg,"SubID"])
    table(bats[match(contid,bats[,2]),3])   
    
    casevl <- paste(case3[,1],case3[,2],sep="_")
    contvl <- paste(cont3[,1],cont3[,2],sep="_")
    
    onev <- "13_25671272"
    caseid <- unique(case3[casevl==onev,"Subject_ID"])
    table(bats[match(caseid,bats[,2]),3])
    contid <- unique(cont3[contvl==onev,"SubID"])
    table(bats[match(contid,bats[,2]),3])   
    
}

