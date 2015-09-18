## family information and population information
Faminfo <- function(){

    ### populations
    bc.pop <- read.delim("PCA/WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    
    ###

}

case_cont <- function(){
    source("Faminfo_9_4_without_log.R")
    load("caselist_8_26")
    load("contlist_8_6")
    
    caselist <- filter_variant(caselist)
    contlist <- filter_variant(contlist)
    
    ### populations
    bc.pop <- read.delim("WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
    ###table(bc.pop[,4])
    
    ## LOF mutations analysis 
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss")                 
    case1 <- caselist[caselist[,"VariantClass"] %in% lof,]    
    cont1 <- contlist[contlist[,"VariantClass"] %in% lof,]
    aa = geneburden_LOF(case1,cont1,bc.pop)
    print(aa)
    
    ## missense mutations in cosmic hotspots analysis
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    case2 <- caselist[caselist[,"VariantClass"] %in% mis,]
    cont2 <- contlist[contlist[,"VariantClass"] %in% mis,]
    hotf <- "hotspots/hotf_cos"
    bb = geneburden_MIS(case2,cont2,bc.pop,hotf)
    print(bb[[1]])
    print(bb[[2]])
    
    ## any LOF and missense mutations in hotspots
    case3 <- rbind(case1,bb$caselist)
    cont3 <- rbind(cont1,bb$contlist)
    
    save(case3,file="case3")
    save(cont3,file="cont3")
    
    ## missense mutations analysis
    hotf <- ""
    cc = geneburden_MIS(case2,cont2,bc.pop,hotf,hotk=2)

} 

filter_variant <- function(onelist){
    
    filters <- c("ExACfreq","VCFPASS","noneSegmentalDup","meta-SVM_PP2","GTEXexp","singletonLOF")
    # EXAC cut <- 0.001 ## rare variants
    Ecut=0.001
    onelist <- onelist[onelist[,"AlleleFrequency.ExAC"]< Ecut,]
    
    ### filtered more details
    subs1 <- onelist[,"FILTER"] == "PASS"
    subs2 <- onelist[,"SegmentalDuplication"] == "none"
    onelist <- onelist[subs1 & subs2,]
    ### missense predicted by meta-SVM and PP2
    mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
    subs1 <- !(onelist[,"MetaSVM"] %in% c("D",".")) & !(onelist[,"PP2prediction"] %in% c("D",".")) ##!!!! 
    subs2 <- onelist[,"VariantClass"] %in% mis
    onelist <- onelist[!(subs1 & subs2), ]
    ### expressed genes in GTEX
    expg <- unlist(read.table("../geneexpression/GTEx/expg.txt"))
    onelist <- onelist[onelist[,"Gene"] %in% expg,]
    
    ### LOF in singleton variants only !!!!
    lof <- c("frameshiftdeletion","frameshiftinsertion","stopgain","stoploss")
    vs <- paste(onelist[,1],onelist[,2],sep="_")
    vsdup <- vs[duplicated(vs)]
    subs1 <- vs %in% vsdup
    subs2 <- onelist[,"VariantClass"] %in% lof
    onelist <- onelist[!(subs1 & subs2), ]
    
    ### no dbsnp
    ### onelist <- onelist[onelist[,"ID"]==".",]
    
    #### write log file for each version
    con <- file("variant_filtering.log","w");
    writeLines(paste("Filtering log: ",date(),sep=""),con)
    writeLines(paste("ExAC frequency: ",Ecut," (ExAc",Ecut,")",sep=""),con)
    writeLines(paste("VCF FILTER: PASS"," (VCFPASS)",sep=""),con)
    writeLines(paste("VCF SegmentalDuplication: none"," (noneSegmentalDup)",sep=""),con)
    writeLines(paste("Meta-SVM or PP2 prediction: D or . for missense mutations"," (meta-SVM_PP2)",sep=""),con)
    writeLines(paste("GTEX expressed genes "," (GTEXexp)",sep=""),con)
    writeLines(paste("Singleton variants for LOF mutations"," (singletonLOF)",sep=""),con)
    close(con)
    
    onelist
}

geneset_burden <- function(){
    ## gene set Panel burden 
    ## cols <- c("Gene","Uni","n.case","case","control","odds_ratio","p_value","case_J","control_J","J_odds_ratio","J_p_value","case_H","control_H","H_odds_ratio","H_p_value","TYPE")
    ts <- unlist(read.table("hotspots/TS_filtered.txt")) ## tumor suppresspor
    drs <- unlist(read.table("hotspots/Driver_filtered.txt")) ## cancer drivers
    
    ### tumor suppressors in rare LGD variants    
    LOFf <- "GeneBurden_LOF.csv"
    write_panel(ts,LOFf, "LGD_TS_burden.csv")
    
    ### cancer dirver genes in rare LGD variants
    LOFf <- "GeneBurden_LOF.csv"
    write_panel(drs,LOFf, "LGD_driver_burden.csv")
    
    ### cancer driver genes in rare missense hotspots
    MISf0 <- "GeneBurden_MIS_hotspot.csv"
    write_panel(drs,MISf0, "MIShotspot_driver_burden.csv")
    
    ### tumor suppressor genes in rare missense hotspots
    MISf0 <- "GeneBurden_MIS_hotspot.csv"
    write_panel(ts,MISf0, "MIShotspot_TS_burden.csv")    
    
    ### cancer driver genes in rare missense with all missense
    MISf <- "GeneBurden_MIS.csv"
    write_panel(drs,MISf, "MIS_driver_burden.csv")

}

filter_gene <- function(){
    ## tumor suppressors 
    ts1 <- unlist(read.table("../genelist/Genelist2.txt"))
    ts2  <- unlist(read.table("../genelist/Tumor_supressor/TS_Vanderbilt.txt")) ## tumor supressors
    ts <- union(ts1,ts2)
    
    ### COSMIC LOF percentage >= 15%
    #cos <- read.delim("hotspots/CosmicMutantExport.tsv",sep="\t")
    #cos1 <- cos[cos[,"Gene.name"] %in% ts,]
    #save(cos1,file="hotspots/TScheck")
    load("hotspots/TScheck")
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
    aa <- read.csv(burdenf)
    cc <- aa[aa[,1] %in% genes,]
    cc <- cc[order(cc[,"p_value"]),]
    write.csv(cc,file=wf,row.names=FALSE)
}

variant_type_burden <- function(){
    
    n.case <- c(345,216,94)
    n.cont <- c(114,59,55)
    
    ## all class: LOF and missense in hotspots
    loff <- "GeneBurden_LOF.csv"
    misf <- "GeneBurden_MIS_hotspot.csv"
    cc <- write_variant(loff,misf,n.case,n.cont)
    print(cc)
    
    ## LOF in tumor suppressors
    loff <- "LGD_TS_burden.csv"
    ## missense in tumor suppressors
    misf <- "MIShotspot_TS_burden.csv"
    cc <- write_variant(loff,misf,n.case,n.cont)
    print(cc)
    
    ## LOF in cancer dirver
    loff <- "LGD_driver_burden.csv"
    ## missense in cancer dirver
    misf <- "MIShotspot_driver_burden.csv"
    cc <- write_variant(loff,misf,n.case,n.cont)
    print(cc)
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
            cc[i,4] <- binom.test(cc[i,1],cc[i,1]+cc[i,2],n.case[j]/(n.case[j]+n.cont[j]),alternative="greater")$p.value
        }
    }
    
    cc
}

check_genes <- function(){
    a <- read.csv("LGD_driver_burden.csv")[,1]
    b <- read.csv("LGD_TS_burden.csv")[,1]
    c <- read.csv("MIShotspot_driver_burden.csv")[,1]
    d <- read.csv("MIShotspot_TS_burden.csv")[,1]
    
    a <- union(a,union(b,c))
    a <- union(a,d)
    print(length(a))
    gene_specific(a)
    
}

gene_specific <- function(genes){
    source("Faminfo.R")  
    ### populations
    bc.pop <- read.delim("WES_BCFR_phenotypic_data-19062015.txt")[,1:5]
    bc.pop[,4] <- paste(bc.pop[,4], bc.pop[,5], sep="")
    bc.pop <- bc.pop[,-5]
  
    ## any LOF and missense mutations in hotspots
    load("case3")
    load("cont3")
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
    write.table(muta10,file="filtered_Qiang/genes_selected.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    
}

## gene level burden test
geneburden_LOF <- function(caselist,contlist,bc.pop){
    
    cols <- c("Gene","UniLOF","n.case","case_LOF","control_LOF","odds_ratio","p_value","case_J_LOF","control_J_LOF","J_odds_ratio","J_p_value","case_H_LOF","control_H_LOF","H_odds_ratio","H_p_value")
    result <- geneburden(caselist,contlist,bc.pop,cols)
    Sta <- result$Sta
    n.case <- result$n.case
    n.cont <- result$n.cont
    write.csv(Sta,file="GeneBurden_LOF.csv",row.names=FALSE)
    
    list(n.case=n.case,n.cont=n.cont)
}

geneburden_MIS <- function(caselist,contlist,bc.pop,hotf,hotk=1){
    
    if(hotk==1){
        load(hotf)
        hotg <- names(hots)
        caselist <- caselist[caselist[,"Gene"] %in% hotg,]
        contlist <- contlist[contlist[,"Gene"] %in% hotg,]
        print(dim(caselist))
        print(dim(contlist))
        ### only missense mutations in hotspots and indels
        p.case <- sapply(1:dim(caselist)[1], function(i) {
            tmp <- unlist(strsplit(caselist[i,"AAchange"],":"))[5];
            if(grepl("_",tmp) | grepl("del",tmp) | grepl("ins",tmp)){ TRUE; }else{ gsub("\\D","",tmp);}
        })
        s.case <- sapply(1:length(p.case),function(i){
            k <- which(hotg==caselist[i,"Gene"])
            p.case[i] %in% hots[[k]]
        })
        caselist <- caselist[s.case | (p.case==TRUE),]
        
        p.cont <- sapply(1:dim(contlist)[1], function(i) {
            tmp <- unlist(strsplit(contlist[i,"AAchange"],":"))[5];
            if(grepl("_",tmp) | grepl("del",tmp) | grepl("ins",tmp)){ TRUE; }else{ gsub("\\D","",tmp);}
        })
        s.cont <- sapply(1:length(p.cont),function(i){
            k <- which(hotg==contlist[i,"Gene"])
            p.cont[i] %in% hots[[k]]
        })
        contlist <- contlist[s.cont | (p.cont==TRUE),]
        print(dim(caselist))
        print(dim(contlist))    
    }
    
    cols <- c("Gene","UniMIS","n.case","case_MIS","control_MIS","odds_ratio","p_value","case_J_MIS","control_J_MIS","J_odds_ratio","J_p_value","case_H_MIS","control_H_MIS","H_odds_ratio","H_p_value")
    result <- geneburden(caselist,contlist,bc.pop,cols)
    Sta <- result$Sta
    n.case <- result$n.case
    n.cont <- result$n.cont
    
    if(hotk==1){
        write.csv(Sta,file="GeneBurden_MIS_hotspot.csv",row.names=FALSE)
    }else{
        write.csv(Sta,file="GeneBurden_MIS.csv",row.names=FALSE)
    }
    list(n.case=n.case,n.cont=n.cont,caselist=caselist,contlist=contlist)
    
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
    n.case <- c(length(caseid),length(intersect(caseid,Jp)),length(intersect(caseid,Hp)))
    n.cont <- c(length(contid),length(intersect(contid,Jp)),length(intersect(contid,Hp)))
    
    Sta[is.na(Sta)] <- 0
    for(i in 1:3){
        ntmp <- sapply(1:length(genes), function(j){
            tmp <- fisher.test(matrix(c(Sta[j,4*i],Sta[j,4*i+1],n.case[i],n.cont[i]),2,2),alternative="greater");
            c(tmp$estimate,tmp$p.value)
        })
        Sta[,c(4*i+2,4*i+3)] <- t(ntmp)
    }
    
    Sta[,1] <- genes
    Sta <- Sta[order(-as.numeric(Sta[,2])),]
    
    list(Sta=Sta,n.case=n.case,n.cont=n.cont)
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
    
}

del <- function(){
# a <- read.csv("filtered_Qiang/PABPC3.csv")
# n <- (dim(a)[2]-26)/3
# m <- dim(a)[1]
# nC <- rep(0,m)
# for(i in 1:m){
#     b <- a[i,27:dim(a)[2]]
#     nC[i] <- n - sum(b=="")/3
# }
# nC
    
}