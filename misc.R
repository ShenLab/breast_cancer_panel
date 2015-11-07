getVariantlist <- function(path,IDfile,namestr=".tsv",savefile){
## get variant lists for a set of samples
    print_log(paste("getVariantlist function is running ...", date(),sep=" "))
    
    IDs <- unlist(read.table(IDfile))
    samf <- paste(IDs,namestr,sep="")
    files <- list.files(path=path,pattern=".tsv$")
    samf <- intersect(samf,files)
    
    print_log(paste("getVariantlist: All samples have variant files in the given path:", all(samf %in% files),sep=" "))
    print_log(paste("getVariantlist: The number of subjects in this variant list are:", length(samf), sep=" "))
    
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
    
    save(onelist,file=savefile)
    print_log(paste("getVariantlist function done!", date(),sep=" "))
    #onelist
}

variant_filtering <- function(onelist,mis,Ecut=0.01,segd=0.95,pp2=TRUE,sig=FALSE,hotf=""){
    ### onelist: variant list
    ### Ecut: ExAC frequency cut off
    ### segd: segment duplication score cut off in VCF
    ### pp2: whether use PolyPhen2 to predict damaging missense or not
    ### sig: whether only consider the singleton variants or not
    ### hotf: filtering missense with hotspot information
    print_log(paste("variant_filtering function is running ...", date(),sep=" "))
    
    filters <- c("filtered","ExACfreq","VCFPASS","noneSegmentalDup","meta-SVM_PP2","singleton","hotspot")
    filS <- matrix(FALSE,dim(onelist)[1],length(filters))
    colnames(filS) <- filters
    onelist <- cbind(onelist,filS)
    
    # EXAC cut <- 0.001 ## rare variants
    onelist[is.na(onelist[,"AlleleFrequency.ExAC"]),"AlleleFrequency.ExAC"] <- 0
    onelist[onelist[,"AlleleFrequency.ExAC"]==".","AlleleFrequency.ExAC"] <- 0
    onelist[as.numeric(onelist[,"AlleleFrequency.ExAC"])< Ecut,"ExACfreq"] <- TRUE
    
    ### filtered more details in VCF
    badvars <- c('QD_Bad_SNP','FS_Bad_SNP','FS_Mid_SNP;QD_Mid_SNP','LowQuality','LowQD_Indel','LowQuality','VQSRTrancheSNP99.90to100.00','VQSRTrancheINDEL99.90to100.00')
    subs1 <- sapply(1:dim(onelist)[1],function(i){
        tmp <- unlist(strsplit(onelist[i,"FILTER"],";"))
        a1 <- intersect(tmp,badvars)
        a2 <- grepl("FS_Mid_SNP;QD_Mid_SNP",onelist[i,"FILTER"])
        (length(a1) > 0) | a2
    })
    onelist[!subs1,"VCFPASS"] <- TRUE
    
    subs2 <- sapply(1:dim(onelist)[1], function(i) {
        if(onelist[i,"SegmentalDuplication"] == "none"){ TRUE;
        }else{ tmp <- unlist(strsplit(onelist[i,"SegmentalDuplication"],","))[1]
               as.numeric(unlist(strsplit(tmp,":"))[2]) < segd
        }
    })
    onelist[subs2,"noneSegmentalDup"] <- TRUE
    
    ### missense predicted by meta-SVM and PP2
    onelist[,"meta-SVM_PP2"] <- TRUE
    subs2 <- onelist[,"VariantClass"] %in% mis
    if(pp2){subs1 <- onelist[,"PP2prediction"]=="D" | onelist[,"MetaSVM"]=="D";}else{subs1 <- onelist[,"MetaSVM"]=="D";}   
    subs3 <- nchar(onelist[,"REF"]) != nchar(onelist[,"ALT"]) ### indels
    subs2 <- subs2 & !subs3
    onelist[subs2 & !subs1,"meta-SVM_PP2"] <- FALSE
    
    onelist[,"singleton"] <- TRUE
    if(sig){
        vs <- paste(onelist[,1],onelist[,2],onelist[,4],onelist[,5],sep="_") ####!!!!
        vsdup <- vs[duplicated(vs)]
        subs1 <- vs %in% vsdup
        onelist[subs1, "singleton"] <- FALSE
    }
    
    onelist[,"hotspot"] <- TRUE
    if(hotf!=""){
        missub <- which( (onelist[,"VariantClass"] %in% mis) & (nchar(onelist[,"REF"]) == nchar(onelist[,"ALT"])) ) ### only missense mutations
        hotspot <- read.table(hotf)
        hotspot[,2] <- paste(":",hotspot[,2],":",sep="")
        subhot <- sapply(1:length(missub), function(i) {
            if(onelist[missub[i],"Gene"] %in% hotspot[,1]){
                tmp <- unlist(strsplit(onelist[missub[i],"AAchange"],":"))[5];
                if(grepl("_",tmp) | grepl("del",tmp) | grepl("ins",tmp)){ TRUE; 
                }else{ 
                    grepl(paste(":",gsub("\\D","",tmp),":",sep=""),hotspot[hotspot[,1]==onelist[missub[i],"Gene"],2])
                }
            }else{
                FALSE; ###
            }
            })
        onelist[missub[!subhot],"hotspot"] <- FALSE
    }
    
    onelist[,"filtered"] <- rowSums(onelist[,filters[-1]])==(length(filters)-1)
    
    ##======================print information to logFile======================================
    print_log(paste("variant_filtering parameters: Allele Frequency ExAC ", Ecut,sep=" "))
    print_log(paste("variant_filtering parameters: Segment duplication score ", segd,sep=" "))
    print_log(paste("variant_filtering parameters: PolyPhen2 used ", pp2,sep=" "))
    print_log(paste("variant_filtering parameters: singleton variant only ",sig,sep=" "))
    print_log(paste("variant_filtering parameters: variant filtered by hotspots ", hotf!="", sep=" "))
    print_log(paste("variant_filtering function is done!", date(),sep=" "))
    ##=========================================================================================
    onelist
}

burden_test <- function(caselist,contlist,testset,testtype,flag,indel=FALSE){
## variant lists: caselist and contlist
## testset: test gene sets or variants 
## testtype: (missense, LOF and indel)
## flag: 1, gene/variant set; 2, single gene test; 3, single variant test
    print_log(paste("burden_test function is running ...", date(),sep=" ")) 

    n.case <- length(caselist[,"Subject_ID"])
    n.cont <- length(contlist[,"Subject_ID"])
    if(testset=="") testset <- union(caselist[,"Gene"],contlist[,"Gene"])
    if(testtype=="") testtype <- unique(caselist[,"VariantClass"])
    caselist <- caselist[caselist[,"VariantClass"] %in% testtype & caselist[,"Gene"] %in% testset, ]
    contlist <- contlist[contlist[,"VariantClass"] %in% testtype & contlist[,"Gene"] %in% testset, ]
    if(indel){
        caselist <- caselist[nchar(caselist[,"REF"]) != nchar(caselist[,"ALT"]), ]
        contlist <- contlist[nchar(contlist[,"REF"]) != nchar(contlist[,"ALT"]), ]
    }
    casevars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_") 
    contvars <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
        
    if(flag == 1){
        a <- dim(caselist)[1]
        b <- dim(contlist)[1]
        oneTable <- matrix(c(0,0,a,b,n.case,n.cont, (a/n.case)/(b/n.cont), ifelse( (a+b)>0, binom.test(a,a+b,n.case/(n.case+n.cont))$p.value,1)),nrow=1,ncol=8)
    }else if(flag == 2){
        genes <- union(caselist[,"Gene"],contlist[,"Gene"])
        oneTable <- sapply(genes,function(gene){
            a <- sum(caselist[,"Gene"] %in% gene)
            b <- sum(contlist[,"Gene"] %in% gene)
            c(gene,length(union(casevars[caselist[,"Gene"] %in% gene],contvars[contlist[,"Gene"] %in% gene])),a,b,n.case,n.cont,(a/n.case)/(b/n.cont),binom.test(a,a+b,n.case/(n.case+n.cont))$p.value)
        })
        oneTable <- t(oneTable)
    }else if(flag == 3){
        vars <- union(casevars,contvars)
        oneTable <- sapply(vars,function(onevar){
            a <- sum(casevars %in% onevar)
            b <- sum(contvars %in% onevar)
            tmpg <- ifelse(onevar %in% casevars, caselist[which(casevars==onevar)[1],"Gene"], contlist[which(contvars==onevar)[1],"Gene"]) 
            c(tmpg,onevar,a,b,n.case,n.cont,(a/n.case)/(b/n.cont),binom.test(a,a+b,n.case/(n.case+n.cont))$p.value)
        })
        oneTable <- t(oneTable)
    }
    cols <- c("Gene","Variant","#in_case","#in_cont","n.case","n.cont","Folds","Pvalue")
    colnames(oneTable) <- cols
    
    ##======================print information to logFile======================================
    print_log(paste("burden_test: there are totally ",dim(oneTable)[1]," tests",sep=""))
    print_log(paste("burden_test function is done!", date(),sep=" "))
    ##========================================================================================
    oneTable
}

print_log <- function(printstr){
    
    ## write log file for each run
    logFile = "misc_run.log"
    cat(printstr, file=logFile, append=TRUE, sep = "\n")

}
