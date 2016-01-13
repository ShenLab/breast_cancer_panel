phenoinfo <- function(){
    phenofile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/WES BCFR phenotypic data.csv"
    variantpath <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"
    
    files <- list.files(path=variantpath,pattern=".tsv$")
    subjects <- gsub(".AllVariants.tsv","",files)
    pheno <- read.csv(phenofile)
    pheno[pheno[,3]=="222357, 222966",3] <- "222357" 
    pheno <- pheno[pheno[,3] %in% subjects,]
    
    pheno
}

excluded_samples <- function(){
    
    outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"  ## outlier samples
    BRCA1_2pathogenicfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/BRCA1_2.txt" ## BRCA1/2 pathogenic samples file
    
    outliers <- read.delim(outlierfile,header=FALSE)[,2]
    tmp <- read.delim(BRCA1_2pathogenicfile)
    pathogenic_sample <- tmp[tmp[,1] %in% c("likely pathogenic","pathogenic"),"Subject_ID"]
    
    exSamples <- union(outliers,pathogenic_sample)
    exSamples
}

getVariantlist <- function(path,IDfile,namestr=".tsv",savefile){
## get variant lists for a set of samples
    print_log(paste("getVariantlist function is running ...", date(),sep=" "))
    
    IDs <- unlist(read.table(IDfile))
    samf <- paste(IDs,namestr,sep="")
    files <- list.files(path=path,pattern=".tsv$")
    samf <- files[match(samf,files)]
	
    if(any(is.na(samf))){
	subs <- sapply(1:length(IDs), function(i) which(grepl(IDs[i], files)) )
	samf <- files[subs]
    }
    
    print_log(paste("getVariantlist: All samples have variant files in the given path:", all(samf %in% files),sep=" "))
    print_log(paste("getVariantlist: The number of subjects in this variant list are:", length(samf), sep=" "))
    
    onelist <- c()
    for(i in 1:length(samf)){
        tmp <- paste(path,samf[i],sep="")
        oner <- read.delim(tmp,check.names=FALSE)
        #subj <- gsub(namestr,"",samf[i])
    	subj <- IDs[i]
	oner <- cbind(oner,subj)
        cols <- colnames(oner)
        colsub <- c(which(grepl(toupper(subj),toupper(cols))), dim(oner)[2])
        colnames(oner)[colsub] <- c("GT","AD","Subject_INFO","Subject_ID")
        onelist <- rbind(onelist,oner)
    }
    
    save(onelist,file=savefile)
    print_log(paste("getVariantlist function done!", date(),sep=" "))
    #onelist
}

getindexcase <- function(phenofile,agem="max"){
    print_log("\n\n")
    print_log(paste("getindexcase function is running ...", date(),sep=" "))
    pheno <- read.csv(phenofile)
    ## one case in one family
    famid <- unique(pheno[,1])
    subs <- rep(FALSE,dim(pheno)[1])
    ages <- sapply(1:dim(pheno)[1], function(i) 114 -  as.numeric(unlist(strsplit(pheno[i,"BIRTHDT"],"/"))[3]) )
    print_log(paste("getindexcase: All subjects' ages are avaiable: ", !any(is.na(ages)),sep=" "))
    for(i in 1:length(famid)){
        casesub <- which(pheno[,"BreastCancer"] == "Yes" & pheno[,1]%in% famid[i])
        if(agem=="max"){
            onesub <- which.max(ages[casesub]) ## older cases
        }else{
            onesub <- which.min(ages[casesub]) ## younger cases
        }
        subs[casesub[onesub]] <- TRUE
    }
    indexcases <- pheno[subs,3] ## subject_IDs
    print_log(paste("getindexcase: The number of index cases is ", length(unique(indexcases)),sep=" "))
    print_log(paste("getindexcase function is done!", date(),sep=" "))
    indexcases
}

VariantSta <- function(casef,contf,cases,controls,outputpath){
    ## print to log
    print_log(paste("VariantSta function is running ...", date(),sep=" "))
    print_log(paste("VariantSta: the number of cases is: ", length(unique(cases)),sep=" "))
    print_log(paste("VariantSta: the number is controls is: ", length(unique(controls)),sep=" "))
    
    ## running
    if(!file.exists(outputpath)){ dir.create(outputpath, showWarnings = TRUE, recursive = FALSE);}
    varT=c("ALL","SNVs","InDels","Frameshift","Missense","Silent")
    varTa <- VariantStaonef(casef,cases,varT,"case")
    varTa1 <- VariantStaonef(contf,controls,varT,"control")
    Ta <- rbind(varTa,varTa1)
    n <- dim(Ta)[2]
    density_plots(Ta[,3:n],Ta[,1],outputpath,VarT=varT)
    
    ## done
    print_log(paste("VariantSta function is done!", date(),sep=" "))
}

VariantStaonef <- function(casef,cases,varT,group){

    tmp <- read.delim(casef,sep="\t",header=FALSE)
    idsub <- which(tmp[,1]=="Samples:")[1]
    if(all(cases %in% tmp[idsub,])){
        subs <- match(cases,tmp[idsub,])
    }else{
        subs <- sapply(1:length(cases),function(i) which(grepl(cases[i],tmp[idsub,])))
    }
    tsub <- match(varT[2:length(varT)],tmp[,1])
    varTa <- t(tmp[tsub,subs])
    atmp <- as.matrix(tmp[tsub[1:2],subs])
    mode(atmp) <- "numeric"
    varTa <- cbind(colSums(atmp),varTa)
    varTa <- cbind(group,cases,varTa)
    colnames(varTa) <- c("Group","Subject_ID",varT)
    
    varTa
}

density_plots <- function(oneVar,g,outputpath,VarT=c("LOF","MIS","indel","synonymous","unknown","ALL")){
    if(!file.exists(outputpath)){
        dir.create(outputpath, showWarnings = TRUE, recursive = FALSE)
    }
    mode(oneVar) <- "numeric"
    
    for(i in 1:dim(oneVar)[2]){
        pdf(paste(outputpath,VarT[i],"_Dis.pdf",sep=""),height=8,width=10)
        par(mai=c(2,1,1,1))
        x1 <- density(oneVar[g=="case",i])
        x2 <- density(oneVar[g=="control",i])
        plot(x1,col=1,type="l",main=paste(VarT[i]," : case-control distribution",sep=""),xlab=paste("Number of ",VarT[i]," variants",sep=""),xlim=c(min(oneVar[,i]),max(oneVar[,i])),ylim=c(0,max(c(x1$y,x2$y))),cex.lab=1.5,cex.axis=1.5)
        lines(x2,col=2,type="l")
        abline(v=mean(oneVar[g=="case",i]),col=1,lty=2)
	abline(v=mean(oneVar[g=="control",i]),col=2,lty=2)
	legend("topright",legend=c("case","control"),lty=rep(1,2),lwd=2,col=1:2,cex=1.5)
        dev.off()
    }
}

variant_filtering <- function(onelist,mis,Ecut=0.01,segd=0.95,pp2=TRUE,hotf="",alleleFrefile=NULL,popcut=0.05){
    ### onelist: variant list
    ### Ecut: ExAC frequency cut off
    ### segd: segment duplication score cut off in VCF
    ### pp2: whether use PolyPhen2 to predict damaging missense or not
    ### hotf: filtering missense with hotspot information
    print_log(paste("variant_filtering function is running ...", date(),sep=" "))
    
    filters <- c("filtered","ExACfreq","VCFPASS","noneSegmentalDup","meta-SVM_PP2","hotspot","alleleFre")
    filS <- matrix(FALSE,dim(onelist)[1],length(filters))
    colnames(filS) <- filters
    onelist <- cbind(onelist,filS)
    
    # EXAC cut <- 0.001 ## rare variants
    onelist[is.na(onelist[,"AlleleFrequency.ExAC"]),"AlleleFrequency.ExAC"] <- 0
    onelist[onelist[,"AlleleFrequency.ExAC"]==".","AlleleFrequency.ExAC"] <- 0
    onelist[as.numeric(onelist[,"AlleleFrequency.ExAC"])< Ecut,"ExACfreq"] <- TRUE
    
    onelist[,"alleleFre"] <- TRUE
    if(!is.null(alleleFrefile)){
        alleleFres <- read.delim(alleleFrefile)
        popvars <- paste(alleleFres[,"CHROM"],alleleFres[,"POS"],alleleFres[,"ALLELE"],sep="_")
        allvars <- paste(onelist[,"Chromosome"],onelist[,"Position"],onelist[,"ALT"],sep="_")
        igvars <- intersect(popvars,allvars)
        filteredvars <- igvars[alleleFres[match(igvars,popvars),"ALLELE_FREQ"] >= popcut]
        onelist[allvars %in% filteredvars,"alleleFre"] <- FALSE
    }
    
    ### filtered more details in VCF
    badvars <- c('QD_Bad_SNP','FS_Bad_SNP','FS_Mid_SNP;QD_Mid_SNP','LowQuality','LowQD_Indel','LowQuality','VQSRTrancheSNP99.90to100.00','VQSRTrancheINDEL99.90to100.00','VQSRTrancheINDEL99.00to99.90')
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
    #subs3 <- nchar(onelist[,"REF"]) != nchar(onelist[,"ALT"]) ### indels
    #subs2 <- subs2 & !subs3
    onelist[subs2 & !subs1,"meta-SVM_PP2"] <- FALSE
        
    onelist[,"hotspot"] <- TRUE
    if(hotf!=""){
        missub <- which(onelist[,"VariantClass"] %in% mis)  ### only missense mutations
        hotspot <- read.delim(hotf,header=FALSE)
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
    print_log(paste("variant_filtering parameters: variant filtered by hotspots ", hotf!="", sep=" "))
    print_log(paste("variant_filtering parameters: population frequency filter used ", !is.null(alleleFrefile), sep=" "))
    if(!is.null(alleleFrefile)) print_log(paste("variant_filtering parameters: population frequency cutoff ", popcut, sep=" "))
    print_log(paste("variant_filtering function is done!", date(),sep=" "))
    ##=========================================================================================
    onelist
}

singleton_variants <- function(casevars,contvars){
    vs <- c(casevars,contvars)
    singletons <- setdiff(vs,vs[duplicated(vs)])
    singletons
}

burden_test <- function(caselist,contlist,testset=NULL,testtype=NULL,flag,sig=FALSE,coe=1){
## variant lists: caselist and contlist
## testset: test gene sets or variants 
## testtype: (missense, LOF and indel)
## flag: 1, gene/variant set; 2, single gene test; 3, single variant test
## sig: singleton variants or not
    
    print_log(paste("burden_test function is running ...", date(),sep=" ")) 

    n.case <- length(unique(caselist[,"Subject_ID"]))
    n.cont <- length(unique(contlist[,"Subject_ID"]))
    if(is.null(testset)) testset <- union(caselist[,"Gene"],contlist[,"Gene"])
    if(is.null(testtype)) testtype <- unique(caselist[,"VariantClass"])
    caselist <- caselist[caselist[,"VariantClass"] %in% testtype & caselist[,"Gene"] %in% testset, ]
    contlist <- contlist[contlist[,"VariantClass"] %in% testtype & contlist[,"Gene"] %in% testset, ]
    casevars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_") 
    contvars <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
    
    if(flag==1){
        if(flag == 1){
            a <- dim(caselist)[1]
            b <- dim(contlist)[1]
            b <- floor(b * coe)
            oneTable <- matrix(c(0,0,a,b,n.case,n.cont, (a/n.case)/(b/n.cont), ifelse( (a+b)>0, binom.test(a,a+b,n.case/(n.case+n.cont))$p.value,1)),nrow=1,ncol=8)
        }else if(flag == 2){
            genes <- union(caselist[,"Gene"],contlist[,"Gene"])
            oneTable <- sapply(genes,function(gene){
                a <- sum(caselist[,"Gene"] %in% gene)
                b <- sum(contlist[,"Gene"] %in% gene)
                b <- floor(b * coe)
                c(gene,length(union(casevars[caselist[,"Gene"] %in% gene],contvars[contlist[,"Gene"] %in% gene])),a,b,n.case,n.cont,(a/n.case)/(b/n.cont),binom.test(a,a+b,n.case/(n.case+n.cont))$p.value)
            })
            oneTable <- t(oneTable)
        }else if(flag == 3){
            vars <- unique(casevars)
            oneTable <- sapply(vars,function(onevar){
                a <- sum(casevars %in% onevar)
                b <- sum(contvars %in% onevar)
                tmpg <- ifelse(onevar %in% casevars, caselist[which(casevars==onevar)[1],"Gene"], contlist[which(contvars==onevar)[1],"Gene"]) 
                c(tmpg,onevar,a,b,n.case,n.cont,(a/n.case)/(b/n.cont),binom.test(a,a+b,n.case/(n.case+n.cont))$p.value)
            })
            oneTable <- t(oneTable)
        }
    }else{
        if(flag == 1){
            a <- length(unique(caselist[,"Subject_ID"]))
            b <- length(unique(contlist[,"Subject_ID"]))
            b <- ifelse(b * coe > n.cont, n.cont, floor(b * coe))
            oneTable <- matrix(c(0,0,a,b,n.case-a,n.cont-b, (a/(n.case-a))/(b/(n.cont-b)), fisher.test(matrix(c(a,b,n.case-a,n.cont-b),2,2))$p.value),nrow=1,ncol=8)
        }else if(flag == 2){
            genes <- union(caselist[,"Gene"],contlist[,"Gene"])
            oneTable <- sapply(genes,function(gene){
                a <- length(unique(caselist[caselist[,"Gene"] %in% gene,"Subject_ID"]))
                b <- length(unique(contlist[contlist[,"Gene"] %in% gene,"Subject_ID"]))
                b <- ifelse(b * coe > n.cont, n.cont, floor(b * coe))
                c(gene,length(union(casevars[caselist[,"Gene"] %in% gene],contvars[contlist[,"Gene"] %in% gene])),a,b,n.case-a,n.cont-b, (a/(n.case-a))/(b/(n.cont-b)), fisher.test(matrix(c(a,b,n.case-a,n.cont-b),2,2))$p.value)
            })
            oneTable <- t(oneTable)
        }else if(flag == 3){
            vars <- unique(casevars)
            oneTable <- sapply(vars,function(onevar){
                a <- length(unique(caselist[casevars %in% onevar,"Subject_ID"]))
                b <- length(unique(contlist[contvars %in% onevar,"Subject_ID"]))
                tmpg <- ifelse(onevar %in% casevars, caselist[which(casevars==onevar)[1],"Gene"], contlist[which(contvars==onevar)[1],"Gene"]) 
                c(tmpg,onevar,a,b,n.case-a,n.cont-b,(a/(n.case-a))/(b/(n.cont-b)), fisher.test(matrix(c(a,b,n.case-a,n.cont-b),2,2))$p.value)
            })
            oneTable <- t(oneTable)
        }
    }

    cols <- c("Gene","Variant","#in_case","#in_cont","n.case","n.cont","Folds","Pvalue")
    colnames(oneTable) <- cols
    
    ##======================print information to logFile======================================
    print_log(paste("burden_test: the number of index cases used is ",n.case,sep=""))
    print_log(paste("burden_test: the number of corresponding controls used is ",n.cont,sep=""))
    print_log(paste("burden_test: there are totally ",dim(oneTable)[1]," tests",sep=""))
    print_log(paste("burden_test function is done!", date(),sep=" "))
    ##========================================================================================
    oneTable
}

qqplot_variantLevel <- function(burdentable,outputpath,Panelg,varT,varstr,caselist){
    
    ## print to log
    print_log(paste("qqplot_variantLevel function is running ...", date(),sep=" "))
   
    if(!file.exists(outputpath)){ dir.create(outputpath, showWarnings = TRUE, recursive = FALSE);}
    casevars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
    pT <- read.delim(burdentable)
    for(i in 1:length(varT)){
        subs <- which((pT[,"Variant"] %in% casevars[caselist[,"VariantClass"] %in% varT[[i]]]) & (pT[,"Folds"] >= 1))
        pval <- pT[subs,"Pvalue"]
        subcol <- which(pT[subs,"Gene"]%in% Panelg)
        pdf(file=paste(outputpath,varstr[i],".pdf",sep=""),width=12,height=10)
        par(mai=c(2,1,1,1))
        PQQ(pval,subcol=subcol,main=paste("QQ plots of ", varstr[i],sep=""),labels=pT[subs,"Gene"])
        dev.off()
    }
    
    ## 
    print_log(paste("qqplot_variantLevel is done!", date(),sep=" "))
    
}

PQQ <- function(pval,ps=0.05,subcol=NULL,main="",labels=NULL,nlab=20){
    ##nlab: number of points with labels
    
    x <- sort(-log10(runif(length(pval),0,1)))
    y <- sort(-log10(pval))
    yix <- sort(-log10(pval), index.return = TRUE)$ix
    
    plot(x,y,xlab="Excepted P-value (-log10 scale)",ylab="Observed P-value (-log10 scale)",main=main,xaxs="i",yaxs="i",xlim=c(0,max(c(x,y))+0.1),ylim=c(0,max(c(x,y))+0.1),cex.lab=1.5,cex.axis=1.5)
    lines(c(0,-log10(ps)),c(-log10(ps),-log10(ps)),type="l",col=1,lty=2)
    lines(c(-log10(ps),-log10(ps)),c(0,-log10(ps)),type="l",col=1,lty=2)
    abline(0,1,lwd=3,col=1)
    if(!is.null(subcol)){
        subs <- match(subcol,yix)
        lines(x[subs],y[subs],col=2,type="p",pch=19)
    }
    if(!is.null(labels)){
        library(wordcloud)
        subs <- which((y > -log10(ps)) & (x > -log10(ps)))
        if(length(y)>=nlab){ tmpsubs <- (length(y)-nlab+1):length(y);
        }else{ tmpsubs <- 1:length(y);}
        subs <- intersect(subs,tmpsubs)
        
        if(length(subs) > 1){
            nc <- wordlayout(x[subs],y[subs],labels[yix[subs]],cex=0.8)
            text(nc[,1],nc[,2],labels[yix[subs]],cex=0.8)
            #textplot(x[subs],y[subs],words=labels[yix[subs]],col=2,cex=0.8)
            if(!is.null(subcol)){
                insub <- match(subcol,yix)
                insub <- intersect(insub,subs)
                if( length(insub)>0 ) text(nc[match(insub,subs),1],nc[match(insub,subs),2],labels[yix[insub]],cex=0.8,col=2)
            }
        }else if(length(subs)==1){
		if(!is.null(subcol)){
			text(x[subs],y[subs],labels[yix[subs]],cex=0.8,col=1+(subs %in% match(subcol,yix)))
		}else{
			text(x[subs],y[subs],labels[yix[subs]],cex=0.8,col=1)
		}
	}
    }
    
}

## compute the corrected coefficient of HISP case and control
coeHisp <- function(caselist,contlist){

    ## based on the number of rare synonymou variant
    n.case <- length(unique(caselist[,"Subject_ID"]))
    n.cont <- length(unique(contlist[,"Subject_ID"]))
    syn <- "synonymousSNV"
    
    syn.case <- sum(caselist[,"VariantClass"] %in% syn)
    syn.cont <- sum(contlist[,"VariantClass"] %in% syn)

    coe <- (syn.case/n.case)/(syn.cont/n.cont)
    
    coe

}


plot_hq <- function(x,y,type,col,main,xlab,ylab,ylim,xlim){
    plot(x,y,type=type,col=col,main=main,xlab=xlab,ylab=ylab,ylim=ylim,xlim=xlim,cex.lab=1.5,cex.axis=1.5)
}

qwt <- function(x,filer,sep="\t",flag=0){
    if(flag==0){
        write.table(x,file=filer,quote=FALSE,row.names=FALSE,col.names=FALSE,sep=sep)
    }
    ## row.names = TRUE
    if(flag==1){
        write.table(x,file=filer,quote=FALSE,row.names=TRUE,col.names=FALSE,sep=sep)
    }
    ## col.names = TRUE
    if(flag==2){
        write.table(x,file=filer,quote=FALSE,row.names=FALSE,col.names=TRUE,sep=sep)
    }
    ## both
    if(flag==3){
        write.table(x,file=filer,quote=FALSE,row.names=TRUE,col.names=TRUE,sep=sep)
    }
    ## ALL 
    if(flag==4){
        write.table(x,file=filer,quote=TRUE,row.names=TRUE,col.names=TRUE,sep=sep)
    }
}

print_log <- function(printstr){
    
    ## write log file for each run
    logFile = "misc_run.log"
    cat(printstr, file=logFile, append=TRUE, sep = "\n")

}
