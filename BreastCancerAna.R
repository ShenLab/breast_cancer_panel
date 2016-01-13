## ===============================================================================
## burden test
## ===============================================================================
source("misc.R")
source("sourcefiles.R")

BreastCancerAna <- function(sig=TRUE,Ecut=0.001,hotf=1,swi=1){
### sig: test for singleton variant only or not
### Ecut: ExAC frequency cutoff
### hotf: 1: hotspots predicted based on HMM (hongjian); 2: hotspots defined by cosmic site with recurrent >=3; 3: not filter missense by hotspot; 4: not filter missense byhotspot but with the hotspot labels (by HMM hongjian)
### swi: 1: Jewish case and control; 2: Hispanic case and contorl
    
## =========================input files===========================================
if(hotf==1){ hotspotfile <- "../data/hotspots/HMM_hotspots_11_12.txt"; ## HongjianPred hotspots file
}else if(hotf==2){ hotspotfile <- "../data/hotspots/cosmic_hotspots_3.txt"; ## COSMIC hotspots file
}else if(hotf==3){ hotspotfile="";
}else if(hotf==4){ hotspotfile <- "../data/hotspots/HMM_hotspots_11_12.txt";}
## switch to the right files for burden test
if(swi==1){
        vburdenfile <- "../resultf/variant_level_burden_anno_Fre_Pseducont.txt"  ### single variant test files with all annotations; look src.R
	    alleleFrefile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/AJ_CONTROLS/combined_variant_call/NonBC_Frequencies.expanded.tsv" ##AJ: not in our cohort 
    	caselistf <- "../data/Rdata/AJcaselist_12_17"
    	contlistf <- "../data/Rdata/AJcontlist_12_17"
	    casestaf <- "../data/contAJ_20151209.hardfiltered.stats_hq.tsv" ##using the variant from AJ case and control joint calling result
	    contstaf <- "../data/contAJ_20151209.hardfiltered.stats_hq.tsv"
}else if(swi==2){
	    vburdenfile = "../resultf/HISP_variant_level_burden_anno_Fre_Pseducont.txt"	
    	alleleFrefile <- NULL
    	caselistf <- "../data/Rdata/HIcaselist_1_5"
    	contlistf <- "../data/Rdata/HIcontlist_1_5"
	    casestaf <- "../data/BR.origin.stats_12_28.tsv"
    	contstaf <- "../data/Hispanic.stats_12_28.tsv"
}

##====================define the output files===================================================================
pop=c("Jewish","Hispanic")
outputpath <- paste("../resultf/burdentest_",sig,"_singleton_",Ecut,"_ExACFre_",gsub(".txt","",basename(hotspotfile)),"_",hotf,"_",pop[swi],"_",format(Sys.Date(),format="%B_%d_%Y"),"/",sep="")
genesetsVarfile <- paste(outputpath,"Panel_genes_variantlist.txt",sep="")
restVarfile <- paste(outputpath,"Variantlist.txt",sep="")
NumOfVarfile <- paste(outputpath,"NumOfVariant.txt",sep="")
if(!file.exists(outputpath)){ dir.create(outputpath, showWarnings = TRUE, recursive = FALSE);}

## =========================variant class definition and test gene sets=================================
print_log(paste("burden_test parameters are singleton only: ",sig,"; ExAC frequency: ",Ecut,"; hotspots file: ",hotspotfile,sep=" "))
VariantClass <- c(".","frameshiftdeletion","frameshiftinsertion","none","nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV","stopgain","stoploss","synonymousSNV","unknown")
lof <- c("frameshiftdeletion","frameshiftinsertion")
mis <- "nonsynonymousSNV"
indel <- c("nonframeshiftdeletion","nonframeshiftinsertion")
syn <- "synonymousSNV"
unknown <- "unknown"
stopins <- c("stopgain","stoploss")
splices <- c("none",".")
singleLOF <- c("stopgain","stoploss","none",".")
### only for QQ plots
varT <- list(lof,mis,indel,syn,stopins,splices,singleLOF)
varTnames <- c("lof_indel","d-mis","mis-indel","Synonymous","Stopgain_loss","splicing","singleLOF")

TSg <- unlist(read.table(TSfile))
DRg <- unlist(read.table(cancerdriverfile)) ## cancer drivers
DNAreg <- unlist(read.table(DNArepairfile))
Panelg <- unlist(read.table(Panelgfile))
Panelg <- setdiff(Panelg,c(TSg,DRg,DNAreg))
genes <- unique(c(TSg,DRg,DNAreg,Panelg))
Gtop <- read.table(GTExfile)
allgenes <- union(Gtop[,1],c(TSg,DRg,DNAreg,Panelg))


## =========================variant list filtering===========================================
##getVariantlist(HIcasevariantpath,HIBRfile,namestr=".AllVariants.tsv","../data/Rdata/HIcaselist_1_5")
##getVariantlist(HIcontvariantpath,HIcontrolfile,namestr=".tsv","../data/Rdata/HIcontlist_1_5")
### AJ lists
##getVariantlist(AJcasevariantpath,AJBRfile,namestr=".tsv","../data/Rdata/AJcaselistBR_12_17")
##getVariantlist(AJcontvariantpath,AJBRfile,namestr=".tsv","../data/Rdata/AJcaselist_12_17")
##getVariantlist(AJcontvariantpath,AJcontrolfile,namestr=".tsv","../data/Rdata/AJcontlist_12_17")
### ====================================getVariantlist save as Rdata===========================
print_log("variant list filtering ...")
load(caselistf)
caselist <- onelist
rm(onelist)
indexcases <- getindexcase(phenofile) ## only index cases
caselist <- caselist[caselist[,"Subject_ID"] %in% indexcases, ]
exSamples <- excluded_samples() ## exclude subjects
caselist <- caselist[!(caselist[,"Subject_ID"] %in% exSamples), ]
## =========================get control variant list
load(contlistf)
contlist <- onelist
rm(onelist)

## variant count statistics: double check for batches effect
VariantSta(casestaf,contstaf,unique(caselist[,"Subject_ID"]),unique(contlist[,"Subject_ID"]),paste(outputpath,"variantSta/",sep=""))

##============================== keep singleton variants only
if(sig){
    print_log(paste("variant_filtering parameters: singleton variant only ",sig,sep=" "))
    casevars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
    contvars <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
    singletons <- singleton_variants(casevars,contvars)
    caselist <- caselist[casevars %in% singletons,]
    contlist <- contlist[contvars %in% singletons,]
}
## variant filtering: filters <- c("filtered","ExACfreq","VCFPASS","noneSegmentalDup","meta-SVM_PP2","hotspot","alleleFre")
caselist <- variant_filtering(caselist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf=hotspotfile,alleleFrefile,popcut=0.05)
print(hotf)
if(hotf==4){ caselist <- caselist[caselist[,"ExACfreq"] & caselist[,"VCFPASS"] & caselist[,"noneSegmentalDup"] & caselist[,"meta-SVM_PP2"] & caselist[,"alleleFre"], ]; }else{ caselist <- caselist[caselist[,"filtered"], ];}
contlist <- variant_filtering(contlist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf=hotspotfile,alleleFrefile,popcut=0.05)
if(hotf==4){ contlist <- contlist[contlist[,"ExACfreq"] & contlist[,"VCFPASS"] & contlist[,"noneSegmentalDup"] & contlist[,"meta-SVM_PP2"] & contlist[,"alleleFre"], ]; }else{ contlist <- contlist[contlist[,"filtered"], ];}


## =========================burden test for gene, variant=========================================
print_log("burden test for gene, variant ...")
## burden test for all genes, tumor suppressors, cancer drivers genes
Ggs <- c("top 25%","top 50%","top 75%","top 100%")
GgL <- list(Gtop[1:floor(0.25*dim(Gtop)[1]),1], Gtop[1:floor(0.5*dim(Gtop)[1]),1], Gtop[1:floor(0.75*dim(Gtop)[1]),1],allgenes)
genesets <- list(TSg,DRg,DNAreg,Panelg)
genesetnames <- c("Tumor suppressors","Cancer drivers","DNA repairs","Panel genes")
vartypes <- list(stopins,splices,singleLOF,lof,mis,indel,NULL,syn,unknown)
vartypenames <- c("stopgain_loss","splicing","singleLOF","indelLOF","D-MIS","Indels","ALL variants","Synonymous","Unknown")

## eliminate batch effect for HISP based on rare synonymous
if(swi==1){coe=1;}else if(swi==2){ coe = coeHisp(caselist,contlist);}
print("Corrected coe: ");
print(coe);

setburdens <- c()
for(i in 1:length(genesets)){
    for(j in 1:length(vartypes)){
        for(k in 1:length(Ggs)){
            tmp <- burden_test(caselist,contlist,testset=intersect(genesets[[i]],GgL[[k]]),testtype=vartypes[[j]],flag=1,sig,coe)
            tmp[1,1:2] <- c(genesetnames[i],vartypenames[j])
            tmp <- cbind(tmp,Ggs[k])
            setburdens <- rbind(setburdens,tmp)
        }
    }
}
## single gene level burden test
tablelist <- list()
for(i in 1:length(vartypes)){
    oneTable <- burden_test(caselist,contlist,testset=NULL,testtype=vartypes[[i]],flag=2,sig=sig,coe=coe)
    oneTable <- oneTable[order(as.numeric(oneTable[,"Pvalue"])), ]
    tablelist[[i]] <- oneTable
}

## =========================output burden test to files=========================================
print_log("Output burden files ...")
### output file names: gene sets and single gene burden files
outputpath1 <- paste(outputpath,"burden/",sep="")
if(!file.exists(outputpath1)){ dir.create(outputpath1, showWarnings = TRUE, recursive = FALSE);}
setburdenfile <- paste(outputpath1,"gene_variant_set.burden.txt",sep="") ### gene set and variant types burden test file
qwt(setburdens,setburdenfile,flag=2)
for(i in 1:length(vartypes)){
    onefile <- paste(outputpath1,vartypenames[i],"_level.burden.txt",sep="")
    qwt(tablelist[[i]],onefile,flag=2)
}

##======output single variant tables == give variant list information==================
vartys <- c(stopins,splices,singleLOF,mis,lof,indel) ## slected variants
cols0 <- colnames(caselist)

Vlist <- caselist[caselist[,"VariantClass"] %in% vartys, ]
vars <- paste(Vlist[,1],Vlist[,2],Vlist[,4],Vlist[,5],sep="_")
univar <- unique(vars)
Vlist <- Vlist[match(univar,vars),]

Vlist[,"Gene_sets"] <- "non-Panel_genes"
Vlist[Vlist[,"Gene"] %in% DNAreg,"Gene_sets"] <- "DNA_repair"
Vlist[Vlist[,"Gene"] %in% Panelg,"Gene_sets"] <- "Panel_genes"
Vlist[Vlist[,"Gene"] %in% TSg,"Gene_sets"] <- "Tumor_suppressor"
Vlist[Vlist[,"Gene"] %in% DRg,"Gene_sets"] <- "Cancer_driver"
Vlist[Vlist[,"Gene"] %in% intersect(TSg,DRg),"Gene_sets"] <- "Cancer_driver_Tumor_suppressor"
Vlist[,"variant_type"] <- ""
Vlist[Vlist[,"VariantClass"] %in% lof,"variant_type"] <- "indels_inframe" 
Vlist[Vlist[,"VariantClass"] %in% c(stopins,splices,singleLOF),"variant_type"] <- "LOF"
Vlist[Vlist[,"VariantClass"] %in% indel,"variant_type"] <- "indels_nonframe"
Vlist[Vlist[,"VariantClass"] %in% mis,"variant_type"] <- "d-mis"

varTable <- read.delim(vburdenfile,check.names=FALSE)
Vlist <- cbind(Vlist,varTable[match(univar,varTable[,"Variant"]) , setdiff(colnames(varTable),colnames(Vlist)) ])

# corrected Pvalue by each variant type
Vlist[,"CorrectedP"] <- 1
NumVar <- c(length(univar[Vlist[,"VariantClass"] %in% lof]), length(univar[Vlist[,"VariantClass"] %in% c(stopins,splices,singleLOF)]), length(univar[Vlist[,"VariantClass"] %in% indel]), length(univar[Vlist[,"VariantClass"] %in% mis]))
Vlist[Vlist[,"variant_type"] == "indels_inframe","CorrectedP"] <-  as.numeric(Vlist[Vlist[,"variant_type"] == "indels_inframe","Pvalue"]) * NumVar[1]
Vlist[Vlist[,"variant_type"] == "LOF","CorrectedP"] <-  as.numeric(Vlist[Vlist[,"variant_type"] == "LOF","Pvalue"]) * NumVar[2]
Vlist[Vlist[,"variant_type"] == "indels_nonframe","CorrectedP"] <-  as.numeric(Vlist[Vlist[,"variant_type"] == "indels_nonframe","Pvalue"]) * NumVar[3]
Vlist[Vlist[,"variant_type"] == "d-mis","CorrectedP"] <-  as.numeric(Vlist[Vlist[,"variant_type"] == "d-mis","Pvalue"]) * NumVar[4]
Vlist[as.numeric(Vlist[,"CorrectedP"]) > 1,"CorrectedP"] <- 1
NumVar <- paste(c("indels_inframe","LOF","indels_nonframe","d-mis")," : ",NumVar,sep="")

#### selected columns
frecols <- c("N.index.AJ","N.pseudoCont.AJ","N.case.AJ","N.non_case.AJ","N.cont.AJ","N.index.HI","N.pseudoCont.HI","N.case.HI","N.non_case.HI","N.cont.HI")
samcols <- c("index.AJ","pseudoCont.AJ","case.AJ","non_case.AJ","cont.AJ","index.HI","pseudoCont.HI","case.HI","non_case.HI","cont.HI")
oddcols <- c("Odds_AJ_pseudo","p_AJ_pseudo","Odds_AJ_cont","p_AJ_cont","Odds_HI_pseudo","p_HI_pseudo","Odds_HI_cont","p_HI_cont")
AJcols <- c("N.index.AJ","N.pseudoCont.AJ","N.case.AJ","N.non_case.AJ","N.cont.AJ","Odds_AJ_pseudo","p_AJ_pseudo","Odds_AJ_cont","p_AJ_cont","index.AJ","pseudoCont.AJ","case.AJ","non_case.AJ","cont.AJ")
HIcols <- c("N.index.HI","N.pseudoCont.HI","N.case.HI","N.non_case.HI","N.cont.HI","Odds_HI_pseudo","p_HI_pseudo","Odds_HI_cont","p_HI_cont","index.HI","pseudoCont.HI","case.HI","non_case.HI","cont.HI")
## ordered columns
if(sig){ 
    cols <- c("Gene_sets","variant_type","Gene","Variant", setdiff(cols0,"Gene"), "hotspot", samcols, frecols)
    Vlist <- Vlist[,cols] 
    Vlist <- Vlist[order(Vlist[,"Gene_sets"]),]
}else{  
    if(swi==1){ cols <- c("Gene_sets","variant_type","Gene","Variant","#in_case","#in_cont","Folds","Pvalue","CorrectedP",setdiff(AJcols,c("N.index.AJ","N.cont.AJ","Odds_AJ_cont","p_AJ_cont")),HIcols,setdiff(cols0,"Gene"));}
    if(swi==2){ cols <- c("Gene_sets","variant_type","Gene","Variant","#in_case","#in_cont","Folds","Pvalue","CorrectedP",setdiff(HIcols,c("N.index.HI","N.cont.HI","Odds_HI_cont","p_HI_cont")),AJcols,setdiff(cols0,"Gene"));}
    Vlist <- Vlist[,cols]
    ### selected columns and re-orders
    Vlist <- Vlist[order(as.numeric(Vlist[,"Pvalue"])),]
}

##=======================write result======================
qwt(Vlist[Vlist[,"Gene"] %in% genes,],file=genesetsVarfile,flag=2)
qwt(Vlist[Vlist[,"Gene_sets"] %in% "non-Panel_genes",],file=restVarfile,flag=2)
qwt(NumVar,file=NumOfVarfile)
## qq plot
if(!sig){ qqplot_variantLevel(vburdenfile,paste(outputpath,"qqplots/",sep=""),genes,varT,varTnames,caselist);}

##========end of test=============
}
