## ===============================================================================
## Jewish burden test
## ===============================================================================
source("misc.R")

BreastCancerAna <- function(sig=TRUE,Ecut=0.001,hotf=1){

## =========================input files===========================================
if(hotf==1){ hotspotfile <- "../data/hotspots/HMM_hotspots_11_12.txt"; ## HongjianPred hotspots file
}else{ hotspotfile <- "../data/hotspots/cosmic_hotspots_3.txt"; }## COSMIC hotspots file
outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"  ## outlier samples
BRCA1_2pathogenicfile <- "../data/phenotype/BRCA1_2.txt" ## BRCA1/2 pathogenic samples file
alleleFrefile <- "/home/local/ARCS/ads2202/data/AJ_PCA/combined_variant_call/NonBC_Frequencies.expanded.tsv" ## non-breast cancer samples allle frequency
TSfile <- "../data/hotspots/Tumor_suppressors_11_11.txt" ## collected tumor suppressors
cancerdriverfile <- "../data/hotspots/Cancer_driver_11_6.txt" ## cancer drivers
DNArepairfile <- "../data/hotspots/DNA_repair_11_6.txt" ## DNA repair genes
Panelgfile <-  "../../genelist/Genelist2.txt" ## Panel genes
GTExfile <- "/home/local/ARCS/qh2159/breast_cancer/geneexpression/GTEx/expg_ranked.txt" ## GTEx expressed ranked gene
phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv" ## phenotype file to get index cases only
dupIDs <- c("222357","222966") ## duplicated with Subject ID 222966 ## duplicated subject
## variant files and case, control samples ## check the log file to see the variant filtering details
casevariantpath <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"
AJcasefile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/tablesandlists/All_AJ_samples.list"
contvariantpath <- "/home/local/ARCS/qh2159/breast_cancer/variants/AJconVariantCalling/"
AJcontrolfile <- "/home/local/ARCS/qh2159/breast_cancer/variants/data/AJs_585.txt"
Cohortfile <- "../resultf/BreastCancer_VariantList_11_12"

print_log(paste("burden_test parameters are singleton only: ",sig,"; ExAC frequency: ",Ecut,"; hotspots file: ",hotspotfile,sep=" "))
## =========================variant class definition and test gene sets=================================
VariantClass <- c(".","frameshiftdeletion","frameshiftinsertion","none","nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV","stopgain","stoploss","synonymousSNV","unknown")
lof <- c("frameshiftdeletion","frameshiftinsertion","none","stopgain","stoploss",".")
mis <- "nonsynonymousSNV"
indel <- c("nonframeshiftdeletion","nonframeshiftinsertion")
syn <- "synonymousSNV"
unknown <- "unknown"

TSg <- unlist(read.table(TSfile))
DRg <- unlist(read.table(cancerdriverfile)) ## cancer drivers
DNAreg <- unlist(read.table(DNArepairfile))
Panelg <- unlist(read.table(Panelgfile))
Panelg <- setdiff(Panelg,c(TSg,DRg,DNAreg))
Gtop <- read.table(GTExfile)
allgenes <- union(Gtop[,1],c(TSg,DRg,DNAreg,Panelg))

print_log("variant list filtering ...")
## =========================variant list filtering===========================================
##getVariantlist(casevariantpath,AJcasefile,namestr=".AllVariants.tsv","../data/Rdata/AJcaselist_11_9")
##getVariantlist(contvariantpath,AJcontrolfile,namestr=".tsv","../data/Rdata/AJcontlist_11_9")
### ====================================getVariantlist save as Rdata==
load("../data/Rdata/AJcaselist_11_9")
caselist <- onelist
rm(onelist)
## only index cases
indexcases <- getindexcase(phenofile)
caselist <- caselist[caselist[,"Subject_ID"] %in% indexcases, ]
## exclude outlier subjects
outliers <- unlist(read.table(outlierfile))
outliers <- sapply(1:length(outliers),function(i) {tmp=unlist(strsplit(outliers[i],"_"));setdiff(gsub("\\D","",tmp),c("",paste("00",1:9,sep="")));})
caselist <- caselist[!(caselist[,"Subject_ID"] %in% outliers), ]
## exclude the duplicated subjects
if(all(dupIDs %in% caselist[,"Subject_ID"])) caselist <- caselist[caselist[,"Subject_ID"]!=dupIDs[1], ]
## ====================================get control variant list==
load("../data/Rdata/AJcontlist_11_9")
contlist <- onelist
rm(onelist)
## remove out the BRCA1/2 pathogenic cases
tmp <- read.delim(BRCA1_2pathogenicfile)
pathogenic_sample <- tmp[tmp[,1] %in% c("likely pathogenic","pathogenic"),"Subject_ID"]
caselist <- caselist[!(caselist[,"Subject_ID"] %in% pathogenic_sample), ]
contlist <- contlist[!(contlist[,"Subject_ID"] %in% pathogenic_sample), ]
## keep singleton variants only
if(sig){
    print_log(paste("variant_filtering parameters: singleton variant only ",sig,sep=" "))
    casevars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
    contvars <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
    singletons <- singleton_variants(casevars,contvars)
    caselist <- caselist[casevars %in% singletons,]
    contlist <- contlist[contvars %in% singletons,]
}
## variant filtering: filters <- c("filtered","ExACfreq","VCFPASS","noneSegmentalDup","meta-SVM_PP2","singleton","hotspot")
caselist <- variant_filtering(caselist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf=hotspotfile,alleleFrefile,popcut=0.05)
caselist <- caselist[caselist[,"filtered"], ]
contlist <- variant_filtering(contlist,mis,Ecut=Ecut,segd=0.95,pp2=TRUE,hotf=hotspotfile,alleleFrefile,popcut=0.05)
contlist <- contlist[contlist[,"filtered"], ]


print_log("burden test for gene, variant ...")
## =========================burden test for gene, variant=========================================
## burden test for all genes, tumor suppressors, cancer drivers genes
Ggs <- c("top 25%","top 50%","top 75%","top 100%")
GgL <- list(Gtop[1:floor(0.25*dim(Gtop)[1]),1], Gtop[1:floor(0.5*dim(Gtop)[1]),1], Gtop[1:floor(0.75*dim(Gtop)[1]),1],allgenes)
genesets <- list(TSg,DRg,DNAreg,Panelg)
genesetnames <- c("Tumor suppressors","Cancer drivers","DNA repairs","Panel genes")
vartypes <- list(lof,mis,indel,NULL,syn,unknown)
vartypenames <- c("LOF","D-MIS","Indels","ALL variants","Synonymous","Unknown")

setburdens <- c()
for(i in 1:length(genesets)){
    for(j in 1:length(vartypes)){
        for(k in 1:length(Ggs)){
            tmp <- burden_test(caselist,contlist,testset=intersect(genesets[[i]],GgL[[k]]),testtype=vartypes[[j]],flag=1,sig)
            tmp[1,1:2] <- c(genesetnames[i],vartypenames[j])
            tmp <- cbind(tmp,Ggs[k])
            setburdens <- rbind(setburdens,tmp)
        }
    }
}
## single gene level burden test
LOFTable <- burden_test(caselist,contlist,testtype=lof,flag=2,sig=sig)
LOFTable <- LOFTable[order(-as.numeric(LOFTable[,"Folds"]),as.numeric(LOFTable[,"Pvalue"])), ]
MISTable <- burden_test(caselist,contlist,testtype=mis,flag=2,sig=sig)
MISTable <- MISTable[order(-as.numeric(MISTable[,"Folds"]),as.numeric(MISTable[,"Pvalue"])), ]
indelTable <- burden_test(caselist,contlist,testtype=indel,flag=2,sig=sig)
indelTable <- indelTable[order(-as.numeric(indelTable[,"Folds"]),as.numeric(indelTable[,"Pvalue"])), ]
synTable <- burden_test(caselist,contlist,testtype=syn,flag=2,sig=sig)
synTable <- synTable[order(-as.numeric(synTable[,"Folds"]),as.numeric(synTable[,"Pvalue"])), ]
unkTable <- burden_test(caselist,contlist,testtype=unknown,flag=2,sig=sig)
unkTable <- unkTable[order(-as.numeric(unkTable[,"Folds"]),as.numeric(unkTable[,"Pvalue"])), ]
## single variant level burden test
variantTable <- burden_test(caselist,contlist,flag=3,sig=sig)
variantTable <- variantTable[order(-as.numeric(variantTable[,"Folds"]),as.numeric(variantTable[,"Pvalue"])), ]

print_log("Output files ...")
## =========================output burden test to files=========================================
outputpath <- paste("../resultf/burdentest_",sig,"_",Ecut,"_",gsub(".txt","",basename(hotspotfile)),"/",sep="")
if(!file.exists(outputpath)){ dir.create(outputpath, showWarnings = TRUE, recursive = FALSE);}
### output file names
setburdenfile <- paste(outputpath,"gene_variant_set.burden.txt",sep="") ### gene set and variant types burden test file
loffile <- paste(outputpath,"LOF_level.burden.txt",sep="") ###LOF: single gene level burden test
misfile <- paste(outputpath,"MIS_level.burden.txt",sep="") ###MIS: single gene level burden test                 
indelfile <- paste(outputpath,"indel_level.burden.txt",sep="") ###indel: single gene level burden test
synfile <- paste(outputpath,"syn_level.burden.txt",sep="") ###synonymous
unkfile <- paste(outputpath,"Unknown_level.burden.txt",sep="") ### unknown
vburdenfile <- paste(outputpath,"variant_level.burden.txt",sep="") ### single variant level burden test
## write to files
qwt(setburdens,setburdenfile,flag=2)
qwt(LOFTable,loffile,flag=2)
qwt(MISTable,misfile,flag=2)
qwt(indelTable,indelfile,flag=2)
qwt(synTable,synfile,flag=2)
qwt(unkTable,unkfile,flag=2)
qwt(variantTable,vburdenfile,flag=2)


## =========================output indels for IGV and variant tables=========================================
ListgoSimfile <- paste(outputpath,"VariantList_Simple.txt",sep="")
indelsfile <- paste(outputpath,"indels_",Ecut,".txt",sep="")
genesetsVarfile <- paste(outputpath,"geneset_variantlist.txt",sep="")
##===============================step 1: give variant list information==================
varTable <- read.delim(vburdenfile)
genes <- unique(c(TSg,DRg,DNAreg,Panelg))
vartys <- c(mis,lof,indel)
variantlist <- caselist[caselist[,"VariantClass"] %in% vartys & caselist[,"Gene"] %in% genes, ]
vars <- paste(variantlist[,1],variantlist[,2],variantlist[,4],variantlist[,5],sep="_")
n.var <- length(vars)
rowTitle <- matrix(0,n.var,10)
rowTitle[variantlist[,"Gene"] %in% DNAreg,1] <- "DNA_repair"
rowTitle[variantlist[,"Gene"] %in% Panelg,1] <- "Panel_genes"
rowTitle[variantlist[,"Gene"] %in% TSg,1] <- "Tumor_suppressor"
rowTitle[variantlist[,"Gene"] %in% DRg,1] <- "Cancer_driver"
rowTitle[variantlist[,"Gene"] %in% intersect(TSg,DRg),1] <- "Cancer_driver_Tumor_suppressor"
rowTitle[variantlist[,"VariantClass"] %in% lof,2] <- "LOF" 
rowTitle[variantlist[,"VariantClass"] %in% indel,2] <- "indels"
rowTitle[variantlist[,"VariantClass"] %in% mis,2] <- "d-mis"
rowTitle[,3] <- variantlist[,"Gene"]
rowTitle[,4:10] <- as.matrix(varTable[match(vars,varTable[,2]),c(2:8)])
variantlist <- variantlist[,colnames(variantlist)!="Gene"]
variantlist <- cbind(rowTitle,variantlist)
variantlist <- variantlist[order(variantlist[,1],variantlist[,2],-as.numeric(variantlist[,9]),as.numeric(variantlist[,10])),]
if(sig){
    colnames(variantlist)[1:10] <- c("Gene_sets","variant_type","Gene","Variant","#cases","#controls","n.case","n.cont","Folds","pvalue")
}else{
    colnames(variantlist)[1:10] <- c("Gene_sets","variant_type","Gene","Variant","#cases carrier","#controls carrier","#case non-carrier","#control non-carrier","OddsRatio","pvalue")
}
indelVars <- variantlist[variantlist[,"VariantClass"] %in% c("frameshiftdeletion","frameshiftinsertion","nonframeshiftdeletion","nonframeshiftinsertion"),c("Chromosome","Position","Subject_ID")]
## ===========step 2: give phenotype information on family and other cohort cases=============
load(Cohortfile)
pheno <- read.csv(phenofile)
varsCheck <- variantlist[,"Variant"]
contvars <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
casevars <- paste(caselist[,1],caselist[,2],caselist[,4],caselist[,5],sep="_")
cohortvars <- paste(onelist[,1],onelist[,2],onelist[,4],onelist[,5],sep="_")
varvariantlist <- c()
for(i in 1:length(varsCheck)){
    subj <- unique(caselist[casevars %in% varsCheck[i],"Subject_ID"]) ## only care about case variants
    famj <- pheno[pheno[,3] %in% subj,1]
    subs1 <- (onelist[,"Subject_ID"] %in% pheno[pheno[,1] %in% famj,3])
    subs2 <- (cohortvars %in% varsCheck[i])
    onefamid <- onelist[subs1 & subs2,"Subject_ID"]
    id1 <- onefamid[pheno[match(onefamid,pheno[,3]),"BreastCancer"]=="Yes"]
    id2 <- onefamid[pheno[match(onefamid,pheno[,3]),"BreastCancer"]=="No"]
    tmp <- onelist[!subs1 & subs2,"Subject_ID"]
    tmp <- intersect(tmp,pheno[,3])
    id3 <- tmp[pheno[match(tmp,pheno[,3]),"BreastCancer"]=="Yes"]
    id4 <- tmp[pheno[match(tmp,pheno[,3]),"BreastCancer"]=="No"]
    id5 <- contlist[contvars %in% varsCheck[i],"Subject_ID"]
    tmp <- c(paste(famj,sep="",collapse=":"), paste(id1,sep="",collapse=":"), paste(id2,sep="",collapse=":"),paste(id3,sep="",collapse=":"),paste(id4,sep="",collapse=":"),paste(id5,sep="",collapse=":"))
    varvariantlist <- rbind(varvariantlist,tmp)
}
colnames(varvariantlist) <- c("Family-ID","Yes-Family","No-Family","Yes-Cohort","No-Cohort","Control")
rownames(varvariantlist) <- rownames(variantlist)
variantlistCom <- cbind(varvariantlist,variantlist)
variantlistCom <- variantlistCom[order(variantlistCom[,"Gene_sets"],variantlistCom[,"variant_type"],variantlistCom[,"Gene"]),]
variantlistSim <- variantlistCom[,c("Gene","Variant","Subject_ID","Family-ID","Yes-Family","No-Family","Yes-Cohort","No-Cohort","Control","variant_type","#cases","#controls","AlleleFrequency.ExAC","GT","AD","Gene_sets")]
## step 3: write result
qwt(variantlistCom,file=genesetsVarfile,flag=2)
qwt(indelVars,file=indelsfile)
qwt(variantlistSim,file=ListgoSimfile,flag=2)


}
