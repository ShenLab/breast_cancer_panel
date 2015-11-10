## ===============================================================================
## Jewish burden test
## ===============================================================================
source("misc.R")

## =========================input files===========================================
## outlier samples
outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"
## BRCA1/2 pathogenic samples file
BRCA1_2pathogenicfile <- "../data/phenotype/BRCA1_2.txt"
## COSMIC hotspots file
hotspotfile <- "../data/hotspots/cosmic_hotspots_2.txt"
## non-breast cancer samples allle frequency
alleleFrefile <- "/home/local/ARCS/ads2202/data/AJ_PCA/combined_variant_call/NonBC_Frequencies.expanded.tsv"
## collected tumor suppressors and cancer drivers and DNA repair genes
TSfile <- "../data/hotspots/Tumor_suppressors_11_6.txt"
cancerdriverfile <- "../data/hotspots/Cancer_driver_11_6.txt" ## cancer drivers
DNArepairfile <- "../data/hotspots/DNA_repair_11_6.txt"
## variant files and case, control samples
casevariantpath <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"
AJcasefile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/tablesandlists/All_AJ_samples.list"
contvariantpath <- "/home/local/ARCS/qh2159/breast_cancer/variants/AJconVariantCalling/"
AJcontrolfile <- "/home/local/ARCS/qh2159/breast_cancer/variants/data/AJs_586.txt"
## check the log file to see the variant filtering details
## phenotype file to get index cases only
phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv"
## duplicated subject
dupIDs <- c("222357","222966") ## duplicated with Subject ID 222966

## =========================variant list filtering===========================================
##getVariantlist(casevariantpath,AJcasefile,namestr=".AllVariants.tsv","../data/Rdata/AJcaselist_11_9")
##getVariantlist(contvariantpath,AJcontrolfile,namestr=".tsv","../data/Rdata/AJcontlist_11_9")
### getVariantlist save as Rdata==================
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
## get control variant list======================
load("../data/Rdata/AJcontlist_11_9")
contlist <- onelist
rm(onelist)
## remove out the BRCA1/2 pathogenic cases
tmp <- read.delim(BRCA1_2pathogenicfile)
pathogenic_sample <- tmp[tmp[,1] %in% c("likely pathogenic","pathogenic"),"Subject_ID"]
caselist <- caselist[!(caselist[,"Subject_ID"] %in% pathogenic_sample), ]
contlist <- contlist[!(contlist[,"Subject_ID"] %in% pathogenic_sample), ]
##  synonynous variant
VariantClass <- c(".","frameshiftdeletion","frameshiftinsertion","none","nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV","stopgain","stoploss","synonymousSNV","unknown")
lof <- c("frameshiftdeletion","frameshiftinsertion","none","stopgain","stoploss")
mis <- c("nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV")
syn <- "synonymousSNV"
print_log(paste("Main: LOF variant defined as : ", paste(lof,sep="",collapse=","), sep=" "))
print_log(paste("Main: D-MIS variant defined as : ", paste(mis,sep="",collapse=","), sep=" "))
print_log(paste("Main: syn variant defined as : ", syn, sep=" "))
#caselist <- caselist[!(caselist[,"VariantClass"] %in% syn), ]
#contlist <- contlist[!(contlist[,"VariantClass"] %in% syn), ]
##  variant filtering
## filters <- c("filtered","ExACfreq","VCFPASS","noneSegmentalDup","meta-SVM_PP2","singleton","hotspot")
caselist <- variant_filtering(caselist,mis,Ecut=0.01,segd=0.95,pp2=TRUE,sig=FALSE,hotf=hotspotfile,alleleFrefile,popcut=0.05)
caselist <- caselist[caselist[,"filtered"], ]
contlist <- variant_filtering(contlist,mis,Ecut=0.01,segd=0.95,pp2=TRUE,sig=FALSE,hotf=hotspotfile,alleleFrefile,popcut=0.05)
contlist <- contlist[contlist[,"filtered"], ]


## =========================burden test for gene, variant=========================================
## burden test for all genes, tumor suppressors, cancer drivers genes
TSg <- unlist(read.table(TSfile))
DRg <- unlist(read.table(cancerdriverfile)) ## cancer drivers
DNAreg <- unlist(read.table(DNArepairfile))
genesets <- list(TSg,DRg,DNAreg,NULL)
genesetnames <- c("Tumor suppressors","Cancer drivers","DNA repairs","ALL genes")
vartypes <- list(lof,mis,mis,NULL,syn)
vartypenames <- c("LOF","D-MIS","Indels","ALL variants","Synonymous")
setburdens <- c()
for(i in 1:length(genesets)){
    for(j in 1:length(vartypes)){
        if(j==3){ indel <- TRUE;}else{indel=FALSE;}
        tmp <- burden_test(caselist,contlist,testset=genesets[[i]],testtype=vartypes[[j]],flag=1,indel)
        tmp[1,1:2] <- c(genesetnames[i],vartypenames[j])
        setburdens <- rbind(setburdens,tmp)
    }
}
## single gene level burden test
geneTable <- burden_test(caselist,contlist,flag=2)
geneTable <- geneTable[order(-as.numeric(geneTable[,"Folds"]),as.numeric(geneTable[,"Pvalue"])), ]
LOFTable <- burden_test(caselist,contlist,testtype=lof,flag=2)
LOFTable <- LOFTable[order(-as.numeric(LOFTable[,"Folds"]),as.numeric(LOFTable[,"Pvalue"])), ]
MISTable <- burden_test(caselist,contlist,testtype=mis,flag=2)
MISTable <- MISTable[order(-as.numeric(MISTable[,"Folds"]),as.numeric(MISTable[,"Pvalue"])), ]
indelTable <- burden_test(caselist,contlist,testtype=mis,flag=2,indel=TRUE)
indelTable <- indelTable[order(-as.numeric(indelTable[,"Folds"]),as.numeric(indelTable[,"Pvalue"])), ]
synTable <- burden_test(caselist,contlist,testtype=syn,flag=2)
synTable <- synTable[order(-as.numeric(synTable[,"Folds"]),as.numeric(synTable[,"Pvalue"])), ]
## single variant level burden test
variantTable <- burden_test(caselist,contlist,flag=3)
variantTable <- variantTable[order(-as.numeric(variantTable[,"Folds"]),as.numeric(variantTable[,"Pvalue"])), ]


## =========================output burden test to files=========================================
### output files
if(!file.exists("../resultf/burdentest/")){
    dir.create("../resultf/burdentest/", showWarnings = TRUE, recursive = FALSE)
}
### gene set and variant types burden test file
setburdenfile <- "../resultf/burdentest/gene_variant_set.burden.txt"
qwt(setburdens,setburdenfile,flag=2)
### single gene level burden test
gburdenfile <- "../resultf/burdentest/gene_level.burden.txt"
qwt(geneTable,gburdenfile,flag=2)
loffile <- "../resultf/burdentest/LOF_level.burden.txt"
qwt(LOFTable,loffile,flag=2)
misfile <- "../resultf/burdentest/MIS_level.burden.txt"
qwt(MISTable,misfile,flag=2)
indelfile <- "../resultf/burdentest/indel_level.burden.txt"
qwt(indelTable,indelfile,flag=2)
synfile <- "../resultf/burdentest/syn_level.burden.txt"
qwt(synTable,synfile,flag=2)
### single variant level burden test
vburdenfile <- "../resultf/burdentest/variant_level.burden.txt"
qwt(variantTable,vburdenfile,flag=2)


## =========================output indels for IGV and variant tables=========================================
varTable <- read.delim("../resultf/burdentest_PP2_hotspots/variant_level.burden.txt")
genes <- union(TSg,union(DRg,DNAreg))
variantlist <- caselist[caselist[,"VariantClass"] %in% c(mis,lof) & caselist[,"Gene"] %in% genes, ]
vars <- paste(variantlist[,1],variantlist[,2],variantlist[,4],variantlist[,5],sep="_")
n.var <- length(vars)
rowTitle <- matrix(0,n.var,8)
rowTitle[variantlist[,"Gene"] %in% TSg,1] <- "Tumor_suppressor"
rowTitle[variantlist[,"Gene"] %in% DRg,1] <- "Cancer_driver"
rowTitle[variantlist[,"Gene"] %in% DNAreg,1] <- "DNA_repair"
rowTitle[variantlist[,"VariantClass"] %in% lof,2] <- "LOF" 
rowTitle[variantlist[,"VariantClass"] %in% mis & nchar(variantlist[,"REF"]) != nchar(variantlist[,"ALT"]),2] <- "indels"
rowTitle[variantlist[,"VariantClass"] %in% mis & nchar(variantlist[,"REF"]) == nchar(variantlist[,"ALT"]),2] <- "d-mis"
rowTitle[,3] <- variantlist[,"Gene"]
rowTitle[,4:8] <- as.matrix(varTable[match(vars,varTable[,2]),c(2:4,7,8)])
variantlist <- variantlist[,colnames(variantlist)!="Gene"]
variantlist <- cbind(rowTitle,variantlist)
variantlist <- variantlist[order(variantlist[,1],variantlist[,2],-as.numeric(variantlist[,7]),as.numeric(variantlist[,8])),]
colnames(variantlist)[1:8] <- c("Gene_sets","variant_type","Gene","Variant","#cases","#controls","Folds","pvalue")
qwt(variantlist,file="../resultf/burdentest/geneset_variantlist.txt",flag=2)
indelsfile <- "../resultf/burdentest/indels.txt"
qwt(variantlist[variantlist[,"variant_type"]=="indels",c("Chromosome","Position","Subject_ID")],file=indelsfile)




