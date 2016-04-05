## ===============================================================================
## burden test
## ===============================================================================
source("misc.R")
source("sourcefiles.R")

BreastCancerAna <- function(sig=TRUE,Ecut=0.001,hotf=1,swi=1,ana=4){
### sig: test for singleton variant only or not
### Ecut: ExAC frequency cutoff
### hotf: 1: hotspots predicted based on HMM (hongjian); 2: hotspots defined by cosmic site with recurrent >=3; 3: not filter missense by hotspot; 4: not filter missense byhotspot but with the hotspot labels (by HMM hongjian)
### swi: 1: Jewish case and control; 2: Hispanic case and contorl
### ana: 1-single gene test only; 2-single variant test only; 3-gene set test only; 4-all related tests; 5-single gene and single variant; 6-single gene and gene sets; 7-single variant and gene sets;
    
## =========================input files===========================================
if(hotf==1){ hotspotfile <- hotHMM; ## HongjianPred hotspots file
}else if(hotf==2){ hotspotfile <- hotCOSMIC; ## COSMIC hotspots file
}else if(hotf==3){ hotspotfile="";
}else if(hotf==4){ hotspotfile <- hotHMM;}
## switch to the right files for burden test
vburdenfile <- vburdenfiles[swi]  ### single variant test files with all annotations; look src.R
alleleFrefile <- alleleFrefiles[swi] 
caselistf <- caselistfs[swi]
contlistf <- contlistfs[swi]
casestaf <- casestafs[swi]   ##using the variant from AJ case and control joint calling result
contstaf <- contstafs[swi]

##====================define the output files===================================================================
pop=c("Jewish","Hispanic")
outputpath <- paste("../resultf/burdentest_",sig,"_singleton_",Ecut,"_ExACFre_",gsub(".txt","",basename(hotspotfile)),"_",hotf,"_",pop[swi],"_",format(Sys.Date(),format="%B_%d_%Y"),"/",sep="")
genesetsVarfile <- paste(outputpath,"Panel_genes_variantlist.txt",sep="")
restVarfile <- paste(outputpath,"Variantlist.txt",sep="")
NumOfVarfile <- paste(outputpath,"NumOfVariant.txt",sep="")
if(!file.exists(outputpath)){ dir.create(outputpath, showWarnings = TRUE, recursive = FALSE);}

## =========================variant class definition and gene sets==========================================
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
vartypes <- list(singleLOF,lof,mis,indel,c(singleLOF,lof),c(mis,indel),NULL,syn,unknown)
vartypenames <- c("singleLOF","indelLOF","singleMIS","indelMIS","LOF","MIS","ALL variants","Synonymous","Unknown")

tmp <- test_genesets()
genesets <- tmp$genesets
genesetnames <- tmp$genesetnames

## =========================variant list filtering===========================================
##getVariantlist(HIcasevariantpath,HIBRfile,namestr=".AllVariants.tsv","../data/Rdata/HIcaselist_1_29")
##getVariantlist(HIcontvariantpath,HIcontrolfile,namestr=".tsv","../data/Rdata/HIcontlist_1_5")
### AJ lists
##getVariantlist(AJcasevariantpath,AJBRfile,namestr=".tsv","../data/Rdata/AJcaselistBR_12_17")
##getVariantlist(AJcontvariantpath,AJBRfile,namestr=".tsv","../data/Rdata/AJcaselist715_12_17")
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
## get control variant list
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

### =========================burden test for single gene, single variant, gene sets=========================================
print_log("burden test for gene, variant, gene sets ...")
## eliminate batch effect for HISP based on rare synonymous
if(swi==1){coe=1;}else if(swi==2){ coe = coeHisp(caselist,contlist);}
print("Corrected coe: ");
print(coe);

##======gene set level burden test=============================================
if(ana==3 | ana==4 | ana==6 | ana==7){
        Gtop <- read.table(GTExfile)
        allgenes <- union(Gtop[,1],unlist(genesets))
        Ggs <- c("top 25%","top 50%","top 75%","top 100%")
        GgL <- list(Gtop[1:floor(0.25*dim(Gtop)[1]),1], Gtop[1:floor(0.5*dim(Gtop)[1]),1], Gtop[1:floor(0.75*dim(Gtop)[1]),1],allgenes)
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
        print_log("Output burden files for gene sets ...")
        outputpath1 <- paste(outputpath,"burden/",sep="")
        if(!file.exists(outputpath1)){ dir.create(outputpath1, showWarnings = TRUE, recursive = FALSE);}
        setburdenfile <- paste(outputpath1,"gene_variant_set.burden.txt",sep="") ### gene set burden test file
        qwt(setburdens,setburdenfile,flag=2)
}

##=====single gene level burden test===========================================
if(ana==1 | ana==4 | ana==5 | ana==6){
        tablelist <- list()
        for(i in 1:length(vartypes)){
                oneTable <- burden_test(caselist,contlist,testset=NULL,testtype=vartypes[[i]],flag=2,sig=sig,coe=coe)
                oneTable <- oneTable[order(as.numeric(oneTable[,"Pvalue"])), ]
                tablelist[[i]] <- oneTable
                onefile <- paste(outputpath1,vartypenames[i],"_level.burden.txt",sep="")
                qwt(tablelist[[i]],onefile,flag=2)
        }
}

##======single variant level burden test=======================================
if(ana==2 | ana==4 | ana==5 | ana==7){
        varSelected <- c(stopins,splices,singleLOF,lof,mis,indel) ## slected variants
        cols0 <- colnames(caselist)
        Vlist <- caselist[caselist[,"VariantClass"] %in% varSelected, ]
        vars <- paste(Vlist[,1],Vlist[,2],Vlist[,4],Vlist[,5],sep="_")
        univar <- unique(vars)
        Vlist <- Vlist[match(univar,vars),]
        
        Vlist[,"Gene_sets"] <- ""
        for(i in 1:length(genesets)){
                Vlist[Vlist[,"Gene"] %in% genesets[[i]],"Gene_sets"] <- genesetnames[i]  
        }
        Vlist[Vlist[,"Gene_sets"]=="","Gene_sets"] <- "non-known_genesets"
        
        Vlist[,"variant_type"] <- ""
        NumVar <- 1:4
        for(i in 1:4){ ### manually changes
                Vlist[Vlist[,"VariantClass"] %in% vartypes[[i]],"variant_type"] <- vartypenames[i]
                NumVar[i] <- length(univar[Vlist[,"VariantClass"] %in% vartypes[[i]]])
        }
        varTable <- read.delim(vburdenfile,check.names=FALSE)
        colnames(varTable)[colnames(varTable)=="Folds"] <- "Odds" ### fix this by update single variant test in src.R 
        Vlist <- cbind(Vlist,varTable[match(univar,varTable[,"Variant"]) , setdiff(colnames(varTable),colnames(Vlist)) ])
        Vlist[,"CorrectedP"] <- 1 # corrected Pvalue by each variant type
        for(i in 1:4){
                Vlist[Vlist[,"variant_type"] == vartypenames[i],"CorrectedP"] <-  p.adjust(as.numeric(Vlist[Vlist[,"variant_type"] == vartypenames[i],"Pvalue"]),method = "BH")
        }
        NumVar <- paste(vartypenames[1:4]," : ",NumVar,sep="")
        
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
                if(swi==1){ cols <- c("Gene_sets","variant_type","Gene","Variant","#in_case","#in_cont","Odds","Pvalue","CorrectedP",setdiff(AJcols,c("N.index.AJ","N.cont.AJ","Odds_AJ_cont","p_AJ_cont")),HIcols,setdiff(cols0,"Gene"));}
                if(swi==2){ cols <- c("Gene_sets","variant_type","Gene","Variant","#in_case","#in_cont","Odds","Pvalue","CorrectedP",setdiff(HIcols,c("N.index.HI","N.cont.HI","Odds_HI_cont","p_HI_cont")),AJcols,setdiff(cols0,"Gene"));}
                Vlist <- Vlist[,cols]
                ### selected columns and re-orders
                Vlist <- Vlist[as.numeric(Vlist[,"Odds"])>1, ] ### keep less variants
                Vlist <- Vlist[order(as.numeric(Vlist[,"Pvalue"])),]
        }
        
        ### label Dominicans with EU-ASN
        if(swi==2){
                domis <- read.table(SubDominif)
                euafr <- domis[domis[,2]=="EU-AFR",1]
                euafs <- paste(euafr,sep="",collapse="|")
                EUASN <- sapply(1:dim(Vlist)[1],function(kk){
                        c(gsub(euafs,"", Vlist[kk,"index.HI"]), gsub(euafs,"", Vlist[kk,"pseudoCont.HI"]),gsub(euafs,"", Vlist[kk,"case.HI"]),gsub(euafs,"", Vlist[kk,"non_case.HI"]))
                })
                EUASN <- t(EUASN)
                colnames(EUASN) <- c("index.EU.ASN","pseudoCont.EU.ASN","case.EU.ASN","non_case.EU.ASN")
                Vlist <- cbind(Vlist,EUASN)
        }
        
        ## write test result files
        qwt(Vlist[Vlist[,"Gene_sets"]!= "non-known_genesets", ],file=genesetsVarfile,flag=2)
        qwt(Vlist[Vlist[,"Gene_sets"]== "non-known_genesets", ],file=restVarfile,flag=2)
        qwt(NumVar,file=NumOfVarfile)
        
        ## QQ plot for not singleton variant tests
        if(!sig){ 
                ### only for QQ plots
                #varT <- list(lof,mis,indel,syn,stopins,splices,singleLOF)
                #varTnames <- c("lof_indel","d-mis","mis-indel","Synonymous","Stopgain_loss","splicing","singleLOF")
                qqplot_variantLevel(vburdenfile,paste(outputpath,"qqplots/",sep=""),unique(unlist(genesets)),vartypes,vartypenames,caselist);
        }
        
}

##========end of test=============
}
