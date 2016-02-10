genelists <- function(){
## tumor suppressors 
ts2  <- unlist(read.table("../../genelist/Tumor_supressor/TS_Cosmic.txt")) ## tumor supressors
ts3 <- unlist(read.table("../data/hotspots/TScell_filtered.txt"))
ts <- union(ts2,ts3)
write.table(ts,file="../data/hotspots/TS_collected.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

ts <- gsub("Sep-05","SEPT5",ts);ts <- gsub("Sep-06","SEPT6",ts);ts <- gsub("Sep-09","SEPT9",ts);
qwt(ts,"../data/hotspots/Tumor_suppressors_11_11.txt")

}

variantDis <- function(){
    source("misc.R")
    source("src.R")
    load("../data/Rdata/AJcaselist_11_9")
    phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv" ## phenotype file to get index cases only
    outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"  ## outlier samples
    caselist <- onelist
    rm(onelist)
    indexcases <- getindexcase(phenofile)
    caselist <- caselist[caselist[,"Subject_ID"] %in% indexcases, ] ## only index cases
    outliers <- unlist(read.table(outlierfile)) ## exclude outlier subjects
    outliers <- sapply(1:length(outliers),function(i) {tmp=unlist(strsplit(outliers[i],"_"));setdiff(gsub("\\D","",tmp),c("",paste("00",1:9,sep="")));})
    caselist <- caselist[!(caselist[,"Subject_ID"] %in% outliers), ]
    load("../data/Rdata/AJcontlist_11_9")
    contlist <- onelist
    varTypes <- variantDis_one(caselist,contlist)
    
    rm(caselist)
    rm(contlist)
    
    load("../resultf/caselist_singleton_0.01")
    load("../resultf/contlist_singleton_0.01")
    varTypes1 <- variantDis_one(caselist,contlist)
    
    #=======================================
    outputpath <- "../resultf/"
    qwt(varTypes,file=paste(outputpath,"Jewish_case_control_variants.txt",sep=""),flag=2)
    qwt(varTypes1,file=paste(outputpath,"Jewish_case_control_variants_filtered_singleton.txt",sep=""),flag=2)
 
    density_plots(varTypes[,3:7],varTypes[,1],"../resultf/plots_all/")
    density_plots(varTypes1[,3:7],varTypes1[,1],"../resultf/plots_singleton_0.01/")
    
}

variantDis_one <- function(caselist,contlist){

    VariantClass <- c(".","frameshiftdeletion","frameshiftinsertion","none","nonframeshiftdeletion","nonframeshiftinsertion","nonsynonymousSNV","stopgain","stoploss","synonymousSNV","unknown")
    lof <- c("frameshiftdeletion","frameshiftinsertion","none","stopgain","stoploss",".")
    mis <- "nonsynonymousSNV"
    indel <- c("nonframeshiftdeletion","nonframeshiftinsertion")
    syn <- "synonymousSNV"
    unknown <- c("unknown")
    
    cases <- unique(caselist[,"Subject_ID"])
    conts <- unique(contlist[,"Subject_ID"])
    n.case <- length(cases)
    n.cont <- length(conts)
    varTypes <- matrix(0,n.case+n.cont,7)
    
    for(i in 1:n.case){
        onecase <- caselist[caselist[,"Subject_ID"]==cases[i], ]
        varTypes[i,3:7] <- c(sum(onecase[,"VariantClass"] %in% lof), sum(onecase[,"VariantClass"] %in% mis), sum(onecase[,"VariantClass"] %in% indel),sum(onecase[,"VariantClass"] %in% syn),sum(onecase[,"VariantClass"] %in% unknown))
    }
    for(i in 1:n.cont){
        onecont <- contlist[contlist[,"Subject_ID"]==conts[i], ]
        varTypes[i+n.case,3:7] <- c(sum(onecont[,"VariantClass"] %in% lof), sum(onecont[,"VariantClass"] %in% mis), sum(onecase[,"VariantClass"] %in% indel), sum(onecont[,"VariantClass"] %in% syn),sum(onecont[,"VariantClass"] %in% unknown))
    }
    
    varTypes[,1] <- "control"
    varTypes[1:n.case,1] <- "case"
    varTypes[1:n.case,2] <- cases
    varTypes[(n.case+1):(n.case+n.cont),2] <- conts
    colnames(varTypes) <- c("case/control","Subject_ID","#LOF","#MIS","#indel","#synonymous","#unknown")
    
    varTypes  
}

variantLargeFamily <- function(){
    # low in control data, high in family cases and low in unaffected memebers
    pedigree <- read.csv("/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/ALL_pedigree.csv")
    onelist <- read.delim("../resultf/LargeFamilyVariants.txt")
    pheno <- read.delim("../resultf/LargeFamilyPhenotype.txt")
    
    subAJ <- pheno[pheno[,"AJFAM"]=="J",3]
    contlistf <- "../data/Rdata/AJcontlist_11_9"
    wfile="../resultf/LargeFamilyJewish.txt"
    Variantcounting(pheno,subAJ,onelist,contlistf,wfile)
    
    subAJ <- pheno[pheno[,"HISPFAM"]=="H",3]
    contlistf <- "../data/Rdata/HIcontlist_11_20"
    wfile="../resultf/LargeFamilyHispanic.txt"
    Variantcounting(pheno,subAJ,onelist,contlistf,wfile)

}

Variantcounting <- function(pheno,subAJ,onelist,contlistf,wfile){
    source("misc.R")
    subAJf <- subAJ[pheno[match(subAJ,pheno[,3]),"BreastCancer"]=="Yes"]
    subAJu <- subAJ[pheno[match(subAJ,pheno[,3]),"BreastCancer"]=="No"]
    caseAJf <- onelist[onelist[,"Subject_ID"] %in% subAJf,]
    caseAJu <- onelist[onelist[,"Subject_ID"] %in% subAJu,]
    load(contlistf)
    contlist <- onelist
    
    fvar <- paste(caseAJf[,1],caseAJf[,2],caseAJf[,4],caseAJf[,5],sep="_")
    uvar <- paste(caseAJu[,1],caseAJu[,2],caseAJu[,4],caseAJu[,5],sep="_")
    convar <- paste(contlist[,1],contlist[,2],contlist[,4],contlist[,5],sep="_")
    
    ffvar <- table(fvar)
    fuvar <- table(uvar)
    fcvar <- table(convar)
    
    utmp <- setdiff(names(fuvar),names(ffvar))
    
    atmp <- setdiff(union(names(ffvar),names(fuvar)),names(fcvar))
    fcvar[atmp] <- 0 
    atmp <- setdiff(names(ffvar),names(fuvar))
    fuvar[atmp] <- 0
    
    
    unifvar <- unique(fvar)
    s <- (ffvar[unifvar]+1)/(fuvar[unifvar]+2*fcvar[unifvar]+1)
    s1 <- - fuvar[utmp] 
    aa <- cbind(c(unifvar,utmp),c(s,s1))  
    aa <- cbind(aa,0,0,0)
    n1 <- length(unifvar)
    n2 <- length(c(unifvar,utmp))
    aa[1:n1,3] <- ffvar[unifvar]
    aa[1:n1,4] <- fuvar[unifvar]
    aa[1:n1,5] <- fcvar[unifvar]
    aa[(n1+1):n2,4] <- fuvar[utmp]
    aa[(n1+1):n2,5] <- fcvar[utmp]

    wvar <- rbind(caseAJf,caseAJu)
    vars <- paste(wvar[,1],wvar[,2],wvar[,4],wvar[,5],sep="_")
    wvar <- cbind(aa[match(vars,aa[,1]),2:5],wvar)
    wvar <- wvar[order(-as.numeric(wvar[,1])),]
    colnames(wvar)[1:4] <- c("score","affect","unaffect","control")
    wvar <- wvar[,!(colnames(wvar) %in% c("filtered","ExACfreq","VCFPASS","noneSegmentalDup","meta.SVM_PP2","hotspot","alleleFre"))]
    wvar <- wvar[,c(dim(wvar)[2],1:(dim(wvar)[2]-1))]
    wvar <- cbind(pheno[match(wvar[,"Subject_ID"],pheno[,3]),c("FAMILYID","BreastCancer")],wvar)
    colnames(wvar)[1] <- "FamilyID"
    
    qwt(wvar,file=wfile,flag=2)

}

doublecheck <- function(){
    source("misc.R")
    source("src.R")
    varfile <- "../resultf/burdentest_FALSE_0.01_HMM_hotspots_11_12/Panel_genes_variantlist.txt"
    onevar <- read.delim(varfile)
    onevar <- onevar[onevar[,"Gene"]=="MSH3",]
    
    ## write indels
    oneids <- c()
    for(i in 2:5){
        tmp <- unlist(strsplit(onevar[1,i],":"))
        oneids <- union(oneids,tmp)
    }
    onevarf <- matrix(0,length(oneids),3)
    onevarf[,1] <- onevar[1,"Chromosome"]
    onevarf[,2] <- onevar[1,"Position"]
    onevarf[,3] <- oneids
    qwt(onevarf,file="../single_check/MSH3_IGV.txt")
        
    pheno <- phenoinfo()
    onephe <- pheno[pheno[,3] %in% oneids,]
    qwt(onephe,file="../single_check/MSH3_pheno.txt",flag=2)
    
    onephe <- pheno[pheno[,3] %in% oneids & pheno[,"BreastCancer"]=="Yes",]
    qwt(onephe,file="../single_check/MSH3affected_pheno.txt",flag=2)
    
    pedigree <- read.csv("/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/ALL_pedigree.csv") # ALLadd_pedigree.csv
    indids <- pheno[pheno[,3] %in% oneids,2]
    oneped <- pedigree[pedigree[,2] %in% indids | pedigree[,3] %in% indids | pedigree[,4] %in% indids,]
    oneped <- cbind(oneped,oneped[,2:4])
    oneped[oneped[,5] %in% pheno[,2],5] <- pheno[match(oneped[oneped[,5] %in% pheno[,2],5],pheno[,2]), 3]
    oneped[oneped[,6] %in% pheno[,2],6] <- pheno[match(oneped[oneped[,6] %in% pheno[,2],6],pheno[,2]), 3]
    oneped[oneped[,7] %in% pheno[,2],7] <- pheno[match(oneped[oneped[,7] %in% pheno[,2],7],pheno[,2]), 3]
    qwt(oneped,file="../single_check/MSH3_pedigree.txt",flag=2)
    
    
    pheno <- phenoinfo()
    affs <- pheno[pheno[,"BreastCancer"]=="Yes",3]
    uffs <- pheno[pheno[,"BreastCancer"]=="No",3]
    AJcasefile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/tablesandlists/All_AJ_samples.list"
    HIcasefile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/HispanicCases548.txt"
    Ajs <- unlist(read.table(AJcasefile))
    His <- unlist(read.table(HIcasefile))
    
    fams <- pheno[match(oneids,pheno[,3]),1]
    allmems <- pheno[pheno[,1] %in% fams,3]
    tyvar <- read.delim("../single_check/MSH3_detail.txt",header=FALSE)
    tys <- c("ins","substitute","complex","substitute/del")
    numT <- matrix(0,4,4)
    conl <- setdiff(allmems,oneids)
    conNum <- c(sum(conl %in% intersect(affs,Ajs)),sum(conl %in% intersect(uffs,Ajs)),sum(conl %in% intersect(affs,His)),sum(conl %in% intersect(uffs,His)))
    for (i in 1:4){
        numT[i,] <- c(sum(oneids %in% intersect(intersect(affs,Ajs),tyvar[tyvar[,4] %in% tys[i],3])),sum(oneids %in% intersect(intersect(uffs,Ajs),tyvar[tyvar[,4] %in% tys[i],3])),sum(oneids %in% intersect(intersect(affs,His),tyvar[tyvar[,4] %in% tys[i],3])),sum(oneids %in% intersect(intersect(uffs,His),tyvar[tyvar[,4] %in% tys[i],3])) )
    }
    
}

depth_variant <- function(){

	AJindexf <- "/home/local/ARCS/qh2159/breast_cancer/variants/AJdepth/indexAJ.ldepth.mean"
	AJcontf <- "/home/local/ARCS/qh2159/breast_cancer/variants/AJdepth/AJ.ldepth.mean"

	HIindexf <- "/home/local/ARCS/qh2159/breast_cancer/variants/HIdepth/indexHI.ldepth.mean"
	HIcontf <- "/home/local/ARCS/qh2159/breast_cancer/variants/HIdepth/HI.ldepth.mean"
	
	onedepthcompare(AJindexf,AJcontf)
	onedepthcompare(HIindexf,HIcontf)
}

onedepthcompare <- function(AJindexf,AJcontf){

  	indexf <- read.delim(AJindexf)
        contf  <- read.delim(AJcontf)
        caseV <- paste(indexf[,1],indexf[,2],sep="_")
        contV <- paste(contf[,1],contf[,2],sep="_")


	interV <- intersect(caseV,contV)
	unicase <- indexf[!(caseV %in% interV),]
	unicont <- contf[!(contV %in% interV),]


        chr=c(1:22,"X")
        for(i in 1:length(chr)){

                a1= unicase[unicase[,1]==chr[i],2]
                a2= unicont[unicont[,1]==chr[i],2]

                a1M <- matrix(a1,ncol=length(a2),nrow=length(a1),byrow=FALSE)
                a2M <- matrix(a2,ncol=length(a2),nrow=length(a1),byrow=TRUE)
                print(sum(abs(a1M-a2M)<=30))

        }
}

twoVariantSta <- function(){

source("misc.R")
outputpath="../resultf/variantSta_badInter/"
contf <- "../data/contAJ_20151209.hardfiltered.stats_hq.tsv"
casef <- contf
cases <- unlist(read.table("../data/AJindexcases265.txt"))
controls <- unlist(read.table("/home/local/ARCS/qh2159/breast_cancer/variants/data/AJs_557.txt"))

    if(!file.exists(outputpath)){ dir.create(outputpath, showWarnings = TRUE, recursive = FALSE);}
    varT=c("synonymous","Missense","indels","Frameshift","ALL")
    
    tmp <- read.table(casef,fill=T)
    if(all(cases %in% tmp[2,])){
        subs <- match(cases,tmp[2,])
    }else{
        subs <- sapply(1:length(cases),function(i) which(grepl(cases[i],tmp[2,])))
    }
    varTa <- t(tmp[c(21,22,4,25,3),subs])
    varTa <- cbind("case",cases,varTa)
    colnames(varTa) <- c("Group","Subject_ID",varT)
    tmp <- read.table(contf,fill=T)
    varTa1 <- t(tmp[c(21,22,4,25,3),match(controls,tmp[2,])])
    varTa1 <- cbind("control",controls,varTa1)
    colnames(varTa1) <- c("Group","Subject_ID",varT)
    Ta <- rbind(varTa,varTa1)
    n <- dim(Ta)[2]
    Ta0 <- Ta

    contf <- "../data/contAJ_20151215.badInter.stats_hq.tsv"
    casef <- contf
    tmp <- read.table(casef,fill=T)    
    if(all(cases %in% tmp[2,])){
        subs <- match(cases,tmp[2,])
    }else{
        subs <- sapply(1:length(cases),function(i) which(grepl(cases[i],tmp[2,])))
    }
    varTa <- t(tmp[c(21,22,4,25,3),subs])
    varTa <- cbind("case",cases,varTa)
    colnames(varTa) <- c("Group","Subject_ID",varT)
    tmp <- read.table(contf,fill=T)
    varTa1 <- t(tmp[c(21,22,4,25,3),match(controls,tmp[2,])])
    varTa1 <- cbind("control",controls,varTa1)
    colnames(varTa1) <- c("Group","Subject_ID",varT)
    Ta <- rbind(varTa,varTa1)

    a1 = Ta0[,3:7]
    mode(a1) <- "numeric"
    a2 = Ta[,3:7]
    mode(a2) <- "numeric"
    Ta[,3:7] <- a1 - a2   
 
    density_plots(Ta[,3:n],Ta[,1],outputpath,VarT=varT)

}

interested_genes <- function(){
    source("src.R")
    source("misc.R")
    source("srcp.R")
    AJf <- "../resultf/variant_level_burden_anno_Fre_Pseducont.txt"
    HIf <- "../resultf/HISP_variant_level_burden_anno_Fre_Pseducont.txt"
    HIg <- "../resultf/Panel_genes_Potential.txt"
    vars <- read.delim(HIf)
    
vars <- vars[vars[,"VariantClass"] %in% c("frameshiftdeletion","frameshiftinsertion","nonsynonymousSNV","nonframeshiftdeletion","nonframeshiftinsertion","stopgain","stoploss","none","."), ]

#genes <- c("MSH3","PARP4","PTPRF","ARID1B","POT1","SETD2","FBXW7","HOXD11")
#genes <- c("BTN3A3","C20orf96","SLC34A2","PTPRF","VAMP5","SHROOM4","RPGR","HLA-A","MTR","PAPLN","NPEPL1","CNOT1","HLA-C","SCYL1","IFI35","ARFGAP3","IFNA7")
#genes <- c("CYP4A22","PER3","CNOT1","AGMO","PCDHGB3","CTBS","GIGYF2","DCAF8L2","LRP1B","SDK1")
#genes <- c("CCDC80","TOM1","NRG1","NEURL4","NLRP14","RAB11FIP5")

tmpT <- read.delim(HIg)

# testg <- c()
# testv <- c()
# for(i in 1:length(genes)){
# 	tmps <- which( vars[,"Gene"]==genes[i] & as.numeric(vars[,"Pvalue"]) < 0.05 )  ##!!!!!!!
# 	testv <- c(testv,paste(vars[tmps,"Chromosome"],vars[tmps,"Position"],sep="_"))
# 	testg <- c(testg,rep(genes[i],length(tmps)))
# }
#testv[2] <- "16_58577315_GAA_GAAAA"

testg <- tmpT[,"Gene"]
testv <- paste(tmpT[,"Chromosome"],tmpT[,"Position"],sep="_")

Pv <- sapply(1:length(testv),function(i) unlist(strsplit(testv[i],"_"))[1:2])
Pv <- t(Pv)
Pv <- paste(Pv[,1],":",as.numeric(Pv[,2]),"-",as.numeric(Pv[,2]),sep="")
Pv <- cbind(Pv, testg)
qwt(Pv,file=paste("../single_check/variant_table",format(Sys.Date(),format="%B_%d_%Y"),".txt",sep=""))

phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv"
load("../data/Rdata/BreastCancer_VariantList_11_12")
pheno <- phenoinfo()
#onevar <- paste(onelist[,1],onelist[,2],onelist[,4],onelist[,5],sep="_")
onevar <- paste(onelist[,1],onelist[,2],sep="_")
for(i in 1:length(testv)){
	oneVariant(testv[i],testg[i],onelist,onevar,pheno)
}

}

interested_genes_inFam <- function(){
    
    Cohortfile <- "../data/Rdata/BreastCancer_VariantList_11_12"
    load(Cohortfile)
    source("src.R")
    pheno <- phenoinfo()
    
    genes <- matrix(c(13,25077802,"PARP4",7,124482897,"POT1",3,47162897,"SETD2",2,176972404,"HOXD11",4,153268224,"FBXW7"),ncol=3,byrow=TRUE)
    genes <- matrix(c(1,7895966,"PER3"),ncol=3,byrow=TRUE)
    for(i in 1:dim(genes)[1]){
        subs <- onelist[,1]==genes[i,1] & onelist[,2]==genes[i,2] & onelist[,"Gene"]==genes[i,3]
        sams <- unique(onelist[subs,"Subject_ID"])
        onephe <- pheno[pheno[,3] %in% sams,]
        qwt(onephe,file=paste("../single_check/Pheno_",genes[i,3],".txt",sep=""),flag=2)
    }
    
    
}

oneVariant <- function(var,gene,onelist,onevar,pheno){
	pops <- paste(pheno[,"AJFAM"],pheno[,"HISPFAM"],sep="")
	samids <- onelist[onevar==var,"Subject_ID"]
	famids <- pheno[match(samids,pheno[,3]),1]
	indids <- pheno[match(samids,pheno[,3]),2]
	popids <- pops[match(samids,pheno[,3])]
	affids <- pheno[match(samids,pheno[,3]),"BreastCancer"]

	oneC <- matrix(0,length(famids),4)
	for(i in 1:length(famids)){
		onefams <- pheno[pheno[,1]==famids[i],3]
		oneC[i,1] <- sum((samids %in% onefams) & affids=="Yes")
		oneC[i,2] <- sum((samids %in% onefams) & affids=="No")
		oneC[i,3] <- sum(pheno[match(onefams,pheno[,3]),"BreastCancer"]=="Yes")
		oneC[i,4] <- sum(pheno[match(onefams,pheno[,3]),"BreastCancer"]=="No")
	}

	oneT <- cbind(famids,indids,samids,popids,affids,oneC)
	oneT[is.na(oneT)] <- ""
	colnames(oneT) <- c("FAMILYID","INDIVIDUALID","SUBJECTID","POP","AFFSTATUS","#case-with","#control-with","#CASE","#CONTROL")
	qwt(oneT,file=paste("../single_check/",gene,"_",var,".txt",sep=""),flag=2)

	## output IGV indels
	tmp <- unlist(strsplit(var,"_"))
	indf <- cbind(tmp[1],tmp[2],samids)
	qwt(indf,file=paste("../single_check/",gene,"_",var,"_IGVs.txt",sep=""))
}

PARP4 <- function(){
	source("src.R")
	pheno <- phenoinfo()
	PARP4sam <- unlist(read.table("../single_check/PARP4.txt"))
	PARP4Phe <- pheno[pheno[,3] %in% PARP4sam,]
	PARP4Phe[is.na(PARP4Phe)] <- ""
	qwt(PARP4Phe,file="../single_check/Pheno_PARP4.txt",flag=2)	
}

filtered_LargeFam <- function(){
	source("misc.R")
	mis <- "nonsynonymousSNV"
	#path= "/home/local/ARCS/qh2159/breast_cancer/variants/families/PriorFamilies29" 
	#path= "/home/local/ARCS/qh2159/breast_cancer/variants/families/L2Families35"
	path="/home/local/ARCS/qh2159/breast_cancer/variants/families/Families78V2"
	## second 5 families
	#sams <- c("200500","300424","200124","200396","200154")
	#sams <- "200352"
	#samf = c(paste(path,"/",sams,".AD.tsv",sep=""), paste(path,"/",sams,".AR.tsv",sep=""))
	## first 9 families
	#samf  <- list.files(path=path,pattern=".tsv$",full.names=TRUE)
	## Prior 29 families
	samf  <- list.files(path=path,pattern=".tsv$",full.names=TRUE)
	for(i in 1:length(samf)){
        	oner <- read.delim(samf[i],check.names=FALSE)
		if(dim(oner)[1]>0){
    		oner <- variant_filtering(oner,mis,Ecut=0.01,segd=0.95,pp2=TRUE,hotf="",alleleFrefile="",popcut=0.05)
		oner <- oner[oner[,"ExACfreq"] & oner[,"VCFPASS"] & oner[,"noneSegmentalDup"],]
		### further filtered by others
		## meta-svm D or CADD >= 15 or polyphen as D
		subs <- rep(TRUE,dim(oner)[1])
		subs[oner[,"VariantClass"] %in% mis] <- FALSE
		oner <- oner[subs | (oner[,"PP2prediction"]=="D" | oner[,"MetaSVM"]=="D" | oner[,"CADDscore"] >= 15), ]
    		qwt(oner,file=gsub(".tsv","filtered.tsv",samf[i]),flag=2)
		}
	}
}

fold1_2Sam <- function(){

gg <- c()
for(c in seq(100,1000,100)){
for(d in seq(100,3000,100))
{
for(b in seq(50,floor(0.5*d),50)){
a <- floor(1.2*(b/d)*c)
if(a < 0.3*c){
tmp <- fisher.test(matrix(c(a,b,c-a,d-b),2,2))$p.value
tmp <- c(a,b,c,d,tmp)
gg <- rbind(gg,tmp)
}
}
}
}
colnames(gg) <- c("carrier","non_carrier","cases","controls","pvalue")
gg <- gg[gg[,5] < 0.05, ]
qwt(gg,file="Sign_1.2.txt",flag=2)


gg <- c()
c=265
a=82
for(d in seq(200,7200,10))
{
b <-  floor(d*(a/c)/1.2)
tmp <- fisher.test(matrix(c(a,b,c-a,d-b),2,2))$p.value
tmp <- c(a,b,c,d,tmp)
gg <- rbind(gg,tmp)
}
plot(seq(200,7200,10),gg[,5],col=2,type='l',xlab="Number of controls",ylab="Pvalue")
abline(h=0.05,lwd=2)

}

### coding variant numbers in four batches of Hispnic 
HispStaBatches <- function(){

phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv" 
HIcasefile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/HispanicCases548.txt"
batches <- "../data/phenotype/Sample_batches.csv"
HIindexfile <- "../data/HIindexcases138.txt"


HIs <- unlist(read.table(HIcasefile))
countT <- read.delim("../data/BR.origin.stats_12_28.tsv",header=FALSE)
bats <- read.csv(batches)

batlab <- unique(bats[,3])
cots <- list()
for(i in 1:4){
	cots[[i]] <- as.numeric(countT[3,match(HIs[HIs %in% bats[bats[,3]==batlab[i],2]],countT[2,])])
	print(length(cots[[i]]))
}

pdf(file="../resultf/HISPSampels.pdf",height=10,width=12)
plot(density(cots[[2]]),xlab="Number of coding variants",main="Coding variants",cex=1.5,col=1,type="l")
dev.off()


countT <- read.delim("../data/BR.origin.stats_12_28.tsv",header=FALSE)
dupsam <- countT[2,which(grepl("cgc",countT[2,]))]
dupsam1 <- gsub("cgc","",dupsam)
aa <- countT[3,match(dupsam,countT[2,])]
bb <- countT[3,match(dupsam1,countT[2,])]
aa <- as.numeric(aa)
bb <- as.numeric(bb)

pdf(file="../resultf/dupSampels.pdf",height=10,width=12)
plot(density(aa),xlab="Number of coding variants",main="Coding variants in duplication samples",cex=1.5,col=1,type="l",ylim=c(0,0.00045))
lines(density(bb),col=2,type="l")
legend("topright",legend=c("Yale","Regeneron"),col=1:2,lwd=rep(1,2),lty=rep(1,2))
dev.off()


}

### MSH3 variant in subtypes: substitutions; indels:different with length; both;
MSH3details <- function(){
        source("sourcefiles.R")
        source("misc.R")
        load(Cohortfile)
        subs <- onelist[,"Gene"]=="MSH3" & onelist[,"Chromosome"]=="5" & onelist[,"Position"]=="79950727"
        aa <- onelist[subs,]
        
        AJs <- unlist(read.table(AJBRfile))
        HIs <- unlist(read.table(HIBRfile))
        samid <- unique(aa[,"Subject_ID"])
        sum(samid %in% AJs)
        sum(samid %in% HIs)
        length(setdiff(samid,c(AJs,HIs)))
        unid <- setdiff(samid,c(AJs,HIs))

        ## three lists   
        cohorts <- list()
        cohorts[[1]] <- aa[aa[,"Subject_ID"] %in% AJs, ]
        cohorts[[2]] <- aa[aa[,"Subject_ID"] %in% HIs, ]
        cohorts[[3]] <- aa[aa[,"Subject_ID"] %in% unid, ]
        
        pheno <- phenoinfo()
        alts <- c("C","GCAGCGCCCCCAGCGCCCC","GCAGCGCCCCCAGCGCCCCCAGCGCCCC","GCAGCGCCCCCAGCGCCCCCAGCGCCCCCAGCGCCCC","G")
        cases <- pheno[pheno[,"BreastCancer"] == "Yes", 3]
        tmp <- matrix(0,3,5)
        tmp1 <- matrix(0,3,5)
        for(i in 1:3){
                for(j in 1:length(alts)){
                        tmp[i,j] <- length(unique(cohorts[[i]][cohorts[[i]][,"ALT"]==alts[j],"Subject_ID"]))
                        tmp1[i,j] <- sum(unique(cohorts[[i]][cohorts[[i]][,"ALT"]==alts[j],"Subject_ID"]) %in% cases)
                }        
        }
        t(tmp)
        t(tmp1)
        
        ## two control lists
        load(contlistfs[1])
        contAJ <- onelist
        rm(onelist)
        
        load(contlistfs[2])
        contHI <- onelist
        rm(onelist)
        
        conts <- list()
        subs <- contAJ[,"Gene"] == "MSH3" & contAJ[,"Chromosome"] == "5" & contAJ[,"Position"] == "79950727"
        conts[[1]] <- contAJ[subs, ]
        subs <- contHI[,"Gene"] == "MSH3" & contHI[,"Chromosome"] == "5" & contHI[,"Position"] == "79950727"
        conts[[2]] <- contHI[subs, ]
        
        tmp <- matrix(0,2,5)
        for(i in 1:3){
                for(j in 1:length(alts)){
                        tmp[i,j] <- length(unique(conts[[i]][conts[[i]][,"ALT"]==alts[j],"Subject_ID"]))
                }        
        }
        t(tmp)
        
        
        subs <- onelist[,"Gene"]=="MSH3" & onelist[,"Chromosome"]=="5" & onelist[,"Position"]=="79950724"
        bb <- onelist[subs,]
}

MSH3lollipop <- function(){
       source("srcp.R")
       onerCaseAJ <- caselist[caselist[,"Gene"]=="MSH3", ]
       onelollipopf(onerCaseAJ,"AJ-Case")
       oner <- contlist[contlist[,"Gene"]=="MSH3", ]
       onelollipopf(oner,"AJ-Control")
      
       source("srcp.R")
       oner <- caselist[caselist[,"Gene"]=="MSH3", ]
       onelollipopf(oner,"HI-Case")
       oner <- contlist[contlist[,"Gene"]=="MSH3", ]
       onelollipopf(oner,"HI-Control")
     
       
       ### make the colors in cbioportal
       lof <- c("frameshiftdeletion","frameshiftinsertion")
       mis <- "nonsynonymousSNV"
       indel <- c("nonframeshiftdeletion","nonframeshiftinsertion")
       stopins <- c("stopgain","stoploss")
       splices <- c("none",".")
       singleLOF <- c("stopgain","stoploss","none",".")
       
       filenames <- paste("../single_check/MSH3/",c("variant_BRAJ-Case.txt","variant_AJ-Case.txt","variant_AJ-Control.txt","variant_HI-Case.txt","variant_HI-Control.txt"),sep="")
       for(i in 1:length(filenames)){
               aa <- read.delim(filenames[i])
               aa[aa[,"Mutation_Type"] %in% c(mis,indel), "Mutation_Type"] <- "Missense_Mutation"
               aa[aa[,"Mutation_Type"] %in% singleLOF, "Mutation_Type"] <- "Truncating_Mutation"
               aa[aa[,"Mutation_Type"] %in% lof, "Mutation_Type"] <- "Inframe_Mutation"
               qwt(aa,file=gsub(".txt","Col.txt",filenames[i]),flag=2)
       }
}

onelollipopf <- function(onerCaseAJ,onestr){
        cls <- c("Hugo_Symbol","Protein_Change","Sample_ID","Mutation_Type","Chromosome","Start_Position","End_Position","Reference_Allele","Variant_Allele","Validation_Status","Mutation_Status","Center")
        n.cl <- length(cls)
        
        n <- dim(onerCaseAJ)[1]
        tmp1 <- matrix(0,n,n.cl)
        tmp1[,1:(n.cl-3)] <- as.matrix(onerCaseAJ[,c("Gene","AAchange","Subject_ID","VariantClass","Chromosome","Position","Position","REF","ALT")])
        tmp1[,n.cl-2] <- "No"
        tmp1[,n.cl-1] <- "Germline"
        tmp1[,n.cl] <- onestr
        
        for(i in 1:n){
                tmp1[i,2] <- unlist(strsplit(tmp1[i,2],":"))[5]
                }

        colnames(tmp1) <- cls
        qwt(tmp1,file=paste("../single_check/MSH3/variant_",onestr,".txt",sep=""),flag = 2)
}

### known risk gene double check
knownRisks <- function(){
        genes <- c("BLM","PALB2","WRN","CHEK2","ATM")

        AJv <- read.delim("../resultf/AJ_variant_level_burden_Pseducont.txt")
        tmp <- AJv[AJv[,"Gene"] %in% genes & !(AJv[,"VariantClass"] %in% c("synonymousSNV","unknown")),]
        qwt(tmp,file="../single_check/knownRisk_AJ.txt",flag=2)
        
        HIv <- read.delim("../resultf/HISP_variant_level_burden_Pseducont.txt")
        tmp <- HIv[HIv[,"Gene"] %in% genes & !(HIv[,"VariantClass"] %in% c("synonymousSNV","unknown")),]
        qwt(tmp,file="../single_check/knownRisk_HI.txt",flag=2)      
}
