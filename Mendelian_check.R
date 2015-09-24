Mendelian_check <- function(){
    a <- read.csv("Famcheck.csv")
    inde <- nchar(a[,"REF"]) != nchar(a[,"ALT"])
    
    vars <- paste(a[,"Gene"],a[,2],a[,3],a[,6],a[,7],sep="_")
    
    varT <- read.delim("gene_mutations_dp4.txt",sep="\t")
    vall <- paste(varT[,"Gene"],varT[,2],varT[,3],varT[,5],varT[,6],sep="_")
    
    vars[168] <- vall[189]
    
    subs <- match(vars,vall)
    subID <- varT[subs,"Subject_IDs"]
    
    subIDs <- c()
    for(i in 1:length(subID)){
        tmp <- unlist(strsplit(subID[i],","))
        subIDs <- union(subIDs,tmp)
    }
    
    pheno <- read.csv("WES BCFR phenotypic data.csv")
    ## excluded outlier samples and undetermined and mismatched sex samples
    outliers <- read.delim("Potential_Outliers.tsv")[,1]
    sexcheck <- read.delim("CUMC_Regeneron.mismatches.sexcheck.tsv")[,1]
    outlfam <- read.delim("Potential_Problem_families.tsv",sep=" ")[,1]
    outliers <- union(outliers,sexcheck)
    pheno <- pheno[!(pheno[,"Subject_ID"] %in% outliers),]
    pheno <- pheno[!(pheno[,"FAMILYID"] %in% outlfam),]
    
    phsub <- match(subIDs,pheno[,3])
    subpheno <- pheno[phsub,]
    write.csv(subpheno,file="allcheckpheno.csv",row.names=FALSE)
    
    ## after add parents individual ids
    a <- read.csv("allcheckpheno.csv")
    subs <- !is.na(a[,"PID"])
    idt <- a[subs,c("FAMILYID","INDIVID","Subject_ID","PID","MID")]
    PID_ID <- pheno[match(idt[,"PID"],pheno[,2]),3]
    MID_ID <- pheno[match(idt[,"MID"],pheno[,2]),3]

    idt <- cbind(idt,PID_ID,MID_ID)
    write.csv(idt,file="Famcheckedids.csv",row.names=FALSE)
    
    
    subs <- !is.na(idt[,"PID_ID"]) & !is.na(idt[,"MID_ID"])
    mencase <- idt[subs,"Subject_ID"]
    
    subs <-c()
    for(i in 1:length(mencase)){
        subs <- union(subs,which(grepl(mencase[i],varT[,"Subject_IDs"])))
    }
    
    menvar <- varT[subs,]
    source("~/.Rprofile")
    qwt(menvar,file="Mencheckvar.txt")
    
}

manual_check <- function(){
        
    ## run in cluster 
    load("alllist_9_10")
    vars <- paste(alllist[,"Gene"],alllist[,1],alllist[,2],sep="_") 
    
    # PID  MID SUBID
    map <- matrix(c(220835,220842,250419,222806,222803,222805,231072,231073,231061,230378,230379,230189,241040,241024,241048),5,3,byrow=T)
    a <- c("OR5AC2_3_97806080","PTPRK_6_128403675","REG3G_2_79254997","BRD4_19_15350575","POLM_7_44119220","TRDN_6_123869623","BRCA2_13_32950808","HMCN1_1_186099172")
    b <- c("250419","222805","222805","231061","231061","230189","241048","241048")
    
    
    tmp <- c()
    for(i in 1:8){
        subs <- which(vars==a[i] & alllist[,"Subject_ID"] %in% map[map[,3]==b[i],])
        tmp <- rbind(tmp,alllist[subs,])
    }
    write.csv(tmp,file="Mencheck.csv",row.names=FALSE)
    
    
    
    
    
    #OR5AC2    3	97806080	rs141477221	C	T	exonic	stopgain	OR5AC2:NM_054106:exon1:c.C64T:p.R22X	0.0008	.	0.0008	.	T	.	.	D	2.2	26.5	none	High	High	T	.	PASS		250419
    
    ##in 220835 other disease
    
    #PTPRK	6	128403675	rs150342971	G	A	exonic	nonsynonymousSNV	PTPRK:NM_001291983:exon9:c.C1297T:p.H433Y	0.0003	0.000599042	0.0007	T	T	D	N	D	5.46	13.04	none	222805
    
    ##in 222806
    
    #REG3G	2	79254997	rs147143164	G	A	exonic	stopgain	REG3G:NM_001270040:exon4:c.G260A:p.W87X	0.0002	0.000199681	7.7e-05	.	T	.	.	A	4.64	14.53	none	High	High	A	.	222805	
    
    # in 222803
    
    #BRD4	19	15350575	.	TCTC	T	exonic	nonframeshiftdeletion	BRD4:NM_058243:exon16:c.3337_3339del:p.1113_1113del	0.0007	.	0.0005	.	.	.	.	.	.	.	none	High	High	231061
    
    # 231073 other
    
    #POLM	7	44119220	.	G	A	exonic	stopgain	POLM:NM_001284330:exon4:c.C592T:p.Q198X	.	.	.	.	D	.	.	A	5.49	26.6	none	High	High	A	.	PASS		231061	
    
    # 231072 other disease
    
    #TRDN	6	123869623	rs201021891	C	T	exonic	nonsynonymousSNV 230189	 
    # 230378 other disease
    
    #BRCA2	13	32950808	.	AAACAC	A	exonic	frameshiftdeletion	BRCA2:NM_000059:exon21:c.8635_8639del:p.N2879fs	.	.	.	.	.	.	.	.	.	.	none	High	High	A	.	PASS		241048	
    
    # 241040 other disease
    
    #HMCN1	1	186099172	.	G	A	exonic	nonsynonymousSNV	HMCN1:NM_031935:exon84:c.G12979A:p.V4327M	0.0001	.	0.0002	T	T	D	L	N	5.62	18.27	none	Med	High	A	241048
    # 241040
    
}

vcf_check_WRN_CBL <- function(){
    varT <- read.delim("gene_mutations_dp4.txt",sep="\t")
    
    #subs <- which(varT[,3]==30945376)
    subs <- which(varT[,3]==119149355)
    a <- matrix(,6,5)
    tmp <- varT[subs,c("Subject_IDs","GT","AD","GT.AD.DP.GQ.PL","alldp4")]
    for(i in 1:5){
        if(i==2 | i==3){
            a0 <- unlist(strsplit(tmp[1,i],";"))
            a1 <- sapply(1:length(a0),function(k) unlist(strsplit(a0[k],":"))[2])
            a[,i] <- a1
        }else if(i==4){
            a0 <- unlist(strsplit(tmp[1,i],";"))
            a1 <- sapply(1:length(a0),function(k) substr(a0[k],8,nchar(a0[k])))
            a[,i] <- a1
        }else if(i==5){
            a0 <- unlist(strsplit(tmp[1,i],";"))
            a1 <- sapply(1:length(a0),function(k)  paste(unlist(strsplit(a0[k],":"))[4],unlist(strsplit(a0[k],":"))[5],sep=":" )     )
            a[,i] <- a1        
        }else{
            a[,i] <- unlist(strsplit(tmp[1,i],","))
        }
    }
    
    write.csv(a,file="WRN_CBL.csv",row.names=FALSE)
}