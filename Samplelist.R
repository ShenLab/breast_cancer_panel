bamfiles1 <- list.files("/mnt/BigData/WENDY/RGN_rare_disease/bams","*.bam$",full.names = TRUE)
bamfiles2 <- list.files("/mnt/BigData/WENDY/BreastCancer/WES_data/Regeneron/bams","*.bam$",full.names = TRUE)
Allbams <- c(bamfiles2,bamfiles1)

AJs <- unlist(read.table("/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/AJ715_samples.list"))
HIs <- unlist(read.table("/home/local/ARCS/qh2159/breast_cancer/Panel/data/HispanicCases550.txt"))
HIs <- unique(c(HIs,"222357","222966","220897","220904"))
bams <- read.table("/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/Sample_batches.txt")
Sams <- matrix(0,1365,8)
Sams[,1] <- bams[,2]
Sams[Sams[,1] %in% AJs,2] <- "AJ"
Sams[Sams[,1] %in% HIs,2] <- "Dominican"
Sams[is.na(Sams[,2]),2] <- "Other"
Sams[,3] <- Allbams[match(bams[,1],basename(Allbams))]
#colnames(Sams) <- c("Sample-ID", "Inferred-Ethnicity","bam-file-location")

## add phenotype information
phenofile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/WES BCFR phenotypic data.csv" ## phenotype file
pheno <- read.csv(phenofile)
sub <- which(pheno[,3]=="222357, 222966")
a1 <- pheno[sub,]
a1[1,3] <- "222966"
pheno[sub,3] <- "222357"
pheno <- rbind(pheno,a1)

sub <- which(pheno[,3]=="220897, 220904")
a1 <- pheno[sub,]
a1[1,3] <- "220904"
pheno[sub,3] <- "220897"
pheno <- rbind(pheno,a1)

Sams[,4] <- pheno[match(Sams[,1],pheno[,3]),"FAMILYID"]
Sams[,5] <- pheno[match(Sams[,1],pheno[,3]),"BreastCancer"]
Sams[,6] <- pheno[match(Sams[,1],pheno[,3]),"CANCODE_1"]
Sams[,7] <- pheno[match(Sams[,1],pheno[,3]),"LiveAge"]
Sams[,8] <- pheno[match(Sams[,1],pheno[,3]),"Sex"]
Sams[is.na(Sams)] <- ""
colnames(Sams) <- c("Sample-ID", "Inferred-Ethnicity","bam-file-location","FAMILYID","BreastCancer","CANCODE_1","LiveAge","Sex")
qwt(Sams[,c(4,1,2,5,6,7,8,3)],file="../data/phenotype/BreastCancer_Samplelist.txt",flag=2)

aa <- cbind("BreastCancer",Sams[,c(1,3)])
colnames(aa) <- c("Study","SampleID","bam-file-location")
qwt(aa,file="../data/phenotype/Exome_Samplelist.txt",flag=2)
