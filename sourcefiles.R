outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"  ## outlier samples
BRCA1_2pathogenicfile <- "../data/phenotype/BRCA1_2.txt" ## BRCA1/2 pathogenic samples file
TSfile <- "../data/hotspots/Tumor_suppressors_11_11.txt" ## collected tumor suppressors
cancerdriverfile <- "../data/hotspots/Cancer_driver_11_6.txt" ## cancer drivers
DNArepairfile <- "../data/hotspots/DNA_repair_11_6.txt" ## DNA repair genes
Panelgfile <-  "../../genelist/Genelist2.txt" ## Panel genes
GTExfile <- "/home/local/ARCS/qh2159/breast_cancer/geneexpression/GTEx/expg_ranked.txt" ## GTEx expressed ranked gene
phenofile <- "../data/phenotype/WES BCFR phenotypic data.csv" ## phenotype file to get index cases only
dupIDs <- c("222357","222966") ## duplicated with Subject ID 222966 ## duplicated subject
## variant files and case, control samples ## check the log file to see the variant filtering details
AJcasevariantpath <- "/home/local/ARCS/qh2159/breast_cancer/variants/FreezeOneCommon/"
HIcasevariantpath <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"
AJcontvariantpath <- "/home/local/ARCS/qh2159/breast_cancer/variants/AJconVariantCalling/"
HIcontvariantpath <- "/home/local/ARCS/qh2159/breast_cancer/variants/HIconVariantCalling/"
AJBRfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/AJ715_samples.list"
HIBRfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/HispanicCases549.txt"
AJcontrolfile <- "/home/local/ARCS/qh2159/breast_cancer/variants/data/AJs_557.txt"
HIcontrolfile <- "/home/local/ARCS/qh2159/breast_cancer/variants/data/HIs_341.txt"
Cohortfile <- "../data/Rdata/BreastCancer_VariantList_11_12"
