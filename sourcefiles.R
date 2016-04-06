### related data files using in our breast cancer cohort analysis
### === common files ==============
outlierfile <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/FreezeOneVariantCalling/Outliers_to_be_excluded.list"  ## outlier samples
BRCA1_2pathogenicfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/BRCA1_2.txt" ## BRCA1/2 pathogenic samples file
phenofile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/WES BCFR phenotypic data.csv" ## phenotype file
dupIDs <- c("222357","222966") ## duplicated with Subject ID 222966 ## duplicated subject

### === BreastCancerAna.R ==============
## gene sets
TSfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/hotspots/Tumor_suppressors_11_11.txt" ## collected tumor suppressors
cancerdriverfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/hotspots/Cancer_driver_11_6.txt" ## cancer drivers
DNArepairfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/hotspots/DNA_repair_11_6.txt" ## DNA repair genes
Panelgfile <-  "/home/local/ARCS/qh2159/breast_cancer/genelist/Genelist2.txt" ## Panel genes
GTExfile <- "/home/local/ARCS/qh2159/breast_cancer/geneexpression/GTEx/expg_ranked.txt" ## GTEx expressed ranked gene
pLIfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"
## interested pathways
MAPKfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/pathways/MAPK_pathways.txt"
PI3K_AKtfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/pathways/PI3K-Akt_pathways.txt"
Metabolismfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/pathways/ReactomeMetabolism_pathways.txt"
mapkgWendyfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/pathways/MAPKgenes_WENDY.txt"
pi3k_aktgWendyfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/pathways/PI3K_AKTgenes_WENDY.txt"
metabolismgWendyfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/pathways/Metabolicgenes_WENDY.txt"
pathwayfiles <- c(MAPKfile,PI3K_AKtfile,Metabolismfile,mapkgWendyfile,pi3k_aktgWendyfile,metabolismgWendyfile)
pathwaynames <- c("MAPK","PI3K_AKt","Metabolism","mapkgWendy","pi3k_aktgWendy","metabolismgWendy")
## hotspot files
hotHMM <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/hotspots/HMM_hotspots_11_12.txt"
hotCOSMIC <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/hotspots/cosmic_hotspots_3.txt"
## variant files and case, control samples ## check the log file to see the variant filtering details
AJcasevariantpath <- "/home/local/ARCS/qh2159/breast_cancer/variants/FreezeOneCommon/"
HIcasevariantpath <- "/home/local/ARCS/yshen/data/WENDY/BreastCancer/Regeneron/Filtering_Oct2015/"
AJcontvariantpath <- "/home/local/ARCS/qh2159/breast_cancer/variants/AJconVariantCalling/"
HIcontvariantpath <- "/home/local/ARCS/qh2159/breast_cancer/variants/HIconVariantCalling/"
AJBRfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/AJ715_samples.list"
HIBRfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/HispanicCases550.txt"
AJcontrolfile <- "/home/local/ARCS/qh2159/breast_cancer/variants/data/AJs_557.txt"
HIcontrolfile <- "/home/local/ARCS/qh2159/breast_cancer/variants/data/HIs_341.txt"
Cohortfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/Rdata/BreastCancer_VariantList_11_12"
SubDominif <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/Dominicans549_mix_label.txt" ## labeled Dominicans with sub-groups
Samplelistfile <- "/home/local/ARCS/qh2159/breast_cancer/Panel/data/phenotype/BreastCancer_Samplelist.txt"
## variant lists and statistical files
vburdenfiles <- c("/home/local/ARCS/qh2159/breast_cancer/Panel/resultf/AJ_variant_level_burden_Pseducont.txt","/home/local/ARCS/qh2159/breast_cancer/Panel/resultf/HISP_variant_level_burden_Pseducont.txt")
alleleFrefiles <- c("/home/local/ARCS/yshen/data/WENDY/BreastCancer/AJ_CONTROLS/combined_variant_call/NonBC_Frequencies.expanded.tsv","") ##AJ: not in our cohort 
caselistfs <- c("/home/local/ARCS/qh2159/breast_cancer/Panel/data/Rdata/AJcaselist715_12_17","/home/local/ARCS/qh2159/breast_cancer/Panel/data/Rdata/HIcaselist_1_29")
contlistfs <- c("/home/local/ARCS/qh2159/breast_cancer/Panel/data/Rdata/AJcontlist_12_17","/home/local/ARCS/qh2159/breast_cancer/Panel/data/Rdata/HIcontlist_1_5")
casestafs <- c("/home/local/ARCS/qh2159/breast_cancer/Panel/data/contAJ_20151209.hardfiltered.stats_hq.tsv","/home/local/ARCS/qh2159/breast_cancer/Panel/data/BR.origin.stats_12_28.tsv")
contstafs <- c("/home/local/ARCS/qh2159/breast_cancer/Panel/data/contAJ_20151209.hardfiltered.stats_hq.tsv","/home/local/ARCS/qh2159/breast_cancer/Panel/data/Hispanic.stats_12_28.tsv")


