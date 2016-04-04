# Breast cancer study
Identification of Germline Breast Cancer Risk Genes based on case-control and family data. All useful information, please read BreastCancerSummary.txt file.

# File list:
BreastCancerSummary.txt:	All the summary related with Breast Cancer Project.

oldR/:	Before getting real control data, we tried to use FB-SKAT, SKAT to do single variant association test and also burden analysis based on index cases with pseudo-controls, and ranking gene based on RWR algorithm. 

BreastCancerAna.R:	Main burden test code for case-control data, including variant filtering, single variant burden test, gene level burden test, and gene set level burden test. Parameter "swi" is used to analysis Jewish or Dominican, other parameters see details in the comment lines.

HISPbatches.R:	Running Dominican analysis with different parameters.

InheritedModels.R:	Automatical generating AD and AR inherited models for families.

Samplelist.R:	All cohort samples information. 

batchesRunBRtest.R:	Running Jewish analysis with different parameters.

misc.R:	Useful functions for analysis.

pedigree.R:	Reading svg pedigree to trios table.

sourcefiles.R:	All data files used in BreastCancerAna.R to make this analysis directory independent.

src.R:	Useful functions for various phenotype and other information statistics and pre-running for singleVariantTest to speed up the BreastCancerAna.R

srcp.R:	 Other functions for more specific tasks, which is privately conserved by myself.

Readme.R:	Other related analysis source location and codes in this project. 


# Version Updates:

Version 0: burden test for single variant, single gene, gene sets. Add single variant test with popG: AJs: index cases, pseudo-controls, cases, non-cases; HIs: index cases, pseudo-controls, cases, non-cases; 

Population Sample Numbers: 271 58 65 313 (AJs); 142 41 23 343 (HIs);
Population Filtered Sample Numbers: 265 58 63 312 (AJs); 138 41 23 340 (HIs);

# Notes and remainning jobs:

For case-control burden analysis, we remove 14 BRCA1/2 pathogenic or likely pathogenic subjects. In family-based SKAT analysis, we don't remove these subjects.
For sample 220673 is re-labeled as Jewish. The original phenotype information is labeled as Hispanic. Which are used in the cohort annotations for both Hispanic and Jewish.In the upudated, we should re-labeled all cohort samples.

Fix Bad read group, Ashley origin code is under Regeneron/tablesandlists : xFixBadReadGroups\*
My code is under: depthOfcoverage/code xFixRGErrors_Piped.sh based on samtools.
This is only causal error when we use GATK to get depth of coverage (DOC). By now, we only analysis the 265 index cases in Jewish case-control study with DOC. We only fixed bad read group for these samples. We should fix for all samples in future. 
There are some read group error for 220365 bam file.

============= old version details==============

For analysis before 10_20, we remove outliers based on three QC files (from Ashley): Potential_Outliers.tsv; CUMC_Regeneron.mismatches.sexcheck.tsv; Potential_Problem_families.tsv; For 10_20 updated, we incldue problem families. There are six subjects with corrected sequencing which are still not included by 10_20 version.

For family analysis before 10_20, we got 813 subjects with pedigree information. 39 trios data, and 285 subjects with at least one parents are sequenced. There are 287 breasr cancers and 526 no breast cancers.
