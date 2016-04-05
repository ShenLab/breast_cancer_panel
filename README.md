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

Version 0: 

1-Update Jewishes and Dominicans.  

		Total	IndexCases	Control

Jewish		715	265		557

Dominican	550	138		341

In details (filtered by outliers and BRCA1/2 positive samples):

	IndexCases	Pseudo-Controls	Controls	Non-indexCases	Non-Cases

AJs	271		58		557		65		313	

Filter	265		58		557		63		312

HIs	142		41		341		23		343

Filter	138		41		341		23		340


2-Splitting single variant test into src.R to save time

Add single variant test with different population sub-groups: 

AJs: index cases, pseudo-controls, cases, non-cases; HIs: index cases, pseudo-controls, cases, non-cases; 

# Notes:

1-For case-control burden analysis, we remove 14 BRCA1/2 pathogenic or likely pathogenic subjects. In family-based SKAT analysis, we don't remove these subjects.

2-For sample 220673 is re-labeled as Jewish. The original phenotype information is labeled as Hispanic. 

3-To fix bad read group, we use  depthOfcoverage/code xFixRGErrors_Piped.sh based on samtools. This is only causal error when we use GATK to get depth of coverage (DOC). By now, we only analysis the 265 index cases in Jewish case-control study with DOC. We only fixed bad read group for these samples. We should fix for all samples in future. 

4-Read group error for 220365 bam file which is different with 3 and not fixed yet by now.

============= old version details==============

For analysis before 10_20, we remove outliers based on three QC files (from Ashley): Potential_Outliers.tsv; CUMC_Regeneron.mismatches.sexcheck.tsv; Potential_Problem_families.tsv; For 10_20 updated, we incldue problem families. There are six subjects with corrected sequencing which are still not included by 10_20 version.

For family analysis before 10_20, we got 813 subjects with pedigree information. 39 trios data, and 285 subjects with at least one parents are sequenced. There are 287 breasr cancers and 526 no breast cancers.
