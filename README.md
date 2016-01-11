# breast_cancer_panel
The Panel gene lists and burden analysis for case-control data and families

Notes:
For case-control burden analysis, we remove 14 BRCA1/2 pathogenic or likely pathogenic subjects. In family-based SKAT analysis, we don't remove these subjects.
For sample 220673 is re-labeled as Jewish. The original phenotype information is labeled as Hispanic. Which are used in the cohort annotations for both Hispanic and Jewish.In the upudated, we should re-labeled all cohort samples.


For analysis before 10_20, we remove outliers based on three QC files (from Ashley): Potential_Outliers.tsv; CUMC_Regeneron.mismatches.sexcheck.tsv; Potential_Problem_families.tsv; For 10_20 updated, we incldue problem families. There are six subjects with corrected sequencing which are still not included by 10_20 version.

For family analysis before 10_20, we got 813 subjects with pedigree information. 39 trios data, and 285 subjects with at least one parents are sequenced. There are 287 breasr cancers and 526 no breast cancers.
