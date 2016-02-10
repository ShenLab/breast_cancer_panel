## Do not run this, this is only for detail information of related analysis in BreastCancer project

Here, we use "~" to present the main folder of breast_cancer.

variants/:	variants calling for our cohort data and controls. 
	pedigree:	All the pedigree information.
	code:	Main code for variant statistics and single sample variant calling and filtering.
	families:	AD and AR inherited models based variant calling for families.
	data:	useful data tables.
	AJconVariantCalling: AJ controls variant tables.
	HIconVariantCalling: HI controls variant tables.
	trios:	trios information.
	FreezeOneCommon: variant calling for AJ cases, we focused on the capture kits overlapped region for AJ cases and controls.
	HIdepth: variant depth from vcf files (HIs).
	AHdepth: variant depth from vcf files (AJs).	
	
genelist/:	collected gene list information, such as, tumor suppressor, cancer driver, DNA repair and others.

pedigree/:	old version pedigree figures.

depthOfcoverage/:	Depth of coverage for different batches samples both case and control based on bam files.
	code:	fix bad read group and limited number jobs in our PowerEdge cluster, and depth calling for different batches.
	data: 	data source and tables for bam files.
	Regeneron, RGN3, RGN4, SSC:	 the case-control used in AJ analysis, depth information for samples.
	AJcontrol and HIcontrol:	controls depth information for AJ and HI.


WES_Sinai/:	Regeneron WES and Sinai data comparsion.

Panel/:	Main analysis based on Panel gene sets.

geneexpression and somatic_mutation: Data source from TCGA.

