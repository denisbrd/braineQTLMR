# braineQTLMR

## Step 1 LD clump the instruments.

clump.instruments.MRBase.R

* clumps the top hit SNPs (P<5x10-8) for each gene.

## Step 2  Standardise the eQTL estimates.

standarise.betas.R

* Applies the standardisation formula to the eQTL SNP estimates from the clumped output.

## Step 3 Perform the MR analysis.

run_neuro_MR.R

* Runs the MR analysis between the eQTLs and neuro/pysch outcome GWASs.

##  Step 4 Run the coloc analysis.

extract.data.coloc.sh

* Extracts the region from the outcome and exposure GWAS to perform the coloc on (1Mb around the instrument.)  Harmonises the SNP effects to produce a single merged file which is used as the input for coloc function.

run.coloc.discrete.R

* Performs Bayesian coloc analysis on the harmonised dataset.

## DOI

### v1.0.1
