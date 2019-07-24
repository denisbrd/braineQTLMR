setwd("")
library("TwoSampleMR")

#Read in the clumped instruments.
output_dir <- "" 

exp_dat <- read_exposure_data(filename="./output/all.concat.run1.run2.run3.run4.txt",
                               snp_col = "SNP",
                               beta_col = "beta.exposure",
                               se_col = "se.exposure",
                               effect_allele_col = "effect_allele.exposure",
                               other_allele_col = "other_allele.exposure",
                               eaf_col = "eaf.exposure",
                               pval_col = "pval.exposure",
                               phenotype_col = "exposure")


#Get the MAF.
exp_dat$MAF <- ifelse(exp_dat$eaf.exposure>0.5, 1-exp_dat$eaf.exposure, exp_dat$eaf.exposure)

#Get the Z.

exp_dat$Z <- exp_dat$beta.exposure/exp_dat$se.exposure

#Get the phenotypic variance.

exp_dat$pheno.var <- 2*exp_dat$MAF*(1-exp_dat$MAF)*(1286 + exp_dat$Z^2)

#Get the standardised beta and SE.
exp_dat$standard.beta <- exp_dat$Z/sqrt(exp_dat$pheno.var)
exp_dat$standard.se <- 1/sqrt(exp_dat$pheno.var)

#Re-calculate the Z score from the standarised value.

exp_dat$z.pval <- exp_dat$standard.beta/exp_dat$standard.se

#Compare to the original Z scores to check transform is OK.

plot(exp_dat$Z, exp_dat$z.pval)
exp_dat$z.pval <- NULL

#Compare the original and standardised betas to check transform is OK.
plot(exp_dat$standard.beta, exp_dat$beta.exposure)

#Write out the table.

write.table(exp_dat, "", quote=F, row.names=F)

