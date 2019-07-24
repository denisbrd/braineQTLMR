#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(coloc)

GWAS.file <- args[1] #GWAS file for the 1Mb region which has merged the exposure and outcome SNPs.
result.file <- args[2] #Name of the file to save the output to.
trait1 <- args[3] #Name of the outcome trait.
trait2 <- args[4] #Name of the exposure trait.
region <- args[5] #Name of region (for example top hit SNP)
total.samplesize <- args[6] #Sample size for the outcome.
number.cases <- args[7] #Number of cases in the outcome.


#Approximate Bayes Factor colocalisation analyses

#H0: neither trait has a genetic association in the region
#H1:only trait 1 has a genetic association in the region.
#H2: only trait 2 has a genetic associaiton in the region.
#H3: both traits are associated, but with different causal variants
#H4: both traits are associated, and share a single causal variant

coloc.analysis.discrete <- function(input_gwas, beta.trait1,beta.trait2,se.trait1,se.trait2,MAF.trait1,MAF.trait2,N.trait1,N.trait1.cases, N.trait2, trait1, trait2, region){


  #type, quant (quantitative) for eQTL study.
  dataset1 <- list(beta=beta.trait2, varbeta=se.trait2^2, MAF=MAF.trait2,type="quant", N=N.trait2)
  #str(dataset1)

  #type, discrete for the outcome data.
  dataset2 <- list(beta=beta.trait1, varbeta=se.trait1^2,MAF=MAF.trait1, type="cc",N=N.trait1, s=N.trait1.cases/N.trait1)
  #str(dataset2)

  #Run the coloc analysis.
  result <- coloc.abf(dataset1, dataset2)

  #Format the data to save out.

  #List into data frame.
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  df$trait1 <- as.character(trait1)
  df$trait2 <- as.character(trait2)
  df$region <- as.character(region)

  #Label the columns in the data frame.
  names(df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf", "trait1", "trait2", "region")

  #Make the filename and save out.
  if (dim(df)[1] != 0){
    write.table(df,file=result.file,sep=" ",col.names=T,row.names=F,quote=F)
  }

}

#Load in the dataset.

GWAS.df <- read.table(GWAS.file, header=T, sep=' ', stringsAsFactors=FALSE)

#Set the sample size column.
GWAS.df$N.trait1<-as.numeric(total.samplesize)
GWAS.df$N.trait1.cases<-as.numeric(number.cases)

#Cast the columns.

#For the outcome (trait1).

GWAS.df$BETA.trait1 <- as.numeric(GWAS.df$BETA.trait1)
GWAS.df$SE.trait1 <- as.numeric(GWAS.df$SE.trait1)
GWAS.df$MAF.trait1 <- as.numeric(GWAS.df$MAF.trait1)
GWAS.df$N.trait1 <- as.numeric(GWAS.df$N.trait1)

#For the exposure (trait2).

GWAS.df$BETA.trait2 <- as.numeric(GWAS.df$BETA.trait2)
GWAS.df$SE.trait2 <- as.numeric(GWAS.df$SE.trait2)
GWAS.df$MAF.trait2 <- as.numeric(GWAS.df$MAF.trait2)
GWAS.df$N.trait2 <- as.numeric(GWAS.df$N.trait2)


#Fill in the missing MAF with the exposure MAF if missing in the outcome.

GWAS.df$MAF.trait1 <- ifelse(is.na(GWAS.df$MAF.trait1), GWAS.df$MAF.trait2, GWAS.df$MAF.trait1)

#Run the coloc analysis.

try(coloc.analysis.discrete(GWAS.file, GWAS.df$BETA.trait1, GWAS.df$BETA.trait2, GWAS.df$SE.trait1, GWAS.df$SE.trait2,GWAS.df$MAF.trait1, GWAS.df$MAF.trait2, GWAS.df$N.trait1, GWAS.df$N.trait1.cases, GWAS.df$N.trait2, trait1, trait2, region))
