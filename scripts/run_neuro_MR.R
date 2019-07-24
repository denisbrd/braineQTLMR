setwd("")

library("TwoSampleMR")

#Sub function to run MR analysis on an exposure and outcome.
run_MR <- function(exposure_dat, outcome_id){
  
  #Get the names of the output files.
  
  harmon_output <- paste0("./harmon/", outcome_id , ".txt")
  default_MR_output <- paste0("./default_MR/", outcome_id , ".txt")
  het_MR_output <- paste0("./het_MR/", outcome_id , ".txt")
  pleio_MR_output <- paste0("./pleio_MR/", outcome_id , ".txt")
  single_SNP_output <- paste0("./singleSNP_MR/", outcome_id , ".txt")
  
  #--------------------------
  #Try to harmonise the files.
  
  #Get a unique list of SNPs to lookup in the outcome.
  snp_list <- unique(exposure_dat$SNP)
  
  print(paste0(as.character(length(snp_list))," SNPs to lookup in id", outcome_id))
  
  #Lookup the instrument SNPs in the outcome.
  
  outcome_dat <- NULL
  try(outcome_dat<-extract_outcome_data(snp_list,outcome_id))
  
  #Continue if instruments are available.
  if(!is.null(outcome_dat)){
    
    #--------------------------
    dat <- NULL
    #Harmonise the SNPs with the outcome, do not drop palindromic SNPs.
    try(dat <- harmonise_data(exposure_dat, outcome_dat, action=1))
    
    #If harmonised SNPs are present then go.
    if(!is.null(dat)){
      
      #Write out the harmonised result file.
      write.table(dat,file=harmon_output,sep="\t",col.names=T,row.names=F,quote=F)
      
      #--------------------------
      #Run the default MR tests. 
      
      mr_results <- NULL
      try(mr_results <- mr(dat))
      
      #Write table.
      if(!is.null(mr_results)){write.table(mr_results,file=default_MR_output,sep="\t",col.names=T,row.names=F,quote=F)}
      
      #--------------------------  
      #Run the single SNP MR (Wald Ratio estimates).
      
      mr_single <- NULL
      try(mr_single <- mr_singlesnp(dat))
      
      #Write table.
      if(!is.null(mr_single)){write.table(mr_single,file=single_SNP_output,sep="\t",col.names=T,row.names=F,quote=F)}
      
      #--------------------------
      #Run the het. analysis
      
      mr_hetero <- NULL
      try(mr_hetero <- mr_heterogeneity(dat))
      
      #Write table.
      if(!is.null(mr_hetero)){write.table(mr_hetero,file=het_MR_output,sep="\t",col.names=T,row.names=F,quote=F)}
      
      #--------------------------
      #Run the pleio analysis.
      
      mr_pleio <- NULL
      try(mr_pleio <- mr_pleiotropy_test(dat)) 
      
      #Write table.
      if(!is.null(mr_pleio)){write.table(mr_pleio,file=pleio_MR_output,sep="\t",col.names=T,row.names=F,quote=F)}
      
    }
    else{print(paste0("no snps extracted for id", outcome_id))}
  }
  else{print(paste0("no harmonised data for id", outcome_id))}
}

#Set up the directories to store the output.
dir.create(paste0("./harmon/"), showWarnings = F)
dir.create(paste0("./default_MR/"), showWarnings = F)
dir.create(paste0("./het_MR/"), showWarnings = F)
dir.create(paste0("./pleio_MR/"), showWarnings = F)
dir.create(paste0("./singleSNP_MR/"),showWarnings = F)

#Load in the instruments.

instrument_file <- ''

#Read in the instrume file.
exp_dat <- read_exposure_data(filename=instrument_file,
snp_col = "SNP",
beta_col = "standard.beta",
se_col = "standard.se",
effect_allele_col = "effect_allele.exposure",
other_allele_col = "other_allele.exposure",
eaf_col = "eaf.exposure",
pval_col = "pval.exposure",
phenotype_col = "exposure")

#Biogen (these are private)

#1175 (AD), 1180 (PD), 1181 (FTD), ALS march 2018

#Public available outcomes.  

#1183 (ADHD), 1186 (Anorexia Nervosa), 1185 (Autism Spectrum Disorder), 801 (Bipolar Disorder)
#990 (Bulimia nervorsa), 1188 (MDD), 1189 (Obessive compulsive disorder), 22 (Sz), 1024 (MS)

traits_to_run <- c('1186', '1185', '801', '990', '1188', '1189', '22', '1024')

#Get the MR-Base GWAS database.
ao <- available_outcomes()

#Subset for the outcomes of interest.
neuro.outcomes.df <- ao[ao$id%in%traits_to_run,]

#Check none of the traits are missing before running.
setdiff(traits_to_run, ao$id)

#Loop through each of the outcomes and run the MR analysis.
for(id in traits_to_run){
  
  print(paste0("Processing id", id, "..."))
  run_MR(exp_dat, id)
}



