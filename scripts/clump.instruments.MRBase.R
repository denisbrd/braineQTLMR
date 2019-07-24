setwd("")

#library(devtools)
#install_github("MRCIEU/TwoSampleMR")

library("TwoSampleMR")

output_dir <- ''

#Sub function to clump the SNPs within each gene transcript.
clump_mydata <- function(exp_data, gene_list){
  
  idx <- 0
  for(gene in gene_list){
    
    idx <- idx + 1
    
    print(paste0("Processing exposure :", gene))
    
    print(paste0("Processing ", as.character(idx), " of ", length(gene_list)))
    
    #Get the data frame.
    
    temp.df <- exp_dat[exp_dat$exposure==gene,]
    print(head(temp.df))
    
    #Clump the data.
    
    clump.df<-NULL
    
    try(clump.df <- clump_data(temp.df))
    
    #Write the output if data was clumped.
    if(!is.null(clump.df)){
      
      output_file <- paste0(output_dir, gene, ".txt")
      write.table(clump.df, output_file, row.names=F, col.names = T, quote=F)
    }
    
  }
  
}


#Read in the exposure data (flat file of all the top hit SNPs from AMP-AD consortium).
exp_dat <- read_exposure_data(filename="AMP_Cortex_CMC_ROSMAP_Mayo_eQTL_metaanalysis_HRC_cleaned_rsid.noMHC.tophits.sorted.withheader.txt",
                                   sep = "\t",
                                   snp_col = "rsid",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   effect_allele_col = "EA",
                                   other_allele_col = "NONEA",
                                   eaf_col = "EAF",
                                   pval_col = "P",
                                   phenotype_col = "ENSG")

#Get the list of genes to clump.
gene_list <- unique(exp_dat$exposure)

#Clump the data.
clump_mydata(exp_dat, gene_list)

