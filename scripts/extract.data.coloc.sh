#!/bin/bash

#Extracts the regions for coloc from the exposure and outcome GWASs and merges together the SNPs.

#System parameters
#$1 - file containing information on the region to be extracted (header is SNP (instrument), chrom, pos, start , stop, outcome, exposure)
#$2 - output directory to place the joined output

#Outcome GWASs in bcf format.
trait1_input_dir=''

#Exposure GWASs, brain eQTL gene expression files split by gene and sorted by rsid.
trait2_input_dir='/panfs/panasas01/sscm/epdab/brain_eQTL/split/'

> $1".lookup.log"

echo "regions to lookup: number=`wc -l $1`"
echo "output files directory: $2"

#Read through the list of regions to be extracted (as specified by $1).
while read line;
do
  echo "Processing "$line

  #Split out the line to get the fields required.
  arr=( $(IFS=" " echo "$line") )

  region=${arr[0]}
  chrom=${arr[1]}
  start=${arr[3]}
  end=${arr[4]}
  trait1=${arr[5]}
  trait2=${arr[6]}

  #echo $region $chrom $start $end $trait1 $trait2

  #Make a directory to store the output for the outcome.
  trait1_dir=$2$trait1"/"

  if [ ! -d "$trait1_dir" ]
  then
    mkdir $trait1_dir
  fi

  #--------------------------------------
  #Extract the data from the outcome based on the start and stop pos.

  #Make a unique id.
  id=$trait1_dir$region"."$trait2
  trait1_output=$id".outcome.lookup.txt"

  #Write in the header.
  echo "SNP chrom pos EA NONEA BETA SE PVAL EAF N MAF trait1" > $trait1_output

  #Lookup region in outcome GWAS, remove sites without rsids, extract out the columns needed, calculate the maf, sort by rsid.
  bcftools view -r "$chrom:$start-$end" $trait1_input_dir$trait1".tab.bcf" | bcftools query -f '%ID %CHROM %POS %ALT %REF %INFO/B %INFO/SE %INFO/PVAL %INFO/AF %INFO/N\n' | awk '{if($1~/^rs/) print}' | awk -v var=$trait1 '{if($9<=0.5) print $0, $9, var; else print $0, 1-$9, var}' | sort -k1,1 >> $trait1_output
  wc1=`wc -l $trait1_output | awk '{print $1}'`
  #--------------------------------------

  #Extract the SNPs from the exposure GWAS.

  #Name the file.
  trait2_output=$id".exposure.lookup.txt"

  #Write the header.
  echo "SNP EA NONEA EAF BETA SE PVAL N MAF trait2"> $trait2_output

  #Extract the sites with rsids from the exposure GWAS, sort by rsid, calculate the MAF.
  awk '{if($4~/^rs/) print $4, $5, $6, $7, $8, $9, $10, "1286"}' $trait2_input_dir$trait2".tab" | awk -v var=$trait2 '{if($4<=0.5) print $0, $4, var; else print $0, 1-$4, var}' | sort -k1,1 >> $trait2_output
  wc2=`wc -l $trait2_output| awk '{print $1}'`

  #--------------------------------------
  #Join the exposure and outcome GWASs.

  #Name the file.
  trait1and2_output=$id".joined.txt"

  #Write the header.
  echo "SNP chrom pos EA.trait1 NONEA.trait1 BETA.trait1 SE.trait1 PVAL.trait1 EAF.trait1 N.trait1 MAF.trait1 trait1 EA.trait2 NONEA.trait2 EAF.trait2 BETA.trait2 SE.trait2 PVAL.trait2 N.trait2 MAF.trait2 trait2" > $trait1and2_output

  #Join the dataset.
  join $trait1_output $trait2_output | grep -v ^"SNP" >> $trait1and2_output

  wc3=`wc -l $trait1and2_output| awk '{print $1}'`

  #Write out the SNP counts for each step to a log file.
  echo $line "$wc1" "$wc2" "$wc3" >> $1".lookup.log"

done < $1
