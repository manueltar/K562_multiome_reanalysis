#!/bin/bash
    
MASTER_ROUTE=$1
analysis=$2
indir=$3


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0


bashrc_file=$(echo "/home/manuel.tardaguila/.bashrc")

source $bashrc_file
eval "$(conda shell.bash hook)"
  

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files
    
#### Larry_counts_and_filter #############################


type=$(echo "Larry_counts_and_filter")
outfile_Larry_counts_and_filter=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_Larry_counts_and_filter
echo -n "" > $outfile_Larry_counts_and_filter
name_Larry_counts_and_filter=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")
 
Rscript_Larry_counts_and_filter=$(echo "$Rscripts_path""402_Larry_counts_CHEK2.R")


samples=$(echo "MCO_1278,MCO_1279,MCO_1280,MCO_1281")


Threshold_attributed_genotypes=$(echo "1")
Threshold_UMIS_per_cell=$(echo "3")
genotypes_string=$(echo "chrGFP_WTA,chrGFP_WTB,chrGFP_WTC,chrGFP_Del_287,chrGFP_Del_235,chrGFP_Del_233,chrGFP_Del_16bp,chrGFP_KI_29,chrGFP_KI_27,chrGFP_KI_13,chrGFP_HET")
mem=$(echo "4000")
processors=$(echo "4")
total_memory=$(( mem * processors ))

 
myjobid_Larry_counts_and_filter=$(sbatch --job-name=$name_Larry_counts_and_filter --output=$outfile_Larry_counts_and_filter --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_Larry_counts_and_filter --indir $indir --samples $samples --Threshold_attributed_genotypes $Threshold_attributed_genotypes --Threshold_UMIS_per_cell $Threshold_UMIS_per_cell --genotypes_string $genotypes_string --type $type --out $output_dir --processors $processors --total_memory $total_memory")
myjobid_seff_Larry_counts_and_filter=$(sbatch --dependency=afterany:$myjobid_Larry_counts_and_filter --open-mode=append --output=$outfile_Larry_counts_and_filter --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Larry_counts_and_filter >> $outfile_Larry_counts_and_filter")



