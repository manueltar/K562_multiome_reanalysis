#!/bin/bash

set -e

eval "$(conda shell.bash hook)"
 
  
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2


##########################################################################

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate multiome_QC_DEF
 
### Recluster_and_peak_calling

type=$(echo "Recluster_and_peak_calling")
outfile_Recluster_and_peak_calling=$(echo "$Log_files""outfile_9_""$type"".log")
touch $outfile_Recluster_and_peak_calling
echo -n "" > $outfile_Recluster_and_peak_calling
name_Recluster_and_peak_calling=$(echo "$type""_job")


Rscript_Recluster_and_peak_calling=$(echo "$Rscripts_path""478_genotyped_clustering_at_low_res_and_peak_calling_vPEAK_old.R")

db_filt_clustered_QCed_genotyped=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_QCed_reclustered_genotyped.rds")
res_param=$(echo '0.5')
frag_file=$(echo "/group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/processing_outputs/""merged.atac_fragments.tsv.gz")


# Define your variables
mem=$(echo "12000") # Memory per CPU in MB
processors=$(echo "30")
nodes=$(echo "1")

# Calculate total_memory using bc for floating-point arithmetic
# Scale=0 rounds to the nearest integer.
# If you want to keep the decimal, remove 'scale=0'.
# However, R's future.globals.maxSize expects an integer (bytes).

total_processors=$(echo "scale=0; ($processors * $nodes)" | bc)
echo "Calculated total_processors: $total_processors"

total_memory=$(echo "scale=0; ($mem / 1) * $processors  * $nodes" | bc)


echo "Calculated total_memory: $total_memory"

myjobid_Recluster_and_peak_calling=$(sbatch --job-name $name_Recluster_and_peak_calling --output=$outfile_Recluster_and_peak_calling --partition=cpuq --time=24:00:00 --nodes=$nodes --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_Recluster_and_peak_calling --db_filt_clustered_QCed_genotyped $db_filt_clustered_QCed_genotyped --res_param $res_param --frag_file $frag_file --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_Recluster_and_peak_calling=$(sbatch --dependency=afterany:$myjobid_Recluster_and_peak_calling --open-mode=append --output=$outfile_Recluster_and_peak_calling --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Recluster_and_peak_calling >> $outfile_Recluster_and_peak_calling")


conda deactivate
