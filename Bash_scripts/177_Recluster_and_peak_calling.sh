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


Rscript_Recluster_and_peak_calling=$(echo "$Rscripts_path""478_genotyped_clustering_at_low_res_and_peak_calling.R")

db_filt_clustered_QCed_genotyped=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_QCed_reclustered_genotyped.rds")
res_param=$(echo '0.5')
frag_file=$(echo "$output_dir""merged.atac_fragments.tsv.gz")


# Define your variables
mem=$(echo "12000") # Memory per CPU in MB
processors=$(echo "30")

# Calculate total_memory using bc for floating-point arithmetic
# Scale=0 rounds to the nearest integer.
# If you want to keep the decimal, remove 'scale=0'.
# However, R's future.globals.maxSize expects an integer (bytes).
total_memory=$(echo "scale=0; ($mem / 1) * $processors" | bc)

echo "Calculated total_memory: $total_memory"



myjobid_Recluster_and_peak_calling=$(sbatch --job-name $name_Recluster_and_peak_calling --output=$outfile_Recluster_and_peak_calling --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_Recluster_and_peak_calling --db_filt_clustered_QCed_genotyped $db_filt_clustered_QCed_genotyped --res_param $res_param --frag_file $frag_file --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_Recluster_and_peak_calling=$(sbatch --dependency=afterany:$myjobid_Recluster_and_peak_calling --open-mode=append --output=$outfile_Recluster_and_peak_calling --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Recluster_and_peak_calling >> $outfile_Recluster_and_peak_calling")

### simba_export

type=$(echo "simba_export")
outfile_simba_export=$(echo "$Log_files""outfile_10_""$type"".log")
touch $outfile_simba_export
echo -n "" > $outfile_simba_export
name_simba_export=$(echo "$type""_job")


Rscript_simba_export=$(echo "$Rscripts_path""479_simba_export.R")

merged_clusters_final_annotated=$(echo "$output_dir""merged_clusters_final_annotated.rds")


# Define your variables
mem=$(echo "8000") # Memory per CPU in MB
processors=$(echo "4")

# Calculate total_memory using bc for floating-point arithmetic
# Scale=0 rounds to the nearest integer.
# If you want to keep the decimal, remove 'scale=0'.
# However, R's future.globals.maxSize expects an integer (bytes).
total_memory=$(echo "scale=0; ($mem / 1) * $processors" | bc)

echo "Calculated total_memory: $total_memory"


# --dependency=afterany:$myjobid_Recluster_and_peak_calling

myjobid_simba_export=$(sbatch --dependency=afterany:$myjobid_Recluster_and_peak_calling --job-name $name_simba_export --output=$outfile_simba_export --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_simba_export --merged_clusters_final_annotated $merged_clusters_final_annotated --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_simba_export=$(sbatch --dependency=afterany:$myjobid_simba_export --open-mode=append --output=$outfile_simba_export --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_simba_export >> $outfile_simba_export")

conda deactivate
