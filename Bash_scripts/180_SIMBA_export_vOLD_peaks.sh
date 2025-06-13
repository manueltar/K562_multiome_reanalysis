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
 
### simba_export

type=$(echo "simba_export")
outfile_simba_export=$(echo "$Log_files""outfile_10_""$type"".log")
touch $outfile_simba_export
echo -n "" > $outfile_simba_export
name_simba_export=$(echo "$type""_job")


Rscript_simba_export=$(echo "$Rscripts_path""479_simba_export_vOLD_PEAKS.R")

merged_clusters_final_annotated=$(echo "$output_dir""merged_clusters_final_annotated.rds")


# Define your variables
mem=$(echo "8000") # Memory per CPU in MB
processors=$(echo "4")
nodes=$(echo "1")
total_memory=$(echo "scale=0; ($mem / 1) * $processors  * $nodes" | bc)

echo "Calculated total_memory: $total_memory"




myjobid_simba_export=$(sbatch --job-name $name_simba_export --output=$outfile_simba_export --partition=cpuq --time=24:00:00 --nodes=$nodes --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_simba_export --merged_clusters_final_annotated $merged_clusters_final_annotated --processors $processors --total_memory $total_memory --type $type --out $output_dir")
myjobid_seff_simba_export=$(sbatch --dependency=afterany:$myjobid_simba_export --open-mode=append --output=$outfile_simba_export --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_simba_export >> $outfile_simba_export")

conda deactivate
