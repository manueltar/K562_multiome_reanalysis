#!/bin/bash

eval "$(conda shell.bash hook)"
 
  
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2


##########################################################################

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate multiome_QC_DEF
 
### Recluster at 0.2

type=$(echo "Recluster_and_export_h5ad")
outfile_Recluster_and_export_h5ad=$(echo "$Log_files""outfile_8_""$type"".log")
touch $outfile_Recluster_and_export_h5ad
echo -n "" > $outfile_Recluster_and_export_h5ad
name_Recluster_and_export_h5ad=$(echo "$type""_job")


Rscript_Recluster_and_export_h5ad=$(echo "$Rscripts_path""448_merged_clustering_at_low_res.R")

db_filt_clustered_QCed=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_QCed.rds")
res_param=$(echo '0.5')

myjobid_Recluster_and_export_h5ad=$(sbatch --job-name $name_Recluster_and_export_h5ad --output=$outfile_Recluster_and_export_h5ad --partition=cpuq --time=24:00:00 --nodes=4 --ntasks-per-node=30 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_Recluster_and_export_h5ad --db_filt_clustered_QCed $db_filt_clustered_QCed --res_param $res_param --type $type --out $output_dir")
myjobid_seff_Recluster_and_export_h5ad=$(sbatch --dependency=afterany:$myjobid_Recluster_and_export_h5ad --open-mode=append --output=$outfile_Recluster_and_export_h5ad --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Recluster_and_export_h5ad >> $outfile_Recluster_and_export_h5ad")

conda deactivate
