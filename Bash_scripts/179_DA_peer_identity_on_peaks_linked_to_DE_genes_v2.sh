#!/bin/bash

set -e

eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
 
MASTER_ROUTE=$1
analysis=$2
SeuratObject=$3
DE_genes=$4

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")


mkdir -p $Log_files


Diff_array=$(echo 'Diff_K562')

a=($(echo "$Diff_array" | tr "," '\n'))

 
declare -a array_2_length

array_2_length=${#a[@]}

for (( i=0; i<${array_2_length}; i=i+1 ));
do

    Diff_array_sel=${a[$i]}
    echo "$Diff_array_sel"
    conda activate /home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis

    ### DA_function

    type=$(echo "DA_function""_""$Diff_array_sel""_""$analysis")
    outfile_DA_function=$(echo "$Log_files""outfile_1_""$type"".out")
    touch $outfile_DA_function
    echo -n "" > $outfile_DA_function
    name_DA_function=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    Rscript_DA_function=$(echo "$Rscripts_path""485_Multiome_DA_per_cluster_Diffs_K562_v2_only_linked_peaks_to_DE_genes_vOLD.R")




    # Define your variables
    mem=$(echo "8000") # Memory per CPU in MB
    processors=$(echo "20")
    nodes=$(echo "1")
    total_memory=$(echo "scale=0; ($mem / 1) * $processors * $nodes"  | bc)


    echo "Calculated total_memory: $total_memory"


    myjobid_DA_function=$(sbatch --job-name=$name_DA_function --output=$outfile_DA_function --partition=cpuq --time=24:00:00 --nodes=$nodes --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_DA_function --SeuratObject $SeuratObject --Diff_sel $Diff_array_sel --DE_genes $DE_genes --type $type --out $output_dir --processors $processors --total_memory $total_memory")
    myjobid_seff_DA_function=$(sbatch --dependency=afterany:$myjobid_DA_function --open-mode=append --output=$outfile_MSigDB_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_DA_function >> $outfile_DA_function")



    conda deactivate


done
