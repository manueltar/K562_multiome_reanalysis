#!/bin/bash

set -e

eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
 
MASTER_ROUTE=$1
analysis=$2
SeuratObject=$3


TF_terms=$(echo 'CUX1,GATA1,GATA2,TET2,RUNX1,RUNX2,FOXM1,SPI1,PU1,FOXP1,GFI1')
search_terms=$(echo "PLATELET,ERYTHROCYTE,CUX1,MEGAKARYOCYTE,GATA1,GATA2,TET2,RUNX1,RUNX2,MITOSIS,ANEUPLOIDY,CYTOKINESIS,MYELOID,AML,HSC,HEMATOPOIETIC,SPHINGOSINE,FOXM1,SPI1,PU1,WP_PI3K_AKT_SIGNALING_PATHWAY,PI3K,AKT,FOXP1,GFI1")   # ADD HSC TERMS
path_to_GMT=$(echo "/home/manuel.tardaguila/GMT_files/msigdb_v2023.1.Hs_files_to_download_locally_ENTREZ/")

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

    ### DE_function

    type=$(echo "DE_function""_""$Diff_array_sel""_""$analysis")
    outfile_DE_function=$(echo "$Log_files""outfile_1_""$type"".out")
    touch $outfile_DE_function
    echo -n "" > $outfile_DE_function
    name_DE_function=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    Rscript_DE_function=$(echo "$Rscripts_path""480_Multiome_DE_per_cluster_Diffs_K562_v2.R")


    # Define your variables
    mem=$(echo "8000") # Memory per CPU in MB
    processors=$(echo "4")
    nodes=$(echo "1")
    total_memory=$(echo "scale=0; ($mem / 1) * $processors * $nodes"  | bc)


    echo "Calculated total_memory: $total_memory"


    myjobid_DE_function=$(sbatch --job-name=$name_DE_function --output=$outfile_DE_function --partition=cpuq --time=24:00:00 --nodes=$nodes --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_DE_function --SeuratObject $SeuratObject --Diff_sel $Diff_array_sel --type $type --out $output_dir --processors $processors --total_memory $total_memory")
    myjobid_seff_DE_function=$(sbatch --dependency=afterany:$myjobid_DE_function --open-mode=append --output=$outfile_MSigDB_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_DE_function >> $outfile_DE_function")





    ### volcano_function

    type=$(echo "volcano_function""_""$Diff_array_sel""_""$analysis")
    outfile_volcano_function=$(echo "$Log_files""outfile_2_""$type"".out")
    touch $outfile_volcano_function
    echo -n "" > $outfile_volcano_function
    name_volcano_function=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    Rscript_volcano_function=$(echo "$Rscripts_path""481_Multiome_DE_volcano_plots_both_Diffs.R")

    DE_results=$(echo "$output_dir""DE_results_""$Diff_array_sel"".rds")
    mem=$(echo "4000") # Memory per CPU in MB
    processors=$(echo "2")
    nodes=$(echo "1")
    total_memory=$(echo "scale=0; ($mem / 1) * $processors * $nodes"  | bc)

    echo "Calculated total_memory: $total_memory"


    # --dependency=afterany:$myjobid_DE_function

    myjobid_volcano_function=$(sbatch --dependency=afterany:$myjobid_DE_function --job-name=$name_volcano_function --output=$outfile_volcano_function --partition=cpuq --time=24:00:00 --nodes=$nodes --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_volcano_function --DE_results $DE_results --Diff_sel $Diff_array_sel --type $type --out $output_dir --processors $processors --total_memory $total_memory")
    myjobid_seff_volcano_function=$(sbatch --dependency=afterany:$myjobid_volcano_function --open-mode=append --output=$outfile_MSigDB_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_volcano_function >> $outfile_MSigDB_ORA")

    ### heatmap_function

    type=$(echo "heatmap_function""_""$Diff_array_sel""_""$analysis")
    outfile_heatmap_function=$(echo "$Log_files""outfile_3_""$type"".out")
    touch $outfile_heatmap_function
    echo -n "" > $outfile_heatmap_function
    name_heatmap_function=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    Rscript_heatmap_function=$(echo "$Rscripts_path""482_Multiome_DE_heatmap_both_Difss.R")

    DE_results=$(echo "$output_dir""DE_results_""$Diff_array_sel"".rds")
    normalised_counts=$(echo "$output_dir""norcounts_FINAL_""$Diff_array_sel"".rds")

    mem=$(echo "4000") # Memory per CPU in MB
    processors=$(echo "2")
    nodes=$(echo "1")
    total_memory=$(echo "scale=0; ($mem / 1) * $processors * $nodes"  | bc)

    echo "Calculated total_memory: $total_memory"

    # --dependency=afterany:$myjobid_DE_function

    myjobid_heatmap_function=$(sbatch --dependency=afterany:$myjobid_DE_function  --job-name=$name_heatmap_function --output=$outfile_heatmap_function --partition=cpuq --time=24:00:00 --nodes=$nodes --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_heatmap_function --DE_results $DE_results --normalised_counts $normalised_counts --Diff_sel $Diff_array_sel --type $type --out $output_dir --processors $processors --total_memory $total_memory")
    myjobid_seff_heatmap_function=$(sbatch --dependency=afterany:$myjobid_heatmap_function --open-mode=append --output=$outfile_MSigDB_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_heatmap_function >> $outfile_heatmap_function")



    conda deactivate



    ############################################################## ################################################################
    ############################################################## ################################################################
    ############################################################## ################################################################
    ############################################################## ################################################################

    conda activate GSEA

    ## MSigDB_GSEA

    type=$(echo "MSigDB_GSEA""_""$Diff_array_sel""_""$analysis")
    outfile_MSigDB_GSEA=$(echo "$Log_files""outfile_4_""$type"".out")
    touch $outfile_MSigDB_GSEA
    echo -n "" > $outfile_MSigDB_GSEA
    name_MSigDB_GSEA=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    Rscript_MSigDB_GSEA=$(echo "$Rscripts_path""483_GSEA_on_identity_both_Diffs.R")

    pval_threshold=$(echo "0.05")
    log2FC_threshold=$(echo "0")
    Threshold_number_of_genes=$(echo '3')

    DE_results=$(echo "$output_dir""DE_results_""$Diff_array_sel"".rds")  
    List_GSEA=$(echo "$output_dir""GSEA_complete_results_""$Diff_array_sel"".rds")
    GSEA_result=$(echo "$output_dir""GSEA_results_significant_""$Diff_array_sel"".rds")
    mem=$(echo "4000") # Memory per CPU in MB
    processors=$(echo "16")
    nodes=$(echo "1")
    total_memory=$(echo "scale=0; ($mem / 1) * $processors * $nodes" | bc)

    echo "Calculated total_memory: $total_memory"
    # --dependency=afterany:$myjobid_DE_function

    myjobid_MSigDB_GSEA=$(sbatch --dependency=afterany:$myjobid_DE_function --job-name=$name_MSigDB_GSEA --output=$outfile_MSigDB_GSEA --partition=cpuq --time=24:00:00 --nodes=$nodes --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_MSigDB_GSEA --path_to_GMT $path_to_GMT --search_terms $search_terms --pval_threshold $pval_threshold --log2FC_threshold $log2FC_threshold --Threshold_number_of_genes $Threshold_number_of_genes --TF_terms $TF_terms --DE_results $DE_results --Diff_sel $Diff_array_sel --List_GSEA $List_GSEA --GSEA_result $GSEA_result --type $type --out $output_dir --processors $processors --total_memory $total_memory")
    myjobid_seff_MSigDB_GSEA=$(sbatch --dependency=afterany:$myjobid_MSigDB_GSEA --open-mode=append --output=$outfile_MSigDB_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MSigDB_GSEA >> $outfile_MSigDB_GSEA")

    ### MSigDB_ORA

    type=$(echo "MSigDB_ORA""_""$Diff_array_sel""_""$analysis")
    outfile_MSigDB_ORA=$(echo "$Log_files""outfile_5_""$type"".out")
    touch $outfile_MSigDB_ORA
    echo -n "" > $outfile_MSigDB_ORA
    name_MSigDB_ORA=$(echo "$type""_job")
    seff_name=$(echo "seff""_""$type")

    Rscript_MSigDB_ORA=$(echo "$Rscripts_path""484_ORA_on_identity_both_Diffs.R")



    pval_threshold=$(echo "0.05")
    log2FC_threshold=$(echo "0")
    Threshold_number_of_genes=$(echo '3')
    DE_results=$(echo "$output_dir""DE_results_""$Diff_array_sel"".rds")  
    ORA_result=$(echo "$output_dir""ORA_results_significant_""$Diff_array_sel"".rds")

    mem=$(echo "4000") # Memory per CPU in MB
    processors=$(echo "16")
    nodes=$(echo "1")
    total_memory=$(echo "scale=0; ($mem / 1) * $processors * $nodes" | bc)

    # --dependency=afterany:$myjobid_DE_function

    echo "Calculated total_memory: $total_memory"

    myjobid_MSigDB_ORA=$(sbatch --dependency=afterany:$myjobid_DE_function  --job-name=$name_MSigDB_ORA --output=$outfile_MSigDB_ORA --partition=cpuq --time=24:00:00 --nodes=$nodes --ntasks-per-node=$processors --mem-per-cpu=$mem --parsable --wrap="Rscript $Rscript_MSigDB_ORA --DE_results $DE_results --path_to_GMT $path_to_GMT --search_terms $search_terms --pval_threshold $pval_threshold --log2FC_threshold $log2FC_threshold --Threshold_number_of_genes $Threshold_number_of_genes --TF_terms $TF_terms --ORA_result $ORA_result --Diff_sel $Diff_array_sel --type  $type --out $output_dir --processors $processors --total_memory $total_memory")
    myjobid_seff_MSigDB_ORA=$(sbatch --dependency=afterany:$myjobid_MSigDB_ORA --open-mode=append --output=$outfile_MSigDB_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MSigDB_ORA >> $outfile_MSigDB_ORA")


    conda deactivate

done
