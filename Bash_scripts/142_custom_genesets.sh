#!/bin/bash
 
MASTER_ROUTE=$1
analysis=$2


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0


bashrc_file=$(echo "/home/manuel.tardaguila/.bashrc")

source $bashrc_file
eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

type=$(echo "Custom_print_gmt""_""$analysis")
outfile_Custom_print_gmt=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_Custom_print_gmt
echo -n "" > $outfile_Custom_print_gmt
name_Custom_print_gmt=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_Custom_print_gmt=$(echo "$Rscripts_path""432_custom_gene_sets_create_gmt_v3.R")


Table_of_gene_sets=$(echo "/group/soranzo/manuel.tardaguila/gene_sets/Michelas_genesets_spreadsheet.txt")
miRNA_table=$(echo "/group/soranzo/manuel.tardaguila/gene_sets/miRPathDB.csv")
new_sets=$(echo "/group/soranzo/manuel.tardaguila/gene_sets/genesets.csv")
c3_collection=$(echo "/group/soranzo/manuel.tardaguila/gene_sets/c3.all.v2024.1.Hs.entrez_selected.rds")
Dorothea_ABC=$(echo "/group/soranzo/manuel.tardaguila/gene_sets/Dorothea_ABC_Hs.entrez_selected.rds")
Dorothea_ABCD=$(echo "/group/soranzo/manuel.tardaguila/gene_sets/Dorothea_ABCD_Hs.entrez_selected.rds")
SIMBA_GRN_targets=$(echo "/group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/simba/result_SIMBA/GRN_target_genes_CUX1_RUNX1.tsv")


myjobid_Custom_print_gmt=$(sbatch --job-name=$name_Custom_print_gmt --output=$outfile_Custom_print_gmt --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Custom_print_gmt --Table_of_gene_sets $Table_of_gene_sets --miRNA_table $miRNA_table --new_sets $new_sets --c3_collection $c3_collection --Dorothea_ABC $Dorothea_ABC --Dorothea_ABCD $Dorothea_ABCD --SIMBA_GRN_targets $SIMBA_GRN_targets --type $type --out $output_dir")
myjobid_seff_Custom_print_gmt=$(sbatch --dependency=afterany:$myjobid_Custom_print_gmt --open-mode=append --output=$outfile_Custom_print_gmt --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Custom_print_gmt >> $outfile_Custom_print_gmt")

