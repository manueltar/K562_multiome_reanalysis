#!/bin/bash

eval "$(conda shell.bash hook)"
  

python_path=$(echo "/home/manuel.tardaguila/Scripts/Python/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

echo "$Log_files"

mkdir -p $Log_files

conda activate /group/soranzo/conda_envs/env_simba

### SIMBA_QC_and_embeddings

type=$(echo "SIMBA_QC_and_embeddings")
outfile_SIMBA_QC_and_embeddings=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_SIMBA_QC_and_embeddings
echo -n "" > $outfile_SIMBA_QC_and_embeddings
name_SIMBA_QC_and_embeddings=$(echo "$type""_job")


python_script_SIMBA_QC_and_embeddings=$(echo "$python_path""1_simba_QC_model_training_all_features_v2.py")

memory=$(echo "12000")
processors=$(echo "30")






### I had to edit the file /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/result_SIMBA/pbg/graph0/model/config.json to include the full path 

cd $MASTER_ROUTE

input_rna=$(echo "/group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/RNA.h5ad")
input_atac=$(echo "/group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/ATAC.h5ad")
input_kmers=$(echo "/group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/simba/result_SIMBA/freq_kmer.h5")
input_motifs=$(echo "/group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/simba/result_SIMBA/freq_motif.h5")
graph_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/simba/result_SIMBA/pbg/graph0/")
path_emb=$(echo "/group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/simba/result_SIMBA/pbg/graph0/model/")



path_entity_alias=$(echo "/group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/simba/result_SIMBA/pbg/graph0/")
path_entity=$(echo "/group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/simba/result_SIMBA/pbg/graph0/input/entity/")

final_dir=$(echo "$output_dir""result_SIMBA/")

myjobid_SIMBA_QC_and_embeddings=$(sbatch --job-name $name_SIMBA_QC_and_embeddings --output=$outfile_SIMBA_QC_and_embeddings --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$memory --parsable --wrap="python3 $python_script_SIMBA_QC_and_embeddings --threads $processors --input-rna $input_rna --input-atac $input_atac --input-kmers $input_kmers --input-motifs $input_motifs --graph-dir $graph_dir --path-emb $path_emb --output $final_dir")
myjobid_seff_SIMBA_QC_and_embeddings=$(sbatch --dependency=afterany:$myjobid_SIMBA_QC_and_embeddings --open-mode=append --output=$outfile_SIMBA_QC_and_embeddings --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_SIMBA_QC_and_embeddings >> $outfile_SIMBA_QC_and_embeddings")

conda deactivate

# --skip-RNA_QC --skip-ATAC_QC --skip-kmers_motifs --skip-graph-generation --preexisting-graph $preexisting_graph --preexisting-model $preexisting_model


# --skip-compare_entities --skip-global_embedding
