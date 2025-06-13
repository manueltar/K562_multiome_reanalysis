# These are the scripts needed to run the alignment and QC of the multiome data of the rs139141690 K562 experiment

## 1. cellranger-arc-2.0.2 mapping

$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ MCO_1278
$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ MCO_1279
$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ MCO_1280
$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ MCO_1281

## 2. First QC steps: rna_min_features = 500, atac_min_fragments = 1000 and maximum percentage of mitochondrial reads = 10%

$ mkdir -p /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/processing_outputs/

$ bash ~/Scripts/Wraper_scripts/120_Seurat_first_v2.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ processing_outputs

## 3. Cellbender correction of ambient RNA

$ sbatch ~/Scripts/sbatch/7_CellBender.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ MCO_1278
$ sbatch ~/Scripts/sbatch/7_CellBender.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ MCO_1279
$ sbatch ~/Scripts/sbatch/7_CellBender.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ MCO_1280
$ sbatch ~/Scripts/sbatch/7_CellBender.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ MCO_1281

## 4. Merge ATAC peaks to create a unified fragments file

$ sbatch ~/Scripts/sbatch/8_merge_atac_peaks_vK562.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/

## 5. Derive the peak matrices of 5 kb bins

$ bash ~/Scripts/Wraper_scripts/122_snATAC_pipeline.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ processing_outputs

## 6. Align unmapped reads to a the reference of the GFP barcodes to genotype the cells

$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_vK562.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ /group/soranzo/paola.benaglio/references/modified_site/GFP_transgenev4.fa MCO_1278
$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_vK562.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ /group/soranzo/paola.benaglio/references/modified_site/GFP_transgenev4.fa MCO_1279
$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_vK562.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ /group/soranzo/paola.benaglio/references/modified_site/GFP_transgenev4.fa MCO_1280
$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_vK562.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ /group/soranzo/paola.benaglio/references/modified_site/GFP_transgenev4.fa MCO_1281

## 7. Filter and keep cells with a concordant barcode asignation from 3 different UMIs

$ bash ~/Scripts/Wraper_scripts/119_Filter_Larry_and_graphs_vK562.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/deconvolute_LARRY/ count_and_filter /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/deconvolute_LARRY/

## 8. Run AMULET to detect doublets in the ATAC modality

$ bash ~/Scripts/Wraper_scripts/123_Amulet.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ processing_outputs

## 9. Run the second pass of Seurat

$ bash ~/Scripts/Wraper_scripts/125_Seurat_second_pass_vK562.sh

## 10. Merge all samples

$ bash ~/Scripts/Wraper_scripts/126_merge_pre_merged_per_sample_vK562.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ processing_outputs /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/processing_outputs/merged.atac_fragments.tsv.gz

## 11. Final QC notebook

----> Jupyter notebook: Final_QC_in_the_merged_object.ipynb

## 12. Recluster and export h5ad for rpca

$ bash ~/Scripts/Wraper_scripts/153_Recluster_and_export_h5ad.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ processing_outputs

## 13. Add genotyping information

----> Jupyter notebook: Post_QC_genotype.ipynb (with Paola's annotation to match the paper)

## 14. Subset, genotype and call peaks

$ bash ~/Scripts/Wraper_scripts/177_Recluster_and_peak_calling.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ processing_outputs_Paola_genotype


====================> new path: /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/

## 15. Characterize the final object

----> Jupyter notebook: Post_genotype_characterization_Downstream_analysis.ipynb
----> Jupyter notebook: Figure_5_and_S5_panels_B_C_and_D_downstream_analysis.ipynb


## 16. SIMBA export

$ bash ~/Scripts/Wraper_scripts/180_SIMBA_export_vOLD_peaks.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/ Downstream_analysis

## 17. DE analysis in Pseudobulks

mkdir -p /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/DE_per_cluster/

$ bash ~/Scripts/Wraper_scripts/178_DE_per_identity_v2.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/ DE_per_cluster /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/merged_clusters_final_annotated.rds

## 18. DA analysis in Pseudobulks

mkdir -p /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/DA_per_cluster/

$ bash ~/Scripts/Wraper_scripts/179_DA_peer_identity_on_peaks_linked_to_DE_genes_v2.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/ DA_per_cluster /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/merged_clusters_final_annotated.rds /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/DE_per_cluster/DE_results_Diff_K562.rds



####################		PASSED TO GITHUB LINE #################### #################### ####################



18. SIMBA find kmers

$ bash ~/Scripts/Wraper_scripts/168_Simba_scan_for_kmers_motifs_v3.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/ simba /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/Peaks.bed

19. SIMBA preprocessing

$ bash ~/Scripts/Wraper_scripts/170_Python_SIMBA_preprocessing_vK562.sh /group/soranzo/manuel.tardaguila/2025_K562_multiome_reanalysis/Downstream_analysis/ simba


==========================================================================> test_DA, find a way to classify peaks TSS vs non TSS
==========================================================================> check

$ awk -F"\t" 'NR == 1 ; {if($1 == "CUX1" && $9 >= 1.3 ) print $0}' DE_per_cluster/DE_results_Diff_K562.tsv| awk -F"\t" 'NR==1;{if($8 == "1"||$8 == "3") print $0}'


$ awk -F"\t" 'NR == 1 ; {if($1 ~ /VOLUME/ && $12 >= 1.3 && $11 >=3) print $0}' ORA_results_significant_Diff_K562.tsv|awk -F"\t" 'NR==1;{if($16 == "1"||$16 == "3" || $16 == "2") print $0}'

$ awk -F"\t" 'NR==1;{if($1 ~ /VOLUME/ && $11 >= 1.3 && $13 != "time") print $0}' GSEA_results_significant_Diff_K562.tsv|awk -F"\t" 'NR==1;{if($14 == "1"||$14 == "3") print $0}'


