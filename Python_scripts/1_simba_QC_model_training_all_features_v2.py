#!/usr/bin/env python3
import os
import simba as si
si.settings.set_figure_params(dpi=80,
                              style='white',
                              fig_size=[5,5],
                              rc={'image.cmap': 'viridis'})

from IPython.display import set_matplotlib_formats
set_matplotlib_formats('retina')
import scanpy as sc
import sys
import argparse
import logging
import subprocess
import time
import numpy as np
import pandas as pd
import pickle

def QC_RNA(args):
	logging.info("RNA QC")
	adata_CG = sc.read(args.input_rna)

	si.pp.filter_genes(adata_CG, min_n_cells=3)
	si.pp.cal_qc_rna(adata_CG)
	si.pp.normalize(adata_CG, method='lib_size')
	si.pp.log_transform(adata_CG)
	
	logging.info("Discretize RNA values in 5 bins")
	si.tl.discretize(adata_CG, n_bins=5)
	si.pl.discretize(adata_CG, kde=False)
	
	workdir = args.output
	
	# 🚫 Clean up reserved '_index' column if present in var
	if '_index' in adata_CG.var.columns:
		adata_CG.var.rename(columns={'_index': 'var_index'}, inplace=True)
	if '_index' in adata_CG.obs.columns:
		adata_CG.obs.rename(columns={'_index': 'obs_index'}, inplace=True)
	# 📌 Set raw after cleanup
	adata_CG.raw = adata_CG  # make sure .raw.var gets the clean version
	# 🔁 Just to be extra safe: check again in raw.var
	if '_index' in adata_CG.raw.var.columns:
		adata_CG.raw.var.rename(columns={'_index': 'raw_var_index'}, inplace=True)
	logging.info("Saving processed RNA")
	adata_CG.write(os.path.join(workdir, 'adata_CG.h5ad'))
	return adata_CG

def QC_ATAC(args):
	logging.info("ATAC QC")
	adata_CP = sc.read(args.input_atac)

	si.pp.filter_genes(adata_CP, min_n_cells=3)
	si.pp.cal_qc_rna(adata_CP)  # assuming ATAC is stored like RNA here
	si.pp.normalize(adata_CP, method='lib_size')
	si.pp.log_transform(adata_CP)

	logging.info("Discretize ATAC values in 5 bins")
	si.tl.discretize(adata_CP, n_bins=5)
	si.pl.discretize(adata_CP, kde=False)

	workdir = args.output

	# 🔒 Clean reserved column names before setting raw
	if '_index' in adata_CP.var.columns:
		adata_CP.var.rename(columns={'_index': 'var_index'}, inplace=True)

	if '_index' in adata_CP.obs.columns:
		adata_CP.obs.rename(columns={'_index': 'obs_index'}, inplace=True)

	adata_CP.raw = adata_CP  # ensure raw gets the cleaned var

	if '_index' in adata_CP.raw.var.columns:
		adata_CP.raw.var.rename(columns={'_index': 'raw_var_index'}, inplace=True)

	logging.info("Saving processed ATAC")
	adata_CP.write(os.path.join(workdir, 'adata_CP.h5ad'))

	return adata_CP

def read_kmers_and_motifs(args):
	logging.info("Read kmers and motifs from Rscript")
	adata_PK = si.read_hdf(args.input_kmers,'mat')
	adata_PM = si.read_hdf(args.input_motifs,'mat')
	logging.info("Convert byte string to string")
	adata_PK.obs.index = [x.decode('utf-8') for x in adata_PK.obs.index]
	adata_PK.var.index = [x.decode('utf-8') for x in adata_PK.var.index]
	adata_PM.obs.index = [x.decode('utf-8') for x in adata_PM.obs.index]
	adata_PM.var.index = [x.decode('utf-8') for x in adata_PM.var.index]
	si.pp.binarize(adata_PK)
	si.pp.binarize(adata_PM)

	# Assuming adata_PM is already loaded

	# Store the original var_names in the .var DataFrame
	adata_PM.var['original_var_names'] = adata_PM.var_names

	original_var_names_list = adata_PM.var_names.tolist() # Get the current (original) names as a list
	transformed_var_names = []

	for var_name in original_var_names_list:
		transformed_name = var_name  # Initialize transformed_name
		if '/' in var_name:
			parts = var_name.split('/')
			first_part = parts[0]
			if '(' in first_part:
				transformed_name = first_part.split('(')[0]
			else:
				transformed_name = first_part
		else:
			if '(' in var_name:
				transformed_name = var_name.split('(')[0]
			else:
				transformed_name = var_name
		# Apply the second substitution (if relevant)
		if ':' in transformed_name and transformed_name.startswith('M_'):
			transformed_name = transformed_name.replace(':', '_')
		# Apply the third substitution (if relevant)
		if 'PU.1' in transformed_name:
			transformed_name = transformed_name.replace('PU.1', 'SPI1')
		transformed_var_names.append(transformed_name)
	# Assign the transformed list as the new var_names (index of .var)
	adata_PM.var_names = pd.Index(transformed_var_names)
	adata_PM.var_names_make_unique() # Add this line
	print("First 5 rows of the .var DataFrame:")
	print(adata_PM.var.head())
	print("\nFirst 10 of the transformed var_names:")
	print(adata_PM.var_names[:10])
	print("\nFirst 700-710 of the original var_names:")
	print(adata_PM.var['original_var_names'][700:710])
	workdir = args.output
	adata_PM.write(os.path.join(workdir,'adata_PM.h5ad'))
	adata_PK.write(os.path.join(workdir,'adata_PK.h5ad'))
	return adata_PK, adata_PM

def graph_generation_and_train_model (args, adata_CG, adata_CP, adata_PK, adata_PM):
	path_to_graph = args.graph_dir
	emb_path = args.path_emb
	processor_number = args.threads
	
	# modify parameters
	dict_config = si.settings.pbg_params.copy()
	dict_config['workers'] = processor_number	
	si.settings.pbg_params = dict_config # Update the settings

	adata_CG.uns['simba'] = {'name': 'C'}
	adata_CP.uns['simba'] = {'name': 'P'}
	adata_PM.uns['simba'] = {'name': 'M'}
	adata_PK.uns['simba'] = {'name': 'K'}
	logging.info(f"graph stage started at {time.strftime('%H:%M:%S')}")
	logging.info(f"adata_CG.obs.index[:5]: {list(adata_CG.obs.index[:5])}") # Print first 5 indices
	logging.info(f"adata_CG.var.index[:5]: {list(adata_CG.var.index[:5])}")
	logging.info(f"adata_CP.obs.index[:5]: {list(adata_CP.obs.index[:5])}")
	logging.info(f"adata_CP.var.index[:5]: {list(adata_CP.var.index[:5])}")
	logging.info(f"adata_PK.obs.index[:5]: {list(adata_PK.obs.index[:5])}")
	logging.info(f"adata_PK.var.index[:5]: {list(adata_PK.var.index[:5])}")
	logging.info(f"adata_PM.obs.index[:5]: {list(adata_PM.obs.index[:5])}")
	logging.info(f"adata_PM.var.index[:5]: {list(adata_PM.var.index[:5])}")
	si.tl.gen_graph(list_CP=[adata_CP],
				list_CG=[adata_CG],
				list_PK=[adata_PK],
				list_PM=[adata_PM],
				copy=False,
				use_highly_variable=False,
				use_top_pcs=False,
				dirname=path_to_graph)
	logging.info(f"graph stage finished at {time.strftime('%H:%M:%S')}")
	si.load_graph_stats(path=path_to_graph)
	logging.info(f"Current si.settings.pbg_params: {si.settings.pbg_params}")
	logging.info(f"Model training started at {time.strftime('%H:%M:%S')}")
	si.tl.pbg_train(auto_wd=True, save_wd=True, output = emb_path)
	logging.info(f"Model training finished at {time.strftime('%H:%M:%S')}")
	# Plot metrics
	si.pl.pbg_metrics(fig_ncol=1)
	return

def extract_embeddings(args):
	path_to_graph = args.graph_dir
	emb_path = args.path_emb
	logging.info("Set path to the directory for pbg")
	si.load_pbg_config(path=emb_path)
	workdir = args.output
	logging.info("Load pre-existing graph and model")
	si.load_graph_stats(path=path_to_graph)# load in graph ('graph0') info	
	si.load_pbg_config(path=emb_path)# load in model info for ('graph0')
	logging.info("Load dictionary")	
	dict_adata = si.read_embedding()	
	keys_list = list(dict_adata.keys())
	print(keys_list)
	logging.info("extracting embeddings")	
	adata_C = dict_adata['C']  # embeddings for cells
	adata_G = dict_adata['G']  # embeddings for genes
	adata_P = dict_adata['P']  # embeddings for peaks
	adata_K = dict_adata['K']  # embeddings for kmers
	adata_M = dict_adata['M']  # embeddings for motifs
	#adata_M.obs.index = 'M_'+adata_M.obs.index #### to distinguish TF motif names from gene names in this case
	adata_M.obs.index = ['M_' + str(idx) for idx in adata_M.obs.index]
	adata_G.write(os.path.join(workdir,'adata_G.h5ad'))
	adata_P.write(os.path.join(workdir,'adata_P.h5ad'))
	adata_K.write(os.path.join(workdir,'adata_K.h5ad'))
	adata_M.write(os.path.join(workdir,'adata_M.h5ad'))
	return adata_C, adata_G, adata_P, adata_K, adata_M

def annotate_and_compare_cell_types(args, adata_CG, adata_C):
	workdir = args.output
	logging.info("Add annotation to cell embedding") ## Add annotation of celltypes (optional)	
	adata_C.obs['seurat_clusters'] = adata_CG[adata_C.obs_names,:].obs['seurat_clusters'].copy()
	adata_C.obs['time_point'] = adata_CG[adata_C.obs_names,:].obs['time_point'].copy()
	adata_C.obs['celltype'] = adata_CG[adata_C.obs_names,:].obs['diff_groups'].copy()
	adata_C.obs['clone_line'] = adata_CG[adata_C.obs_names,:].obs['clone_line'].copy()
	adata_C.obs['Genotype'] = adata_CG[adata_C.obs_names,:].obs['Genotype'].copy()
	logging.info("UMAP of cell embedding") ## Add annotation of celltypes (optional)	
	si.tl.umap(adata_C,n_neighbors=15,n_components=2)
	adata_C.write(os.path.join(workdir,'adata_C.h5ad'))
	return adata_C

def gene_compare_entities(args, adata_C, adata_G):
	logging.info("compare gene and  cell entities") ## Add annotation of celltypes (optional)
	adata_cmp_CG = si.tl.compare_entities(adata_ref=adata_C, adata_query=adata_G)
	workdir = args.output
	adata_cmp_CG.write(os.path.join(workdir,'adata_cmp_CG.h5ad'))
	return adata_cmp_CG

def motif_compare_entities(args, adata_C, adata_M):
	logging.info("compare motif and  cell entities") ## Add annotation of celltypes (optional)
	adata_cmp_CM = si.tl.compare_entities(adata_ref=adata_C, adata_query=adata_M)
	workdir = args.output
	adata_cmp_CM.write(os.path.join(workdir,'adata_cmp_CM.h5ad'))
	return adata_cmp_CM

def peak_compare_entities(args, adata_C, adata_P):
	logging.info("compare peaks and  cell entities") ## Add annotation of celltypes (optional)
	adata_cmp_CP = si.tl.compare_entities(adata_ref=adata_C, adata_query=adata_P)
	workdir = args.output
	adata_cmp_CP.write(os.path.join(workdir,'adata_cmp_CP.h5ad'))
	return adata_cmp_CP

def global_embedding(args, adata_C, adata_G, adata_M, adata_K, adata_P):
	workdir = args.output
	adata_all = si.tl.embed(adata_ref=adata_C, list_adata_query=[adata_G, adata_M, adata_K, adata_P])
	adata_all.write(os.path.join(workdir,'adata_all.h5ad'))
	return adata_all

def main(args):
	logging.info('Start.')

	# Initialize critical variables

	adata_CG = adata_CP = adata_PK = adata_PM = None

	if not args.skip_RNA_QC:
		logging.info(f'Loading RNA data from {args.input_rna}')
		adata_CG = QC_RNA(args)
	if not args.skip_ATAC_QC:
		logging.info(f'Loading ATAC data from {args.input_atac}')
		adata_CP = QC_ATAC(args)
	if not args.skip_kmers_motifs:
		logging.info(f'Loading kmers data from {args.input_kmers}')
		logging.info(f'Loading motifs data from {args.input_motifs}')		
		adata_PK, adata_PM = read_kmers_and_motifs(args)
	if not args.skip_graph_generation:
		workdir = args.output
		path_to_graph = args.graph_dir
		adata_CG = si.read_h5ad(os.path.join(workdir,'adata_CG.h5ad'))
		adata_CP = si.read_h5ad(os.path.join(workdir,'adata_CP.h5ad'))
		adata_PK = si.read_h5ad(os.path.join(workdir,'adata_PK.h5ad'))		
		adata_PM = si.read_h5ad(os.path.join(workdir,'adata_PM.h5ad'))
		# Check if all required data are available for graph generation
		if adata_CG is None or adata_CP is None or adata_PK is None or adata_PM is None:
			logging.error("Missing required input data for graph generation.")
			sys.exit(1)
		graph_generation_and_train_model(args, adata_CG, adata_CP, adata_PK, adata_PM)

	if not args.skip_compare_entities:		
		adata_C, adata_G, adata_P, adata_K, adata_M = extract_embeddings(args)
		workdir = args.output
		path_to_graph = args.graph_dir
		adata_CG = si.read_h5ad(os.path.join(workdir,'adata_CG.h5ad'))
		# Check if RNA data (adata_CG) is available for annotation
		if adata_CG is None:
			logging.error("Annotation requires RNA data (adata_CG). Please provide RNA input file.")
			sys.exit(1)	
		logging.info('annotate cells with labels')
		adata_C = annotate_and_compare_cell_types(args, adata_CG, adata_C)
		adata_cmp_CG = gene_compare_entities(args, adata_C, adata_G)
		adata_cmp_CM = motif_compare_entities(args, adata_C, adata_M)
		adata_cmp_CP = peak_compare_entities(args, adata_C, adata_P)
	if not args.skip_global_embedding:
		workdir = args.output
		adata_C = si.read_h5ad(os.path.join(workdir,'adata_C.h5ad'))
		adata_G = si.read_h5ad(os.path.join(workdir,'adata_G.h5ad'))
		adata_M = si.read_h5ad(os.path.join(workdir,'adata_M.h5ad'))		
		adata_P = si.read_h5ad(os.path.join(workdir,'adata_P.h5ad'))
		adata_K = si.read_h5ad(os.path.join(workdir,'adata_K.h5ad'))
		adata_all = global_embedding(args, adata_C, adata_G, adata_M, adata_K, adata_P)
		
	logging.info('Finish.')
	return

def process_args():
	parser = argparse.ArgumentParser(description='Simba initial steps.')
	io_group = parser.add_argument_group('I/O arguments')
	io_group.add_argument('-rna', '--input-rna', required=True, type=str, metavar='FILE', help='Path to input RNA (.h5ad) file')
	io_group.add_argument('-atac', '--input-atac', required=True, type=str, metavar='FILE', help='Path to input ATAC (.h5ad) file')
	io_group.add_argument('-kmers', '--input-kmers', required=True, type=str, metavar='FILE', help='Path to kmer prediction file')
	io_group.add_argument('-motifs', '--input-motifs', required=True, type=str, metavar='FILE', help='Path to motif prediction file')	
	io_group.add_argument('-graph', '--preexisting-graph', required=False, type=str, help='graph preexisting dir')
	io_group.add_argument('-model', '--preexisting-model', required=False, type=str, help='model preexisting dir')	
	io_group.add_argument('-emb', '--path-emb', required=False, type=str, help='path to pbg model')
	io_group.add_argument('-dir_graph', '--graph-dir', required=False, type=str, help='path to store graph')
	io_group.add_argument('-o', '--output', required=True, type=str, help='Output directory to store processed files')
	



	model_group = parser.add_argument_group('Model arguments')
	model_group.add_argument('-t', '--threads', required=True, type=int, help='Number of threads to use for alignment')
		

	skip_group = parser.add_argument_group('Skip steps')
	skip_group.add_argument('--skip-RNA_QC', required=False, action='store_true', default=False, help='Skip RNA_QC step')
	skip_group.add_argument('--skip-ATAC_QC', required=False, action='store_true', default=False, help='Skip ATAC_QC step')
	skip_group.add_argument('--skip-kmers_motifs', required=False, action='store_true', default=False, help='Skip kmers and motifs step')
	skip_group.add_argument('--skip-graph-generation', required=False, action='store_true', default=False, help='Skip graph generation and training model steps')
	skip_group.add_argument('--skip-compare-entities', required=False, action='store_true', default=False, help='Skip compare entities step')
	skip_group.add_argument('--skip-global_embedding', required=False, action='store_true', default=False, help='Skip global embedding step')		
	return parser.parse_args()

if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.INFO)
	args = process_args()
	main(args)

