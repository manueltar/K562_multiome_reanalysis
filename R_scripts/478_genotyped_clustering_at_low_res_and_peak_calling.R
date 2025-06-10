
.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/lib/R/library"))
.libPaths()
# sessionInfo()

Sys.setenv(RETICULATE_PYTHON="/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/bin/python")
library(reticulate)
reticulate::use_python("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/bin/python")
reticulate::use_condaenv("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF")
reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg')
suppressMessages(library("optparse"))
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(scDblFinder))
suppressMessages(library("tidyr"))
suppressMessages(library("tibble"))
suppressMessages(library("biovizBase"))
suppressMessages(library("patchwork"))
suppressMessages(library(glmGamPoi))
suppressMessages(library(SeuratData))
suppressMessages(library(SeuratDisk))
suppressMessages(library(future))
suppressMessages(library(plyr))



opt = NULL

options(warn = 1)

log_info_simple <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(timestamp, "INFO:", message, "\n")
}

cluster_macs2_simba = function(option_list)
{
  
 
  
  cat("All options:\n")
  printList(opt)
  
  
 
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform res_param ----
  
  res_param = opt$res_param
  
  cat("res_param_\n")
  cat(sprintf(as.character(res_param)))
  cat("\n")
  
  #### READ and transform processors ----
  
  processors = opt$processors
  
  cat("processors\n")
  cat(sprintf(as.character(processors)))
  cat("\n")
  
  #### READ and transform memory ----
  
  #### READ and transform total_memory (memory in MB) ----
  total_memory = opt$total_memory # Corrected variable name to match bash script
  cat("Total Memory (MB) for global objects:", as.character(total_memory), "\n") # Improved log message
  
  #### Assign resources -------------
  
  log_info_simple("plan stage")
  
  # Set up parallel processing: 'multiprocess' works on a single machine across cores.
  # 'total_memory' is expected in MB from the bash script, convert to bytes for future.globals.maxSize.
  plan("multicore", workers = processors)
  options(future.globals.maxSize = total_memory * 1024^2) # Corrected: Convert MB to bytes
  
  #### Read filtered object by doublets -----
  
  
  adata2<-readRDS(file=opt$db_filt_clustered_QCed_genotyped)
  
  log_info_simple("Subsetting to genotyped cells")
  
  adata_sub<-subset(adata2, Assignation_GEX_not_amplified != "NA")
  
  # cat("adata_sub_0\n")
  # cat(str(adata_sub))
  # cat("\n")
  
 
  metadata<-adata_sub@meta.data
  
  cat("metadata_0\n")
  cat(str(metadata))
  cat("\n")
  
  
  #### Cluster without integration and check if any cluster looks specially bad for QC metrics -------------
  
  log_info_simple("RNA clusters")
  
  DefaultAssay(adata_sub) <- 'RNA'
  
 
  
  adata_sub <- SCTransform(adata_sub, verbose = FALSE) 
  adata_sub <- RunPCA(adata_sub) 
  adata_sub <- RunUMAP(adata_sub, dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_')
  
  #### ATAC modality -------------
  # We exclude the first dimension as this is typically correlated with sequencing depth
  
  log_info_simple("ATAC clusters")
  
  
  DefaultAssay(adata_sub) <- 'ATAC'

  adata_sub <- RunTFIDF(adata_sub)
  adata_sub <- FindTopFeatures(adata_sub, min.cutoff='q0')
  adata_sub <- RunSVD(adata_sub)
  adata_sub <- RunUMAP(adata_sub, reduction='lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')
  
  
  #### WNN ATAC+RNA modality -------------
  
  log_info_simple("WNN")
  
  
  adata_sub <- FindMultiModalNeighbors(adata_sub, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:50))
  adata_sub <- RunUMAP(adata_sub, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
  adata_sub <- FindClusters(adata_sub, graph.name='wsnn', algorithm=4, resolution = res_param, verbose=FALSE, method = "igraph")
  
  
  
  
  
  ## Call Peaks and make new peak matrix----
  
  log_info_simple("call peaks based on seurat clusters")
  
  DefaultAssay(adata_sub) <- 'ATAC'
  peaks <- CallPeaks(
    object = adata_sub,
    group.by = "seurat_clusters",    
    macs2.path = "/group/soranzo/conda_envs/Manuel_macs2/bin/macs2")
  
  frag_file<-opt$frag_file
  
  Fragmobj <- CreateFragmentObject(frag_file,cells =Cells(adata_sub))
  
  
  
  log_info_simple("peakmat")
  
  peakmat = FeatureMatrix(fragments = Fragmobj, features = peaks, cells = Cells(adata_sub), process_n = 10000,
                          sep = c(":", "-"), verbose = TRUE)
  
  norm_chr = rownames(peakmat)[stringr::str_split_fixed(rownames(peakmat), ":",2)[,1] %in% 
                                 paste0("chr", c(1:22, "X", "Y"))]
  
  peakmat=peakmat[norm_chr,]
  
  suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
  seqlevelsStyle(annotations)  <- 'UCSC'
  genome(annotations)          <- 'hg38'
  
  
  
  
  suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=peakmat, sep=c(':', '-'), 
                                                       genome='hg38', fragments=Fragmobj, 
                                                       min.cells=-1, min.features=-1, 
                                                       annotation=annotations))
  
  
  
  adata_sub[['ATAC_by_seurat_clusters']] <- chrom_assay
  
  #### the new factors ----
  
  log_info_simple("new factors")
  
  adata_sub@meta.data$clone_line<-revalue(adata_sub@meta.data$Assignation_GEX_not_amplified,
                                      c('chrGFP_WTA' = 'wt_1',
                                        'chrGFP_WTB' = 'wt_2',
                                        'chrGFP_WTC' = 'wt_3',
                                        'chrGFP_KI_13' = 'rs139141690_1',
                                        'chrGFP_KI_27' = 'rs139141690_2',
                                        'chrGFP_KI_29' = 'rs139141690_3',
                                        'chrGFP_HET' = 'rs139141690_HET_1',
                                        'chrGFP_Del_16bp' = 'Del_16bp_1',
                                        'chrGFP_Del_233' = 'Del_80bp_1',
                                        'chrGFP_Del_235' = 'Del_80bp_2',
                                        'chrGFP_Del_287' = 'Del_80bp_3'))
  
  
  adata_sub@meta.data$Genotype<-NA
  
  
  adata_sub@meta.data$Genotype[which(adata_sub@meta.data$clone_line%in%c('wt_1','wt_2','wt_3'))]<-"wt"
  adata_sub@meta.data$Genotype[which(adata_sub@meta.data$clone_line%in%c('rs139141690_1','rs139141690_2','rs139141690_3'))]<-"rs139141690"
  adata_sub@meta.data$Genotype[which(adata_sub@meta.data$clone_line%in%c('rs139141690_HET_1'))]<-"rs139141690_HET"
  adata_sub@meta.data$Genotype[which(adata_sub@meta.data$clone_line%in%c('Del_16bp_1'))]<-"Del_16bp"
  adata_sub@meta.data$Genotype[which(adata_sub@meta.data$clone_line%in%c('Del_80bp_1','Del_80bp_2','Del_80bp_3'))]<-"Del_80bp"
  
  
  
  #### Save object  ------------------
  
  log_info_simple("save Seurat object")
  
  setwd(out)
  

  saveRDS(adata_sub,file="merged_clusters_final.rds")
  
  
 
}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--db_filt_clustered_QCed_genotyped"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--frag_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--res_param"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--processors"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--total_memory"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  cluster_macs2_simba(opt)
 

}


###########################################################################

system.time( main() )
  