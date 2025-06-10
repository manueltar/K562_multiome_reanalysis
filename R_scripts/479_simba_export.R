
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

export_simba = function(option_list)
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
  
  
  adata<-readRDS(file=opt$merged_clusters_final_annotated)
  

  metadata<-adata@meta.data
  
  cat("metadata_0\n")
  cat(str(metadata))
  cat("\n")
  
  
  #### export for SIMBA ----
  
  log_info_simple("export for SIMBA")
  
  adata@meta.data$seurat_clusters<-as.character(adata@meta.data$seurat_clusters)
  adata@meta.data$clone_line<-as.character(adata@meta.data$clone_line)
  adata@meta.data$Genotype<-as.character(adata@meta.data$Genotype)
  adata@meta.data$time_point<-as.character(adata@meta.data$time_point)
  adata@meta.data$celltype<-as.character(adata@meta.data$celltype)
  
  
  ##### Reduce the Seurat object to h5ad with RNA counts not corrected by CellBender It doesn't work with CellBender corrected counts------
  
  log_info_simple("SIMBA RNA")
  
  DefaultAssay(adata)<-'RNA_raw'
  RNA_only<-DietSeurat(adata, assays = "RNA_raw")
  
  
  setwd(out)
  
  unlink(c("RNA.h5Seurat","RNA.h5ad"))
  
  SaveH5Seurat(RNA_only, filename = "RNA.h5Seurat")
  Convert("RNA.h5Seurat", dest = "h5ad")
  
  
  ##### Reduce the Seurat object to h5ad with ATAC_by_seurat_clusters counts not corrected by CellBender It doesn't work with CellBender corrected counts------
  
  log_info_simple("SIMBA ATAC_by_seurat_clusters")
  
  
  DefaultAssay(adata)<-'ATAC_by_seurat_clusters'
  
  
  
  
  Peaks<-Features(adata)
  
  cat("Peaks_0\n")
  cat(str(Peaks))
  cat("\n")
  
  Peaks<-unique(Peaks)
  
  cat("Peaks_1\n")
  cat(str(Peaks))
  cat("\n")
  
  
  tmp.gather<- data.frame(matrix(vector(), length(Peaks), 4,
                                 dimnames=list(c(),
                                               c("chr","start","end","name"))),
                          stringsAsFactors=F)
  
  
  tmp.gather$name<-Peaks
  tmp.gather$chr<-gsub("-.+$","",Peaks)
  
  tmp.gather$start<-gsub("^[^-]+-","",Peaks)
  tmp.gather$start<-as.integer(gsub("-.+$","",tmp.gather$start))
  tmp.gather$end<-as.integer(gsub("^[^-]+-[^-]+-","",Peaks))
  
  
  cat("tmp.gather_0\n")
  cat(str(tmp.gather))
  cat("\n")
  
  tmp.gather$chr<-factor(tmp.gather$chr,
                         levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                  "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                  "chr22","chr23","chrX","chrY"), ordered=T)
  
  tmp.gather<-tmp.gather[order(tmp.gather$chr, tmp.gather$start, decreasing = F),]
  
  cat("tmp.gather_1\n")
  cat(str(tmp.gather))
  cat("\n")
  
  
  # Peaks$chr<-gsub("-.+$","",)
  
  
  
  #### remove motifs before writing h5ad ----
  
  if (!is.null(adata[["ATAC_by_seurat_clusters"]]@motifs)) {
    cat("Stripping motif information from ATAC_by_seurat_clusters assay to prevent HDF5 write errors\n")
    adata[["ATAC_by_seurat_clusters"]]@motifs <- NULL
  }
  
  ATAC_only <- DietSeurat(adata, assays = "ATAC_by_seurat_clusters", counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, dimreducs = NULL, graphs = NULL, misc = FALSE)
  
  setwd(out)
  
  unlink(c("ATAC_by_seurat_clusters.h5Seurat","ATAC_by_seurat_clusters.h5ad"))
  
  SaveH5Seurat(ATAC_only, filename = "ATAC_by_seurat_clusters.h5Seurat")
  Convert("ATAC_by_seurat_clusters.h5Seurat", dest = "h5ad")
  
  
  write.table(tmp.gather, file="Peaks.bed", sep="\t", row.names = F,quote = F,col.names = F)
  
  
  
 
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
    make_option(c("--merged_clusters_final_annotated"), type="character", default=NULL, 
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
  
  export_simba(opt)
 

}


###########################################################################

system.time( main() )
  