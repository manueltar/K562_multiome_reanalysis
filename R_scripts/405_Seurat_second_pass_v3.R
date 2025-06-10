
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




opt = NULL

options(warn = -1)

create_the_premerged_Seurat_object = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform sample_name ----
  
  sample_name = opt$sample_name
  
  cat("sample_name_\n")
  cat(sprintf(as.character(sample_name)))
  cat("\n")
  
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
  
  path_processing_outputs = opt$path_processing_outputs
  
  cat("path_processing_outputs_0\n")
  cat(sprintf(as.character(path_processing_outputs)))
  cat("\n")
  
  intermediate_dir = opt$intermediate_dir
  
  cat("intermediate_dir_0\n")
  cat(sprintf(as.character(intermediate_dir)))
  cat("\n")
  
  snATAC_dir = opt$snATAC_dir
  
  cat("snATAC_dir_0\n")
  cat(sprintf(as.character(snATAC_dir)))
  cat("\n")
  
  crange_dir = opt$crange_dir
  
  cat("crange_dir_0\n")
  cat(sprintf(as.character(crange_dir)))
  cat("\n")
  
  premerge_dir = opt$premerge_dir
  
  cat("premerge_dir_0\n")
  cat(sprintf(as.character(premerge_dir)))
  cat("\n")
  
  fragfile = opt$fragfile
  
  cat("fragfile_0\n")
  cat(sprintf(as.character(fragfile)))
  cat("\n")
  
  ############### READ 5 kB ATAC windows per sample ------------------------
  
  filename<-paste(sample_name, '_snATAC_pipeline_job.long_fmt_mtx.txt.gz', sep='')
  
  cat("filename_\n")
  cat(sprintf(as.character(filename)))
  cat("\n")
  
  atac_lfmtx = read.table(paste0(snATAC_dir, filename))
  
  cat("atac_lfmtx_0\n")
  cat(str(atac_lfmtx))
  cat("\n")
  
  #### load preliminary_filtered seurat object ------------------------
  
  
  adata <- readRDS(file = opt$preliminary_filtered)
  
  # cat("adata_0\n")
  # cat(str(adata))
  # cat("\n")
  
  
  #### Read in cell bender counts -- filter them for previous bcs -- uses it as the main RNA modality ------------------
  
  cb = Read10X_h5(file.path(path_processing_outputs, 'cellbender_gex_seurat.h5'))
  cb_counts   <- cb$'Gene Expression'
  cb_counts   <- cb_counts[,colnames(adata)]
  
  cat("cb_counts_0\n")
  cat(str(cb_counts))
  cat("\n")
  
  adata2 = CreateSeuratObject(counts = cb_counts)
  
  adata2@meta.data$orig.ident = sample_name
  
  adata2[['percent.mt']] <- PercentageFeatureSet(adata2, pattern = '^MT-')
  
  cat("adata2_0\n")
  # cat(str(adata2))
  cat("\n")
  
  #add in previous raw RNA data as another assay (RNA_raw) for comparison
  
  DefaultAssay(adata) <- 'RNA'
  raw_rna <-  GetAssayData(object = adata, slot = "counts")
  raw_rna_assay <- CreateAssayObject(counts = raw_rna)
  adata2[['RNA_raw']] <- raw_rna_assay
  
  
  cat("adata2_1\n")
  # cat(str(adata2))
  cat("\n")
  
  
  
  ####### Add ATAC modality using previously generated 5kb windows matrix ---------------------------------
  
  FALSE %in% (colnames(adata2) %in%  atac_lfmtx$V1)
  FALSE %in% ( atac_lfmtx$V1 %in% colnames(adata2))
  
  atac_lfmtx$V1 <- factor(atac_lfmtx$V1, levels=colnames(adata2))
  reordered_lfm <- atac_lfmtx[order(atac_lfmtx$V1),]
  
  cat("atac_lfmtx_1\n")
  cat(str(atac_lfmtx))
  cat("\n")
  
  atac_sm <- with(reordered_lfm,
                  sparseMatrix(i=as.numeric(as.factor(V2)), j=as.numeric(V1),
                               x=V3, dimnames=list(levels(as.factor(V2)), levels(V1))))
  
  
  cat("atac_sm_0\n")
  cat(str(atac_sm))
  cat("\n")
  
  
  #create the new chromatin assay object and add to Seurat object
  
  atac_sm       <- atac_sm[,colnames(adata2)]
  grange.counts <- StringToGRanges(rownames(atac_sm), sep = c(':', '-'))
  grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_sm       <- atac_sm[as.vector(grange.use), ]
  suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
  seqlevelsStyle(annotations)  <- 'UCSC'
  genome(annotations)          <- 'hg38'
  
  cat("atac_sm_1\n")
  cat(str(atac_sm))
  cat("\n")
  
  suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=atac_sm, sep=c(':', '-'),
                                                       genome='hg38', fragments=fragfile,
                                                       min.cells=-1, min.features=-1,
                                                       annotation=annotations))
  adata2[['ATAC']] <- chrom_assay
  
  invisible(gc())
  
  #### Compute new  metrics for cellbender rna and 5 kb windows atac -----------------------
  
  qc <- read.table(file.path(crange_dir, 'per_barcode_metrics.csv'), sep=',', header=TRUE, stringsAsFactors=1)
  qc <- as.data.frame(qc)
  rownames(qc) <- qc$gex_barcode
  qc <- qc[Cells(adata2), 6:length(colnames(qc))]
  adata2 <- AddMetaData(adata2, qc)
  
  DefaultAssay(adata2) <- 'ATAC'
  adata2 <- TSSEnrichment(adata2)
  
  
  ############ add previous metadata ------------------------------------------
  
  old_meta = adata@meta.data
  
  colkeep = c('scDblFinder.class','scDblFinder.score','scDblFinder.weighted','scDblFinder.cxds_score',
              'scDblFinder.class_atac','scDblFinder.score_atac','scDblFinder.weighted_atac','scDblFinder.cxds_score_atac',
              'DBL_comb')
  
  adata2@meta.data = cbind(adata2@meta.data,old_meta[,colkeep])
  
  ################# Read in Amulet barcodes and add to metadata ------------------------------
  
  amures           = read.table(file.path(intermediate_dir, "Amulet_selected_bc.tsv"))
  colnames(amures) = paste0("amulet_",colnames(amures))
  amures           = amures[rownames(adata2@meta.data), ]
  adata2@meta.data = cbind (adata2@meta.data, amures)
  adata2@meta.data$doublet_amulet = adata2@meta.data$amulet_q.value <0.05
  
  
  metadata_check<-adata2[[]]
  
  cat("metadata_check_0\n")
  cat(str(metadata_check))
  cat("\n")
  
  ################ CLUSTERIZATION ----------------------------------------------------------------------------
  
  adata2 <- adata2[, unname(which( colSums(GetAssayData(adata2, slot = "counts", assay = "RNA"))!=0))]
  
  ###### SAVE -----
  
  
  
  saveRDS(adata2, file = file.path(premerge_dir,'pre_merged.rds'))
  
  
  
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
    make_option(c("--sample_name"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--fragfile"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--path_processing_outputs"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--intermediate_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--snATAC_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--crange_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--premerge_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--preliminary_filtered"), type="character", default=NULL, 
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
  
  create_the_premerged_Seurat_object(opt)

}


###########################################################################

system.time( main() )
