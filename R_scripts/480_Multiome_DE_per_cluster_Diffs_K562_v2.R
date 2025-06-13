.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis/lib/R/library"))
.libPaths()
# sessionInfo()
Sys.setenv(RETICULATE_PYTHON="/home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis/bin/python")

suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(dplyr)) 
suppressMessages(library(ggplot2)) 
suppressMessages(library(Matrix)) 
suppressMessages(library(data.table)) 
suppressMessages(library(ggpubr)) 
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(presto))
suppressMessages(library("qlcMatrix"))
suppressMessages(library("cowplot"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("plyr"))
suppressMessages(library("forcats"))
suppressMessages(library('ggeasy'))
suppressMessages(library('dplyr'))
suppressMessages(library("svglite"))
suppressMessages(library("ape"))
suppressMessages(library("ggforce"))
suppressMessages(library("tidyr"))
suppressMessages(library("edgeR"))
suppressMessages(library("apeglm"))
suppressMessages(library("DESeq2"))
suppressMessages(library("tibble")) 
suppressMessages(library(future))

library("ggrepel")

library("optparse")


opt = NULL

options(warn = 1)

refine_metadata_levels <- function(seurat_data){
  for (i in base::colnames(seurat_data@meta.data)){
    if (base::is.factor(seurat_data@meta.data[[i]])){
      base::print(base::paste("Re-evaluating levels for a factor column", i))
      base::print(
        base::paste(
          "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
      seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])  # need to drop levels of the removed values
      base::print(
        base::paste(
          "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
    }
  }
  return (seurat_data)
}


multiVals <- function(x) paste(x,collapse=";")

options(warn = 1)

log_info_simple <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(timestamp, "INFO:", message, "\n")
}

DE_function = function(option_list)
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
  
  cat("OUT_\n")
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
  
  #### READ and transform Diff_sel ----
  
  Diff_sel = opt$Diff_sel
  
  cat("Diff_sel_\n")
  cat(sprintf(as.character(Diff_sel)))
  cat("\n")
  

  #### READ and transform adata ----
  
  
  adata<-readRDS(file=opt$SeuratObject)
  
  cat("adata_0\n")
  # cat(str(adata))
  cat("\n")
  
  #### Subset the Diff to investigate----
  
  adata_sub_Diff<-subset(adata, Diff == Diff_sel)
  
  
  cat("adata_sub_Diff_0\n")
  # cat(str(adata_sub_Diff))
  cat("\n")
  
  cat("SC unused levels global drop\n")
  
  adata_sub_Diff<-refine_metadata_levels(adata_sub_Diff)
 
 
  
  ## Create a variable to correct for clone line ------------------------------------
  
  
  adata_sub_Diff$sample_id<-droplevels(interaction(adata_sub_Diff$time_point, adata_sub_Diff$clone_line, sep="_"))
  
  cat(sprintf(as.character(names(summary(adata_sub_Diff$sample_id)))))
  cat("\n")
  cat(sprintf(as.character(summary(adata_sub_Diff$sample_id))))
  cat("\n")
  
  ## Extract RNA counts corrected by cell bender but not normalized -----------------------------------------------
  
  
  matrix_RNA<-LayerData(object = adata_sub_Diff, assay = "RNA", layer = "counts")
  
  cat("matrix_RNA\n")
  cat(str(matrix_RNA))
  cat("\n")
  
  ## Extract metadata ------------------------------------------------
  
  
  metadata<-adata_sub_Diff[[]]
  
  cat(sprintf(as.character(names(summary(metadata$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata$time_point))))
  cat("\n")
  
  
  cat(sprintf(as.character(names(summary(metadata$seurat_clusters)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata$seurat_clusters))))
  cat("\n")
  
  
  cat(sprintf(as.character(names(summary(metadata$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata$Genotype))))
  cat("\n")
  
  
  cat(sprintf(as.character(names(summary(metadata$clone_line)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata$clone_line))))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(metadata$sample_id))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(metadata$sample_id)))))
  cat("\n")
  
  
  ## Create a new Seurat object with the RNA and the metadata ---------------------------------------
  
  RNA_object <- CreateSeuratObject(counts = matrix_RNA, assay = "RNA",
                                   meta.data=metadata)
  
  ##### Pseudobulk per sample_id and identity of interest -------------------------------------------------------
  
  cluster_names <- levels(metadata[,which(colnames(metadata) == 'seurat_clusters')])
  
  cat("cluster_names_0\n")
  cat(str(cluster_names))
  cat("\n")
  
  sample_names <- levels(metadata[,which(colnames(metadata) == 'sample_id')])
  
  cat("sample_names_0\n")
  cat(str(sample_names))
  cat("\n")
  
  groups <- metadata[,c(which(colnames(metadata) == 'sample_id'),which(colnames(metadata) == 'seurat_clusters'))]
  
  cat("groups_0\n")
  cat(str(groups))
  cat("\n")
  
  aggr_counts <- Seurat2PB(RNA_object, sample="sample_id", cluster="seurat_clusters")
  
  cat("aggr_counts_0\n")
  cat(str(aggr_counts))
  cat("\n")
  
  ## Populate the List of counts per identity and clone line ---------------------------------------------
  
  ## Initiate empty list
  counts_ls <- list()
  
  DEBUG<-0
  
  
  for (i in 1:length(cluster_names)) {
    
    cluster_names[i]
    
    ## Extract indexes of columns in the global matrix that match a given cluster
    column_idx <- which(tstrsplit(colnames(aggr_counts), "_cluster")[[2]] == cluster_names[i])
    
    sub_aggr<- aggr_counts[, column_idx]
    
    if(DEBUG == 1)
    {
      cat("sub_aggr_0\n")
      cat(str(sub_aggr))
      cat("\n")
      
    }
    
    ## Store corresponding sub-matrix as one element of a list
    counts_ls[[i]] <-sub_aggr
    names(counts_ls)[i] <- cluster_names[i]
    
    #break
    
  }
  
  
  #### subset and prepare metadata -----
  
  # Extract sample-level variables
  metadata_NEW <- metadata %>% 
    as.data.frame() %>% 
    dplyr::select(Genotype, clone_line, time_point, sample_id)
  
  cat("metadata_NEW_0\n")
  cat(str(metadata_NEW))
  cat("\n")
  
  # Exclude duplicated rows
  metadata_NEW <- metadata_NEW[!duplicated(metadata_NEW), ]
  
  
  cat("metadata_NEW_0.5\n")
  cat(str(metadata_NEW))
  cat("\n")
  
  
  # Rename rows
  rownames(metadata_NEW) <- metadata_NEW$sample_id
  
  cat("metadata_NEW_1\n")
  cat(str(metadata_NEW))
  cat("\n")
  
 
  cat(sprintf(as.character(names(summary(metadata_NEW$Genotype)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata_NEW$Genotype))))
  cat("\n")
  
  
  cat(sprintf(as.character(names(summary(metadata_NEW$clone_line)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata_NEW$clone_line))))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(metadata_NEW$sample_id))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(metadata_NEW$sample_id)))))
  cat("\n")
  
  
  
  #### Cell count table -----
  
  t <- table(metadata$sample_id,
             metadata$seurat_clusters)
  
  cat("t_0\n")
  cat(str(t))
  cat("\n")
  
  ### Populate metadata list  --------------------------------------
  
  

  ## Initiate empty list
  
  metadata_ls <- list()
  
  DEBUG<-0
  
  for (i in 1:length(counts_ls)) {
    
    # cat("counts_ls\n")
    # cat(sprintf(as.character(i)))
    # cat("\n")
    
    ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    df <- data.frame(seurat_clusters_sample_id = colnames(counts_ls[[i]]))
    
    if(DEBUG == 1){
      
      cat("df_0\n")
      cat(str(df))
      cat("\n")
    }
    
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$seurat_clusters_id <- tstrsplit(df$seurat_clusters_sample_id, "_cluster")[[2]]
    df$sample_id  <- tstrsplit(df$seurat_clusters_sample_id, "_cluster")[[1]]
    
    identity_sel<-unique(df$seurat_clusters_id)
    
    
    if(DEBUG == 1){
      
      cat("df_1\n")
      cat(str(df))
      cat("\n")
    }
    
    
    ## Retrieve cell count information for this cluster from global cell count table
    idx <- which(colnames(t) == unique(df$seurat_clusters_id))
    
    if(DEBUG == 1){
      
      cat("idx_0\n")
      cat(str(idx))
      cat("\n")          
    }
    
    cell_counts <- t[, idx]
    
    if(DEBUG == 1){
      
      cat("cell_counts_0\n")
      cat(str(cell_counts))
      cat("\n")          
    }
    
    
    
    
    
    
    ## Remove samples with zero cell contributing to the cluster
    
    cell_counts <- cell_counts[cell_counts > 0]
    
    if(DEBUG == 1){
      
      cat("cell_counts_1\n")
      cat(str(cell_counts))
      cat("\n")          
    }
    
    names_of_vector<-names(cell_counts)
    
    if(DEBUG == 1){
      
      cat("names_of_vector_0\n")
      cat(str(names_of_vector))
      cat("\n")          
    }
    
    
    ## Match order of cell_counts and sample_ids
    sample_order <- which(df$sample_id%in%names_of_vector)
    
    if(DEBUG == 1){
      
      cat("sample_order_0\n")
      cat(str(sample_order))
      cat("\n")          
    }
    cell_counts <- cell_counts[sample_order]
    
    if(DEBUG == 1){
      
      cat("cell_counts_2\n")
      cat(str(cell_counts))
      cat("\n")          
    }
    
    cell_counts_df<-as.data.frame(cbind(cell_counts,names(cell_counts)))
    
    colnames(cell_counts_df)<-c("cell_count","sample_id")
    
    cell_counts_df$cell_count<-as.integer(cell_counts_df$cell_count)
    
    if(DEBUG == 1){
      
      cat("cell_counts_df_0\n")
      cat(str(cell_counts_df))
      cat("\n")          
    }
    
    
    
    ## Merge cell_counts to data frame
    
    df<-merge(df,
              cell_counts_df,
              by="sample_id")
    
    if(DEBUG == 1){
      
      cat("df_2\n")
      cat(str(df))
      cat("\n")
    }
    
    
    ## Join data frame (capturing metadata specific to cluster) to generic metadata
    df <- plyr::join(df, metadata_NEW, 
                     by = intersect(names(df), names(metadata_NEW)))
    
    if(DEBUG == 1){
      
      cat("df_3\n")
      cat(str(df))
      cat("\n")
    }
    
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df) <- df$current_anot_sample_id
    
    if(DEBUG == 1){
      
      cat("df_4\n")
      cat(str(df))
      cat("\n")
    }
    
    ## Store complete metadata for cluster i in list
    metadata_ls[[i]] <- df
    names(metadata_ls)[i] <- identity_sel
    
    
    
  }
  
  
  cat("metadata_ls\n")
  cat(str(metadata_ls))
  cat("\n")
  
  
  
  
  
  ## Double-check that both lists have same names --------------------------------
  
  
  # Double-check that both lists have same names
  check_1<-all(names(counts_ls) == names(metadata_ls))
  
  cat("check_1\n")
  cat(sprintf(as.character(check_1)))
  cat("\n")
  
  
  
  #################################### LOOP DE -------------------------------------------------------------------------------------------------------
  
  
  array_identities<-levels(adata_sub_Diff@meta.data$seurat_clusters)
  
  
  cat("array_identities_0\n")
  cat(str(array_identities))
  cat("\n")
  
  array_identities<-names(counts_ls)
  
  
  cat("array_identities_1\n")
  cat(str(array_identities))
  cat("\n")
  
  
  array_identities<-names(metadata_ls)
  
  
  cat("array_identities_2\n")
  cat(str(array_identities))
  cat("\n")
  
  CLASS_identities<-NULL
  

  
  
  DEBUG<-0
  
  DE_results<-data.frame()
  
  norcounts_FINAL<-data.frame()
  
  contrasts<-vector()
  
  for(i in 1:length(array_identities)){
    
    identity_sel<-gsub("\\s+","_",array_identities[i])
    
    log_info_simple(paste("------------------------------------------------------------------------------------------------>Subsetting cluster",identity_sel, collapse=" "))
    
    
    # cat("------------------------------------------------------------------------------------------------>\t")
    # cat(sprintf(as.character(array_identities[i])))
    # cat("\t")
    # cat(sprintf(as.character(identity_sel)))
    # cat("\n")
    
    CLASS_identities[i]<-identity_sel
    
    ### Subset the metadata_ls
    
    idx<-which(names(metadata_ls) == array_identities[i])
    
    if(DEBUG == 1){
      
      cat("idx_0\n")
      cat(str(idx))
      cat("\n")
    }
    
    metadata_df<-metadata_ls[[idx]]
    
    if(DEBUG == 1){
      
      cat("metadata_df_0\n")
      cat(str(metadata_df))
      cat("\n")
    }
   
    ##################################### Cell count plot ------------------------------------------------------
    
    fill_colours<-c(brewer.pal(9, "Greens")[c(5,6,7)],
                    brewer.pal(9, "YlOrRd")[c(2)],
                    brewer.pal(9, "Reds")[c(5,6,7)],
                    brewer.pal(9, "Purples")[c(7)],
                    brewer.pal(9, "Blues")[c(4,5,6)])
    
    barplot<-ggplot(data=metadata_df,
                    aes(x=clone_line, y=cell_count,
                        fill=clone_line)) +
      geom_bar(stat="identity",colour='white')+
      scale_y_continuous(name=paste("Number of cells in each pseudobulk",sep=" "), trans='log10')+
      scale_fill_manual(values=fill_colours,
                        drop=F,
                        name="clone line")
    barplot<-barplot+
      theme_cowplot(font_size = 2,
                    font_family = "sans")+
      facet_grid(. ~ time_point, scales='free_x', space='free_x', switch="y")+
      theme_classic()+
      theme(axis.title.y=element_text(size=12, color="black", family="sans"),
            axis.title.x=element_blank(),
            axis.text.y=element_text(angle=0,size=8, color="black", family="sans"),
            axis.text.x=element_text(angle=45,size=12,vjust=1,hjust=1,color="black", family="sans"),
            axis.line.x = element_line(size = 0.2),
            axis.ticks.x = element_line(size = 0.2),
            axis.ticks.y = element_line(size = 0.2),
            axis.line.y = element_line(size = 0.2))+
      theme(legend.title = element_blank(),
            legend.text = element_text(size=8, color="black", family="sans"),
            legend.key.size = unit(0.5, 'cm'), #change legend key size
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.position="hidden")
    
    
    path_identity_sel<-paste(out,gsub("\\/","_",gsub("\\s+","_",identity_sel)),'_',Diff_sel,'/',sep='')
    
    if (file.exists(path_identity_sel)){
      
      unlink(path_identity_sel, recursive =T)
      dir.create(path_identity_sel)
      
      
    }else{
      
      dir.create(path_identity_sel)
    }
    
    
    setwd(path_identity_sel)
    
    svgname<-paste("barplot_cell_counts",".svg",sep='')
    
    ggsave(svgname,plot=barplot, device ='svg', height=5, width=5)
    
    
    ### Subset the counts list ----------------------------------------------------------------------------------
    
    
    idx<-which(names(counts_ls) == array_identities[i])
    
    if(DEBUG == 1){
      
      cat("idx_0\n")
      cat(str(idx))
      cat("\n")
    }
    
    counts_ls_subset<-counts_ls[[idx]]
    
    
    
    if(DEBUG == 1){
      
      cat("counts_ls_subset_0\n")
      cat(str(counts_ls_subset))
      cat("\n")
    }
    
    ### Reference the genotype factor and make time a continuous variable ---------------------------------------------------
    
    
    metadata_df$Genotype<-as.character(metadata_df$Genotype)
    
    FLAG_ref_present<-unique(metadata_df$Genotype[which(metadata_df$Genotype == 'wt')])
    
    cat("FLAG_ref_present_0\n")
    cat(str(FLAG_ref_present))
    cat("\n")
    
    if(length(FLAG_ref_present) >0){
      
      metadata_df$Genotype<-factor(metadata_df$Genotype)
      metadata_df$Genotype<-relevel(metadata_df$Genotype, ref='wt') ### Nothing works with ordered factors
      
      if(DEBUG == 1){
        
        cat("metadata_df_1\n")
        cat(str(metadata_df))
        cat("\n")
      }
      
      metadata_df$time_point<-as.character(metadata_df$time_point)
      
      metadata_df$time<-NA
      
      metadata_df$time[which(metadata_df$time_point == '24_hours')]<-24
      metadata_df$time[which(metadata_df$time_point == '48_hours')]<-48
      metadata_df$time[which(metadata_df$time_point == '72_hours')]<-72
      metadata_df$time[which(metadata_df$time_point == '96_hours')]<-96
      
      
      if(DEBUG == 1){
        
        cat("metadata_df_2\n")
        cat(str(metadata_df))
        cat("\n")
      }
      
      ### Second check ---------------------------------------------------
      
      row.names(metadata_df)<-metadata_df$seurat_clusters_sample_id
      
      if(DEBUG == 1){
        
        cat("metadata_df_3\n")
        cat(str(metadata_df))
        cat("\n")
      }
      
      
      check_2<-all(colnames(counts_ls_subset) == rownames(metadata_df))
      
      if(DEBUG == 1){
        cat("check_2\n")
        cat(sprintf(as.character(check_2)))
        cat("\n")
        cat(sprintf(as.character(rownames(metadata_df))))
        cat("\n")
      }
      
      #### DESeqDataSetFromMatrix ----
      
      
      dds_NEW <- DESeqDataSetFromMatrix(counts_ls_subset, 
                                        colData = metadata_df, 
                                        design = ~ Genotype +time)
      
      
      ########### estimateSizeFactors----------------------------------
      
      
      dds_NEW<-estimateSizeFactors(dds_NEW)
      
      #### Filter out low count genes more than 5 counts in 3 samples------------------------------
      
      cat("PRE_filter_genes:\n")
      cat(sprintf(as.character(dim(dds_NEW)[1])))
      cat("\n")
      
      
      idx <- rowSums( counts(dds_NEW, normalized=TRUE) >= 10 ) >= 3
      
      
      dds_NEW <- dds_NEW[idx,]
      
      
      cat("POST_filter_genes:\n")
      cat(sprintf(as.character(dim(dds_NEW)[1])))
      cat("\n")
      
      ##### LRT comparing to the reduced model ---------------------------------------
      
      dds_NEW_lrt <- DESeq(dds_NEW, test = "LRT", reduced = ~ time)
      
      #### Wald test to extract specific contrasts -------------------------------------
      
      Results_per_identity<-data.frame()
      
      possible_contrasts<-colnames(dds_NEW_lrt@modelMatrix)[-1] # -1 because 1 is the Intercept term
      possible_contrasts<-possible_contrasts[-which(possible_contrasts == 'time')] # we are not interested in the effect of time alone
      
      
      for(iteration_contrasts in 1:length(possible_contrasts)){
        
        contrast_sel<-possible_contrasts[iteration_contrasts]
        
        cat("------------------------------------->\t")
        cat(sprintf(as.character(contrast_sel)))
        cat("\n")
        
        
        ############ obtain the result for this contrast with its specific pvalue coming from Wald test see ?results
        
        tmp_results<-results(dds_NEW_lrt, test= "Wald", name=contrast_sel, independentFiltering=FALSE)
        
        #### expand LogFC #########################################
        
        tmp_results <- lfcShrink(dds_NEW_lrt, 
                                 coef = contrast_sel,
                                 res=tmp_results,
                                 type = "apeglm")
        
        
        #### obtain data frame #########################################
        
        tmp_tb <- as.data.frame(tmp_results %>%
                                  data.frame() %>%
                                  rownames_to_column(var = "gene") %>%
                                  as_tibble() %>%
                                  arrange(padj), stringsAsFactors=F)
        
        colnames(tmp_results)[which(colnames(tmp_results) == 'gene')]<-'Symbol'
        
        tmp_tb$contrast<-contrast_sel
        
        
        Results_per_identity<-rbind(tmp_tb,Results_per_identity)
        
        
        
      }#iteration_contrasts in 1:length(possible_contrasts)
      
      
      Results_per_identity$identity<-identity_sel
      
      Results_per_identity$minuslog10padj <- -log10(Results_per_identity$padj)
      
      Results_per_identity$abslogfc<-abs(Results_per_identity$log2FoldChange)
      
      
      DE_results<-rbind(Results_per_identity, DE_results)
      
      #### Extract normalized counts to plot heatmaps later ----------------------------------------
      
      
      nor_counts<-as.data.frame(counts(dds_NEW, normalized=TRUE))
      
      nor_counts$gene<-row.names(nor_counts)
      row.names(nor_counts)<-NULL
      
      if(DEBUG == 1){
        cat("nor_counts_0\n")
        str(nor_counts)
        cat("\n")
      }
      
      
      
      nor_counts.m<-reshape2::melt(nor_counts, id.vars='gene', variable.name="sample_id_seurat_identity_id", value.name="count")
      nor_counts.m$identity<-gsub("\\s+","_",gsub("^.+_cluster","",nor_counts.m$sample_id_seurat_identity_id))
      nor_counts.m$clone_line_string<-gsub("_cluster.+$","",nor_counts.m$sample_id_seurat_identity_id)    
      
      nor_counts.m$time_point<-paste(gsub("_hours.+$","",nor_counts.m$clone_line_string),"_hours",sep="")
      nor_counts.m$clone_line<-gsub(paste(unique(nor_counts.m$time_point), collapse="|"), "", nor_counts.m$clone_line_string)
      nor_counts.m$clone_line<-gsub("^_","",nor_counts.m$clone_line)
      
      nor_counts.m<-nor_counts.m[,-which(colnames(nor_counts.m)%in%c('sample_id_seurat_identity_id','clone_line_string'))]
      
      
      if(DEBUG == 1){
        cat("nor_counts.m_0\n")
        str(nor_counts.m)
        cat("\n")
      }
      
      
      norcounts_FINAL<-rbind(nor_counts.m,norcounts_FINAL)
      
      
      
      
      
      
      
    }else{
      
      cat("---------------------------->skip for lack of ref cells\n")
      
    }#length(FLAG_ref_present) >0
  }#i in 1:length(array_identities)
  
  
  cat("DE_results_0\n")
  cat(str(DE_results))
  cat("\n")
  
  
  DE_results$identity<-factor(DE_results$identity,
                              levels=CLASS_identities,
                              ordered=TRUE)
  
  
  DE_results$contrast<-factor(DE_results$contrast,
                              levels=c('Genotype_rs139141690_HET_vs_wt','Genotype_rs139141690_vs_wt','Genotype_Del_16bp_vs_wt','Genotype_Del_80bp_vs_wt','time'),
                              ordered=TRUE)
  

  
  cat("DE_results_1\n")
  cat(str(DE_results))
  cat("\n")
  
  
  cat("norcounts_FINAL_0\n")
  cat(str(norcounts_FINAL))
  cat("\n")

  norcounts_FINAL$identity<-factor(norcounts_FINAL$identity,
                              levels=CLASS_identities,
                              ordered=TRUE)


  norcounts_FINAL$clone_line<-factor(norcounts_FINAL$clone_line,
                              levels=levels(metadata$clone_line),
                              ordered=TRUE)
  
  norcounts_FINAL$time_point<-factor(norcounts_FINAL$time_point,
                                     levels=levels(metadata$time_point),
                                     ordered=TRUE)


  cat("norcounts_FINAL_1\n")
  cat(str(norcounts_FINAL))
  cat("\n")

  ##### SAVE RESULTS ----------------------------------
  
  
  setwd(out)
  
  write.table(DE_results, file=paste("DE_results_",Diff_sel,".tsv",sep=''), sep="\t", quote=F, row.names = F)
  saveRDS(DE_results, file=paste("DE_results_",Diff_sel,".rds",sep=''))
  
  
  setwd(out)
  
  write.table(norcounts_FINAL, file=paste("norcounts_FINAL_",Diff_sel,".tsv",sep=''), sep="\t", quote=F, row.names = F)
  saveRDS(norcounts_FINAL, file=paste("norcounts_FINAL_",Diff_sel,".rds",sep=''))
  
  
  
  
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
    make_option(c("--SeuratObject"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Diff_sel"), type="character", default=NULL, 
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
                metavar="filename", 
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
  
  DE_function(opt)
  
  
  
  
}


###########################################################################

system.time( main() )