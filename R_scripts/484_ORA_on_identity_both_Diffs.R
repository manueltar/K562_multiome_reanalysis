
.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/GSEA/lib/R/library"))
.libPaths()

suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(tidyverse))
suppressMessages(library(ggupset))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library('org.Hs.eg.db'))
suppressMessages(library(DOSE))
suppressMessages(library("splitstackshape"))
suppressMessages(library('data.table'))
suppressMessages(library("optparse"))
suppressMessages(library('cowplot'))
suppressMessages(library(future))



opt = NULL

options(warn = 1)


multiVals <- function(x) paste(x,collapse=";")

log_info_simple <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(timestamp, "INFO:", message, "\n")
}

ORA_function = function(option_list)
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
  
  #### READ and transform Diff_sel ----
  
  Diff_sel = opt$Diff_sel
  
  cat("Diff_sel_\n")
  cat(sprintf(as.character(Diff_sel)))
  cat("\n")
  

  
  out_path <- paste0(out,'/','ORA','_','background_adapted','_',Diff_sel,'/') # output path, where you want your results exported to

  if(file.exists(out_path)){

    unlink(out_path, recursive =T)

    dir.create(out_path)
  }else{

    dir.create(out_path)
  }


  cat("out_path_\n")
  cat(sprintf(as.character(out_path)))
  cat("\n")

  
  
  #### READ and transform path_to_GMT ----
  
  path_to_GMT = opt$path_to_GMT
  
  cat("path_to_GMT_0\n")
  cat(sprintf(as.character(path_to_GMT)))
  cat("\n")
  
  #### READ and transform search_terms ----
  
  search_terms = unlist(strsplit(opt$search_terms, split=","))
  
  cat("search_terms_\n")
  cat(sprintf(as.character(search_terms)))
  cat("\n")
  
  
  #### READ and transform pval_threshold ----
  
  pval_threshold = opt$pval_threshold
  
  cat("pval_threshold_\n")
  cat(sprintf(as.character(pval_threshold)))
  cat("\n")
  
  #### READ and transform log2FC_threshold ----
  
  log2FC_threshold = opt$log2FC_threshold
  
  cat("log2FC_threshold_\n")
  cat(sprintf(as.character(log2FC_threshold)))
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
  
  #### READ and transform DE_results ----
  
  
  DE_results<-readRDS(file=opt$DE_results)
  
  cat("DE_results_0\n")
  cat(str(DE_results))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(DE_results$contrast))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DE_results$contrast)))))
  cat("\n")
  
  
  ### Filter out NA padj -------------------------------
  
  DE_results_NO_NA<-DE_results[!is.na(DE_results$padj),]
  
  cat("DE_results_NO_NA_0\n")
  cat(str(DE_results_NO_NA))
  cat("\n")
  
  DE_results_NO_NA$ENTREZID <- mapIds(org.Hs.eg.db, keys=DE_results_NO_NA$gene, keytype="SYMBOL",
                                          column="ENTREZID", multiVals=multiVals)
  
  cat("DE_results_NO_NA_1\n")
  str(DE_results_NO_NA)
  cat("\n")
  
  DE_results_NO_NA_with_ENTREZ<-DE_results_NO_NA[-which(DE_results_NO_NA$ENTREZID == "NA"),]
  
  cat("DE_results_NO_NA_with_ENTREZ_0\n")
  str(DE_results_NO_NA_with_ENTREZ)
  cat("\n")
  
  ############ create  diffexpressed column and prepare dge df -----------------
  
  
  DE_results_NO_NA_with_ENTREZ <- DE_results_NO_NA_with_ENTREZ %>% mutate(diffexpressed = case_when(
    log2FoldChange > log2FC_threshold & padj < pval_threshold ~ 'UP',
    log2FoldChange < log2FC_threshold & padj < pval_threshold ~ 'DOWN',
    padj > 0.05 ~ 'NO'
  ))
  
  cat("DE_results_NO_NA_with_ENTREZ_1\n")
  cat(str(DE_results_NO_NA_with_ENTREZ))
  cat("\n")
  
  ## Get the genes that are present in your dataframe -----------------------------------------
  
  dataset_genes_df<-as.data.frame(unique(DE_results_NO_NA_with_ENTREZ$ENTREZID), stringsAsFactors=F)
  colnames(dataset_genes_df)<-'ENTREZID'
  
  cat("dataset_genes_df_0\n")
  cat(str(dataset_genes_df))
  cat("\n")
  
  gmt_files <- list.files(path = path_to_GMT, pattern = '.gmt$', full.names = TRUE)
  
  cat("gmt_files_0\n")
  cat(str(gmt_files))
  cat("\n")
  
  ####### LOOP 1 save all the background files a-----------------------------------------------------------------------------------
  
  DEBUG<-0
  
  universe_set<-'c5.all.v2024.1.Hs.entrez.gmt'
  
  for (iteration_gmt_files in 1:length(gmt_files)){
    
    file <- gmt_files[iteration_gmt_files]
    filename<-gsub(paste0(path_to_GMT,'/'),"",file)
    
    cat("------------------------------------------->\t")
    cat(sprintf(as.character(filename)))
    cat("\n")
    
    
    complete_collection <- read.gmt(file)
    
    if(DEBUG ==1){
      cat("complete_collection_0\n")
      cat(str(complete_collection))
      cat("\n")
    }
    
    
    
    complete_collection_gene_sel <- complete_collection[which(complete_collection$gene %in% dataset_genes_df$ENTREZID),]
    
    if(DEBUG ==1){
      cat("complete_collection_gene_sel_0\n")
      cat(str(complete_collection_gene_sel))
      cat("\n")
    }
    
    
    ### set universe -------------------
    
    
    FLAG_universe_set<-sum(filename == universe_set)
    
    if(FLAG_universe_set == 1){
      
      universe<-complete_collection_gene_sel
      
      filename_universe<-gsub("\\.gmt$","_universe.rds",filename)
    
      ## SAVE universe
      
      setwd(out_path)
      
      saveRDS(universe, file=filename_universe) 
      
    }else{
      
      
    }#FLAG_universe_set == 1
    
    
    FLAG_custom<-length(grep("Custom_Soranzo_",filename))
    
    cat("FLAG_custom\t")
    cat(sprintf(as.character(FLAG_custom)))
    cat("\n")
    
    if(FLAG_custom == 0){
      
      
      matches <- grep(paste(search_terms,collapse="|"),complete_collection_gene_sel$term)
      
      toMatch<-tolower(search_terms)
      
      
      matches_lc <- grep(paste(toMatch,collapse="|"),complete_collection_gene_sel$term)
      
      total_matches<-unique(c(matches,matches_lc))               
      
      
      if(length(total_matches) > 1){
        
        complete_collection_gene_sel_term_sel<-droplevels(complete_collection_gene_sel[total_matches,])
        
        
        if(DEBUG ==1){
          
          cat("complete_collection_gene_sel_term_sel_0\n")
          str(complete_collection_gene_sel_term_sel)
          cat("\n")
        }
        
        
        filename<-gsub("\\.gmt$","_selected.rds",filename)
        
        
        
        #filename <- paste('test','_',i,'.rds', sep='') #paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RDS', sep = '')
        
        cat(sprintf(as.character(filename)))
        cat("\n")
        cat(sprintf(as.character(out_path)))
        cat("\n")
        
        ## SAVE 1
        
        setwd(out_path)
        
        saveRDS(complete_collection_gene_sel_term_sel, file=filename) 
        
        
        equivalence_df<-complete_collection_gene_sel_term_sel
        
        ### 2 equivalence
        
        equivalence_df$SYMBOL <- mapIds(org.Hs.eg.db, keys=equivalence_df$gene, keytype="ENTREZID",
                                        column="SYMBOL", multiVals=multiVals)
        
        
        equivalence_df<-equivalence_df[,-which(colnames(equivalence_df) == 'gene')]
        
        colnames(equivalence_df)[which(colnames(equivalence_df) == 'SYMBOL')]<-'gene'
        
        
        equivalence_df.dt<-data.table(equivalence_df, key="term")
        
        
        equivalence_df_collapsed<-as.data.frame(equivalence_df.dt[,.(genes=paste(gene, collapse=",")), by=key(equivalence_df.dt)], stringsAsFactors=F)
        
        if(DEBUG ==1){
          
          cat("equivalence_df_collapsed_0\n")
          str(equivalence_df_collapsed)
          cat("\n")
        }
        
        
        filename_2<-gsub("_selected\\.rds$","_selected_equivalence.tsv",filename)
        
        
        #filename_2 <- paste('test','_',i,'.rds', sep='') #paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RDS', sep = '')
        
        cat(sprintf(as.character(filename_2)))
        cat("\n")
        cat(sprintf(as.character(out_path)))
        cat("\n")
        
        
        setwd(out_path)
        
        write.table(equivalence_df_collapsed, file=filename_2, sep = "\t",quote=F,row.names = F) 
        
      }#length(total_matches) > 1
      
    }else{
      
      complete_collection_gene_sel_term_sel<-complete_collection_gene_sel
      
      
      if(DEBUG ==1){
        
        cat("complete_collection_gene_sel_term_sel_0\n")
        str(complete_collection_gene_sel_term_sel)
        cat("\n")
      }
      
      
      filename<-gsub("\\.gmt$","_selected.rds",filename)
      
      
      
      #filename <- paste('test','_',i,'.rds', sep='') #paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RDS', sep = '')
      
      cat(sprintf(as.character(filename)))
      cat("\n")
      cat(sprintf(as.character(out_path)))
      cat("\n")
      
      ## SAVE 1
      
      setwd(out_path)
      
      saveRDS(complete_collection_gene_sel_term_sel, file=filename) 
      
      
      equivalence_df<-complete_collection_gene_sel_term_sel
      
      ### 2 equivalence
      
      equivalence_df$SYMBOL <- mapIds(org.Hs.eg.db, keys=equivalence_df$gene, keytype="ENTREZID",
                                      column="SYMBOL", multiVals=multiVals)
      
      
      equivalence_df<-equivalence_df[,-which(colnames(equivalence_df) == 'gene')]
      
      colnames(equivalence_df)[which(colnames(equivalence_df) == 'SYMBOL')]<-'gene'
      
      
      equivalence_df.dt<-data.table(equivalence_df, key="term")
      
      
      equivalence_df_collapsed<-as.data.frame(equivalence_df.dt[,.(genes=paste(gene, collapse=",")), by=key(equivalence_df.dt)], stringsAsFactors=F)
      
      if(DEBUG ==1){
        
        cat("equivalence_df_collapsed_0\n")
        str(equivalence_df_collapsed)
        cat("\n")
      }
      
      
      filename_2<-gsub("_selected\\.rds$","_selected_equivalence.tsv",filename)
      
      
      #filename_2 <- paste('test','_',i,'.rds', sep='') #paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RDS', sep = '')
      
      cat(sprintf(as.character(filename_2)))
      cat("\n")
      cat(sprintf(as.character(out_path)))
      cat("\n")
      
      
      setwd(out_path)
      
      write.table(equivalence_df_collapsed, file=filename_2, sep = "\t",quote=F,row.names = F) 
      
    }#FLAG_custom == 0
    
   
    
    
    
 
  }#iteration_gmt_files in 1:length(gmt_files)
  
  
  #### Restrict to DE genes ----

  DE_results_NO_NA_SIG<-DE_results_NO_NA_with_ENTREZ[which(DE_results_NO_NA_with_ENTREZ$diffexpressed !='NO'), ]
  
  cat("DE_results_NO_NA_SIG_0\n")
  str(DE_results_NO_NA_SIG)
  cat("\n")
 
  
 
  
 
  #### LOOP of identities  -----
  
  #======> key option in the ORA https://github.com/YuLab-SMU/clusterProfiler/issues/283 
  
  options(enrichment_force_universe=TRUE)
  
  array_identities<-levels(DE_results_NO_NA_SIG$identity)
  
  
  cat("array_identities_0\n")
  cat(str(array_identities))
  cat("\n")
  
  DEBUG<-0
  
  ORA_df<-data.frame()

  universe_set<-'c5.all.v2024.1.Hs.entrez_universe.rds'
  
  
  for(i in 1:length(array_identities)){
    
    identity_sel<-array_identities[i]
    
    cat("-------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(identity_sel)))
    cat("\n")
    
    identity_df<-data.frame()

    
    SIG_with_ENTREZ_sel<-droplevels(DE_results_NO_NA_SIG[which(DE_results_NO_NA_SIG$identity == identity_sel),])
    
    if(DEBUG == 1){
      
      cat("DE_result_sel_0\n")
      cat(str(SIG_with_ENTREZ_sel))
      cat("\n")
    }
    
    
    # LOOP per contrast and cell type -----------------------------------------------------------------
    
    
    array_contrasts<-unique(SIG_with_ENTREZ_sel$contrast)
    
    cat("array_contrasts\n")
    str(array_contrasts)
    cat("\n")
    
    
    
    
    DEBUG<-0
    
    bg_files <- list.files(path = out_path, pattern = '_selected.rds$', full.names = TRUE)
    
    cat("bg_files_0\n")
    str(bg_files)
    cat("\n")
    
   
    
    contrast_list<-list()
    
    for(k in 1:length(array_contrasts)){
      
      contrast_df<-data.frame()
      
      
      contrast_sel<-array_contrasts[k]
      
      cat("-------------------------------------------------------------------------------------------------------------------------------------------->\t")
      cat(sprintf(as.character(k)))
      cat("\t")
      cat(sprintf(as.character(contrast_sel)))
      cat("\n")
      
      SIG_with_ENTREZ_sel_contrast_sel<-SIG_with_ENTREZ_sel[which(SIG_with_ENTREZ_sel$contrast == contrast_sel),]
      
      if(dim(SIG_with_ENTREZ_sel_contrast_sel)[1] > 0){
        
        if(DEBUG ==1){
          
          cat("SIG_with_ENTREZ_sel_contrast_sel_0\n")
          cat(str(SIG_with_ENTREZ_sel_contrast_sel))
          cat("\n")
        }
        
        # Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
        
        deg_results_list <- split(SIG_with_ENTREZ_sel_contrast_sel, SIG_with_ENTREZ_sel_contrast_sel$diffexpressed)
        
        if(DEBUG ==1){
          
          cat("deg_results_list_0\n")
          cat(str(deg_results_list))
          cat("\n")
        }
        
        
        
        for(iteration_bg_files in 1:length(bg_files)){
          
          
          
          selected_collection<-bg_files[iteration_bg_files]
          
          cat("--------------------------------------------------------------------->\t")
          cat(sprintf(as.character(iteration_bg_files)))
          cat("\t")
          cat(sprintf(as.character(selected_collection)))
          cat("\n")
          
          selected_collection_df<-readRDS(file=selected_collection)
          
          if(DEBUG ==1){
            
            cat("selected_collection_df_0\n")
            # str(selected_collection_df)
            cat("\n")
          }
          
          setwd(out_path)
          
          universe_df<-readRDS(file=universe_set)
          universe<-unique(universe_df$gene)
          
            
          cat("universe_set_of_unique_genes\n")
          str(unique(universe_df$gene))
          cat("\n")
        
          
          
          FLAG_custom<-length(grep("Custom_Soranzo_", selected_collection))
          
          
          
          cat("FLAG_custom_0\n")
          str(FLAG_custom)
          cat("\n")
          
          if(FLAG_custom == 0){
            
            minGSSize_spec<-10
            maxGSSize_spec<-500
            DEBUG<-0
            
            
          }else{
            
            setwd(out_path)
            
            
           
            minGSSize_spec<-1
            maxGSSize_spec<-500
            
            DEBUG<-1
            
            if(DEBUG ==1){
              
              cat("selected_collection_df_0\n")
              str(selected_collection_df)
              cat("\n")
            }
            # universe<-unique(selected_collection_df$gene)
            
          }
          
          
          
           #### Key function of ORA -----------------------------------------------
          

          
          res<-lapply(names(deg_results_list),
                      function(x) enricher(gene= deg_results_list[[x]]$ENTREZID, 
                                           TERM2GENE = selected_collection_df, 
                                           universe=universe,
                                           maxGSSize=maxGSSize_spec,
                                           minGSSize=minGSSize_spec,
                                           pvalueCutoff = 0.05,
                                           pAdjustMethod = "BH",
                                           qvalueCutoff = 0.25))
          
          names(res) <- names(deg_results_list)
            
            
            if(DEBUG ==1){
              
              cat("res_0\n")
              # str(res)
              cat("\n")
            }
          
          ### Filter null results ----------------------
          
          filtered<-Filter(Negate(is.null), res)
          
          if(DEBUG ==1){
            cat("filtered_0\n")
            #str(filtered)
            cat("\n")
          }
          
          if(length(filtered) > 0){
            
            res_df <- lapply(names(filtered), 
                             function(x) rbind(filtered[[x]]@result))
            
            if(DEBUG ==1){
              
              cat("res_df_0\n")
              str(res_df)
              cat("\n")
            }
            
            names(res_df) <- names(filtered)
            
            if(DEBUG ==1){
              
              cat("res_df_1\n")
              str(res_df)
              cat("\n")
            }
            
            res_df <- do.call(rbind, res_df)
            
            if(DEBUG ==1){
              
              cat("res_df_2\n")
              str(res_df)
              cat("\n")
            }
            
            res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                                        diffexpressed = gsub('\\..+$', '', rownames(res_df)))                    
            
            if(DEBUG ==1){
              
              cat("res_df_3\n")
              str(res_df)
              cat("\n")
            }
            
            
            FLAG_ZERO<-length(which(unique(res_df$geneID) == ""))
            
            
            cat("FLAG_ZERO_0\n")
            str(FLAG_ZERO)
            cat("\n")
            
            
            if(FLAG_ZERO == 0){
              
              res_df_long<-unique(as.data.frame(cSplit(res_df,sep = '/', direction = "long",
                                                       splitCols = "geneID"),stringsAsFactors=F))
              
              res_df_long$geneID<-as.character(res_df_long$geneID)
              
              if(DEBUG ==1){
                str(res_df_long)
                cat("\n")
              }
              
              res_df_long$geneID <- mapIds(org.Hs.eg.db, keys=res_df_long$geneID, keytype="ENTREZID",
                                           column="SYMBOL", multiVals=multiVals)
              
              if(DEBUG ==1){
                str(res_df_long)
                cat("\n")
              }
              
              res_df_long.dt<-data.table(res_df_long, key=colnames(res_df_long)[-which(colnames(res_df_long) == 'geneID')])
              
              
              res_df_long_collapsed<-as.data.frame(res_df_long.dt[,.(geneID=paste(geneID, collapse='/')), by=key(res_df_long.dt)], stringsAsFactors=F)
              
              if(DEBUG ==1){
                str(res_df_long_collapsed)
                cat("\n")
              }
              
              
              
              contrast_df<-rbind(res_df_long_collapsed, contrast_df)     
              
            }else{
              
              contrast_df<-rbind(res_df, contrast_df)          
              
              
              
            }# FLAG_ZERO == 0
            
            
            
          }# length(filtered) > 0     
          
          
          
        }#iteration_bg_files in 1:length(bg_files)
        
        
        
        if(dim(contrast_df)[1] >0){
          
          contrast_df$contrast<-contrast_sel        
          
          
          if(DEBUG ==1){
            cat("contrast_df-----------------------------------------------\n")
            str(contrast_df)
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(contrast_df$contrast))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(contrast_df$contrast)))))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(contrast_df$ID))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(contrast_df$ID)))))
            cat("\n")
          }
          
          
          identity_df<-rbind(contrast_df,identity_df)
          
          
        }# dim(contrast_df)[1] >0    
        
        
        
      }# dim(SIG_with_ENTREZ_sel_contrast_sel)[1] > 0
      
    }# k in 1:length(array_contrasts)
    
    
    if(dim(identity_df)[1] >0){             
      
      identity_df$identity<-identity_sel
      
      if(DEBUG ==1){
        cat("identity_df\n")
        str(identity_df)
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor(identity_df$identity))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor(identity_df$identity)))))
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor(identity_df$contrast))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor(identity_df$contrast)))))
        cat("\n")
      }
      
      ORA_df<-rbind(identity_df,ORA_df)
      
    }#dim(identity_df)[1] >0
    
  }# i in 1:length(array_identities)
  
  
  if(dim(ORA_df)[1] >0){
    
    cat("GSEA_df_0\n")
    str(ORA_df)
    cat("\n")
    
    ORA_df_SIG<-ORA_df[which(ORA_df$p.adjust <= pval_threshold),]

    cat("GSEA_df_SIG_0\n")
    str(ORA_df_SIG)
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(ORA_df_SIG$identity))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ORA_df_SIG$identity)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(ORA_df_SIG$contrast))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ORA_df_SIG$contrast)))))
    cat("\n")
    
    ORA_df_SIG$contrast<-factor(ORA_df_SIG$contrast,
                                 levels=rev(levels(DE_results_NO_NA_SIG$contrast)),
                                 ordered=T)
    
    ORA_df_SIG$identity<-factor(ORA_df_SIG$identity,
                                 levels=rev(levels(DE_results_NO_NA_SIG$identity)),
                                 ordered=T)
    
    ORA_df_SIG$minuslog10padj<--log10(ORA_df_SIG$p.adjust)
    
    cat("GSEA_df_SIG_1\n")
    str(ORA_df_SIG)
    cat("\n")
    cat(sprintf(as.character(names(summary(ORA_df_SIG$identity)))))
    cat("\n")
    cat(sprintf(as.character(summary(ORA_df_SIG$identity))))
    cat("\n")
    cat(sprintf(as.character(names(summary(ORA_df_SIG$contrast)))))
    cat("\n")
    cat(sprintf(as.character(summary(ORA_df_SIG$contrast))))
    cat("\n")
    
    ##### SAVE RESULTS ----------------------------------
    
    
    setwd(out)
    
    
    write.table(ORA_df_SIG, file=paste("ORA_results_significant_",Diff_sel,".tsv",sep=''), sep="\t", quote=F, row.names = F)
    

    saveRDS(ORA_df_SIG, file=paste("ORA_results_significant_",Diff_sel,".rds",sep=''))
    
    
    
  }#dim(ORA_df)[1] >0
   
}


Lolliplot_and_gene_annotation = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
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
  
  #### READ and transform Diff_sel ----
  
  Diff_sel = opt$Diff_sel
  
  cat("Diff_sel_\n")
  cat(sprintf(as.character(Diff_sel)))
  cat("\n")
  
  #### READ and transform Threshold_number_of_genes ----
  
  Threshold_number_of_genes = opt$Threshold_number_of_genes
  
  cat("Threshold_number_of_genes_\n")
  cat(sprintf(as.character(Threshold_number_of_genes)))
  cat("\n")
  
  

  ORA_result<-readRDS(file=opt$ORA_result)
  
  cat("ORA_result_0\n")
  str(ORA_result)
  cat("\n")
  
  array_identities<-levels(ORA_result$identity)
  
  
  cat("array_identities_0\n")
  cat(str(array_identities))
  cat("\n")
  
  array_contrasts<-levels(ORA_result$contrast)
  
  
  cat("array_contrasts_0\n")
  cat(str(array_contrasts))
  cat("\n")
  
  ORA_result$contrast<-factor(ORA_result$contrast,
                               levels = rev(array_contrasts),
                               ordered = TRUE)
  
  
  cat("ORA_result_1\n")
  str(ORA_result)
  cat("\n")
  
  # Maximize pvalue by ID, identity, contrast -----------------------------------
  
  
  ORA_result.dt<-data.table(ORA_result, key=c("ID","identity","contrast","diffexpressed"))
  
  ORA_result_MAX<-as.data.frame(ORA_result.dt[,.SD[which.max(minuslog10padj)],by=key(ORA_result.dt)], stringsAsFactors=F)
  
  cat("ORA_result_MAX_0\n")
  str(ORA_result_MAX)
  cat("\n")
  
  
  
  ############### Read TF terms -------------------------
  
  TF_terms<-unique(unlist(strsplit(opt$TF_terms, split=',')))
  
  cat("TF_terms_0\n")
  str(TF_terms)
  cat("\n")
  
  TF_terms_to_match<-paste(TF_terms, collapse="|")
  
  cat("TF_terms_to_match_0\n")
  str(TF_terms_to_match)
  cat("\n")
  
  
  ## Classify pathways ------------------------------------
  
  ORA_result_MAX$PATH_class<-NA
  
  ORA_result_MAX$PATH_class[grep(TF_terms_to_match,ORA_result_MAX$ID)]<-'TF_targets'
  
  ORA_result_MAX$PATH_class[-grep(TF_terms_to_match,ORA_result_MAX$ID)]<-'other'
  
  ORA_result_MAX$PATH_class<-factor(ORA_result_MAX$PATH_class, levels=c('TF_targets','other'), ordered=T)
  
  
  cat("ORA_result_MAX_1\n")
  str(ORA_result_MAX)
  cat("\n")
  
  ## Annotate genes in pathways ------------------------------------------
  
  indx<-c(which(colnames(ORA_result_MAX)%in%c('ID','diffexpressed','identity','geneID','PATH_class')))
  
  ORA_result_MAX_sub<-unique(ORA_result_MAX[,indx])
  
  cat("ORA_result_MAX_sub_0\n")
  str(ORA_result_MAX_sub)
  cat("\n")
  
  ORA_result_MAX_sub_long<-unique(as.data.frame(cSplit(ORA_result_MAX_sub,sep = '/', direction = "long",
                                                       splitCols = "geneID") , stringsAsFactors=F))
  
  colnames(ORA_result_MAX_sub_long)[which(colnames(ORA_result_MAX_sub_long) == 'geneID')]<-'gene'
  
  cat("ORA_result_MAX_sub_long_0\n")
  str(ORA_result_MAX_sub_long)
  cat("\n")
  
  ORA_result_MAX_sub_long.dt<-data.table(ORA_result_MAX_sub_long, key=c('gene','identity','diffexpressed','PATH_class'))
  
  gene_annotation<-as.data.frame(ORA_result_MAX_sub_long.dt[,.(string_ID=paste(unique(ID), collapse='|')), by=key(ORA_result_MAX_sub_long.dt)], stringsAsFactors=F)
  
  cat("gene_annotation_0\n")
  str(gene_annotation)
  cat("\n")
  
  gene_annotation_wide<-as.data.frame(pivot_wider(gene_annotation, id_cols=c('gene','identity','diffexpressed'), names_from=PATH_class, values_from=string_ID), stringsAsFactors=F)
  
  cat("gene_annotation_wide_0\n")
  str(gene_annotation_wide)
  cat("\n")
  
  setwd(out)
  
  write.table(gene_annotation_wide, file=paste("genes_ORA_annotated_",Diff_sel,".tsv",sep=''), sep="\t", quote=F, row.names = F)
  
  ### Lolliplot -----------------------
  
  ORA_result_MAX$PATH_class<-factor(as.character(ORA_result_MAX$PATH_class), levels=rev(c('TF_targets','other')), ordered=T)
  
  
  cat("ORA_result_MAX_REMEMBER\n")
  str(ORA_result_MAX)
  cat("\n")
  
  

  ORA_result_MAX_Thresholded<-ORA_result_MAX[which(ORA_result_MAX$Count >= Threshold_number_of_genes),]
  
  
  cat("ORA_result_MAX_Thresholded_0\n")
  str(ORA_result_MAX_Thresholded)
  cat("\n")
  
  
  
  #### LOOP of identities  -----
  
  
  DEBUG<-1
 
  
  for(i in 1:length(array_identities)){
    
    identity_sel<-array_identities[i]
    
    cat("-------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(identity_sel)))
    cat("\n")
    
    path_identity_sel<-paste(out,gsub("\\/","_",gsub("\\s+","_",identity_sel)),'_',Diff_sel,'/',sep='')
    
    if (file.exists(path_identity_sel)){
      
      
      
    }else{
      
      dir.create(path_identity_sel)
    }
    
    
    ORA_result_MAX_Thresholded_sel<-ORA_result_MAX_Thresholded[which(ORA_result_MAX_Thresholded$identity == identity_sel),]
    
    
    if(dim(ORA_result_MAX_Thresholded_sel)[1] >0){
      
      if(DEBUG == 1){
        
        cat("ORA_result_MAX_Thresholded_sel_0\n")
        cat(str(ORA_result_MAX_Thresholded_sel))
        cat("\n")
      }
      
      array_diffexpressed<-unique(ORA_result_MAX_Thresholded_sel$diffexpressed)
      
      if(DEBUG == 1){
        cat("array_diffexpressed_0\n")
        str(array_diffexpressed)
        cat("\n")
      }
      
      if(length(array_diffexpressed) >0){
        
        for(iteration_diffexpressed in 1:length(array_diffexpressed)){
          
          direction_sel<-array_diffexpressed[iteration_diffexpressed]
          
          cat("------------------------------->\t")
          cat(sprintf(as.character(direction_sel)))
          cat("\n")
          
          color_selected<-NA
          
          if(direction_sel == 'UP'){
            
            color_selected<-'red'
            
          }else{
            
            if(direction_sel == 'DOWN'){
              
              color_selected<-'blue'
              
            }
          }#direction_sel == 'UP'
          
          if(DEBUG ==1){
            
            cat("color_selected_0\n")
            str(color_selected)
            cat("\n")
          }
          
          directionality_df<-ORA_result_MAX_Thresholded_sel[which(ORA_result_MAX_Thresholded_sel$diffexpressed == direction_sel),]
          
          
          
          
          
          if(dim(directionality_df)[1] >0){
            
            directionality_df<-directionality_df[order(directionality_df$PATH_class),]
            
            levels_ID<-unique(as.character(directionality_df$ID))
            
            directionality_df$DUMMY<-factor(directionality_df$ID, levels=levels_ID, ordered=T)
            
            if(DEBUG ==1){
              
              cat("directionality_df_0\n")
              str(directionality_df)
              cat("\n")
            }
            
            
            breaks_gene_sets<-as.numeric(directionality_df$DUMMY)
            labels_gene_sets<-as.character(gsub("\\..+$","",directionality_df$DUMMY))
            
            
            ORA_lolliplot<-ggplot(data=directionality_df, 
                                  aes(y=as.numeric(DUMMY),
                                      x=minuslog10padj)) +
              geom_segment(data=directionality_df,
                           aes(y=as.numeric(DUMMY),
                               yend=as.numeric(DUMMY),
                               x=0,
                               xend=minuslog10padj),
                           color=color_selected,
                           size=0.8)+
              geom_point(size=5, stroke=1, shape=21, color=color_selected, fill="white")+
              geom_text(data=directionality_df,
                        aes(x=minuslog10padj, y=as.numeric(DUMMY), label=Count),color="black",size=2, family="sans",fontface="bold")
            
            
            ORA_lolliplot <-ORA_lolliplot+
              theme_cowplot(font_size = 2,
                            font_family = "sans")+
              facet_grid(. ~ contrast+identity, scales='free_x', space='free_x', switch="y", drop=TRUE)+
              theme( strip.background = element_blank(),
                     strip.placement = "outside",
                     strip.text = element_text(size=5,color="black", family="sans"),
                     panel.spacing = unit(0.2, "lines"),
                     panel.background=element_rect(fill="white"),
                     panel.border=element_rect(colour="white",size=0,5),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
              scale_x_continuous(name='-log10pval')+
              scale_y_continuous(name=NULL, breaks=breaks_gene_sets,
                                 labels=labels_gene_sets)+
              theme_classic()+
              theme(axis.title=element_blank(),
                    axis.title.y=element_blank(),
                    axis.title.x=element_text(size=8,color="black", family="sans"),
                    axis.text.y=element_text(size=6,color="black", family="sans", face='bold'),
                    axis.text.x=element_text(size=6,color="black", family="sans"),
                    axis.line.x = element_line(size = 0.4),
                    axis.ticks.x = element_line(size = 0.4),
                    axis.ticks.y = element_line(size = 0.4),
                    axis.line.y = element_line(size = 0.4))+
              theme(legend.title = element_blank(),
                    legend.text = element_text(size=6),
                    legend.key.size = unit(0.5, 'cm'), #change legend key size
                    legend.key.height = unit(0.5, 'cm'), #change legend key height
                    legend.key.width = unit(0.5, 'cm'), #change legend key width
                    legend.position="bottom")+
              guides(fill=guide_legend(nrow=1,byrow=TRUE))
            
            
            
            
            
            setwd(path_identity_sel)
            
            svgname<-paste(paste("Lolliplot",'ORA',direction_sel, sep='_'),".svg",sep='')
            makesvg = TRUE
            
            if (makesvg == TRUE)
            {
              ggsave(svgname, plot= ORA_lolliplot,
                     device="svg", width=13)
            }
            
            
          }# dim(directionality_df)[1] >0
          
        }# iteration_diffexpressed in 1:length(array_diffexpressed)
        
        
      }#length(array_diffexpressed) >0
      
    }#dim(ORA_result_MAX_Thresholded_sel)[1] >0
    
   
    
  }#i in 1:length(array_identities)
 
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
    make_option(c("--path_to_GMT"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--search_terms"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--pval_threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--log2FC_threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DE_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ORA_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_terms"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_number_of_genes"), type="numeric", default=NULL, 
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
  
  ORA_function(opt)
  Lolliplot_and_gene_annotation(opt)
  
  
  
}


###########################################################################

system.time( main() )