
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


GSEA_function = function(option_list)
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

  

  out_path <- paste0(out,'/','GSEA','_','background_adapted','_',Diff_sel,'/') # output path, where you want your results exported to

  if(file.exists(out_path)){

    unlink(out_path, recursive =T)

    dir.create(out_path)
  }else{

    dir.create(out_path)
  }


  cat("out_path_\n")
  cat(sprintf(as.character(out_path)))
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
  
  ## Get the genes that are present in your dataframe -----------------------------------------
  
  dataset_genes_df<-as.data.frame(unique(DE_results_NO_NA$gene), stringsAsFactors=F)
  colnames(dataset_genes_df)<-'gene'
  
  
  dataset_genes_df$ENTREZID <- mapIds(org.Hs.eg.db, keys=dataset_genes_df$gene, keytype="SYMBOL",
                                      column="ENTREZID", multiVals=multiVals)
  
  cat("dataset_genes_df_0\n")
  cat(str(dataset_genes_df))
  cat("\n")
  
  gmt_files <- list.files(path = path_to_GMT, pattern = '.gmt$', full.names = TRUE)
  
  cat("gmt_files_0\n")
  cat(str(gmt_files))
  cat("\n")
  
  ####### LOOP 1 save all the background files a-----------------------------------------------------------------------------------
  
  DEBUG<-0
  
  for (iteration_gmt_files in 1:length(gmt_files)){
    file <- gmt_files[iteration_gmt_files]
    
    cat("------------------------------------------->\t")
    cat(sprintf(as.character(file)))
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
    
    FLAG_custom<-length(grep("Custom_Soranzo_",file))
    
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
        
        filename<-gsub(paste0(path_to_GMT,'/'),"",file)
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
      
      filename<-gsub(paste0(path_to_GMT,'/'),"",file)
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
  
 
  #### LOOP of identities  -----
  
  
  array_identities<-levels(DE_results_NO_NA$identity)
  
  
  cat("array_identities_0\n")
  cat(str(array_identities))
  cat("\n")
  
  DEBUG<-0
  
  GSEA_df<-data.frame()
  Leading_edge_list<-list()
  
  START<-which(array_identities == 'early_MK')
  START<-1
  
  
  for(i in START:length(array_identities)){
    
    identity_sel<-array_identities[i]
    
    cat("----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(identity_sel)))
    cat("\n")
    
    identity_df<-data.frame()

    
    DE_result_sel<-droplevels(DE_results_NO_NA[which(DE_results_NO_NA$identity == identity_sel),])
    
    if(DEBUG == 1){
      
      cat("DE_result_sel_0\n")
      cat(str(DE_result_sel))
      cat("\n")
    }
    
    if(dim(DE_result_sel)[1] >0){
      # Prepare ALL genes results -----------------------------------------------
      
      
      
      DE_result_sel$ENTREZID <- mapIds(org.Hs.eg.db, keys=DE_result_sel$gene, keytype="SYMBOL",
                                       column="ENTREZID", multiVals=multiVals)
      
      if(DEBUG == 1){
        cat("DE_result_sel_2\n")
        str(DE_result_sel)
        cat("\n")
      }
      
      
      DE_result_sel_with_ENTREZ<-DE_result_sel[-which(DE_result_sel$ENTREZID == "NA"),]
      
      if(DEBUG == 1){
        cat("DE_result_sel_with_ENTREZ_0\n")
        str(DE_result_sel_with_ENTREZ)
        cat("\n")
      }
      
      
      # LOOP per contrast and cell type -----------------------------------------------------------------
      
      
      array_contrasts<-unique(DE_result_sel_with_ENTREZ$contrast)
      
      if(DEBUG == 1){
        cat("array_contrasts\n")
        str(array_contrasts)
        cat("\n")
      }
      
      
      
      
      
      
      bg_files <- list.files(path = out_path, pattern = '_selected.rds$', full.names = TRUE)
      
      if(DEBUG == 1){
        cat("bg_files_0\n")
        str(bg_files)
        cat("\n")
      }
      
      
      contrast_list<-list()
      
      for(k in 1:length(array_contrasts)){
        
        contrast_df<-data.frame()
        
        
        contrast_sel<-array_contrasts[k]
        
        cat("-------------------------------------------------------------------------------------------------------------------------------------------->\t")
        cat(sprintf(as.character(k)))
        cat("\t")
        cat(sprintf(as.character(contrast_sel)))
        cat("\n")
        
        DE_result_sel_with_ENTREZ_sel<-DE_result_sel_with_ENTREZ[which(DE_result_sel_with_ENTREZ$contrast == contrast_sel),]
        
        if(dim(DE_result_sel_with_ENTREZ_sel)[1] > 0){
          
          if(DEBUG ==1){
            
            cat("DE_result_sel_with_ENTREZ_sel_0\n")
            cat(str(DE_result_sel_with_ENTREZ_sel))
            cat("\n")
          }
          
          
          
          
          DE_result_sel_with_ENTREZ_sel_ordered<-DE_result_sel_with_ENTREZ_sel[order(DE_result_sel_with_ENTREZ_sel$log2FoldChange, decreasing=TRUE),]
          
          if(DEBUG ==1){
            
            cat("DE_result_sel_with_ENTREZ_sel_ordered_0\n")
            cat(str(DE_result_sel_with_ENTREZ_sel_ordered))
            cat("\n")
          }
          
          
          List_bg_files<-list()
          
          
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
              str(selected_collection_df)
              cat("\n")
            }
            
            
            FLAG_custom<-length(grep("Custom_Soranzo_",selected_collection))
            
            
            
            cat("FLAG_custom_0\n")
            str(FLAG_custom)
            cat("\n")
            
            if(FLAG_custom == 0){
              
              minGSSize_spec<-10
              maxGSSize_spec<-500
              
              
            }else{
              
              setwd(out_path)
              
              
              
              minGSSize_spec<-10
              maxGSSize_spec<-500
              
              
              if(DEBUG ==1){
                
                cat("selected_collection_df_0\n")
                str(selected_collection_df)
                cat("\n")
              }
              # universe<-unique(selected_collection_df$gene)
              
            }
            
            
            
            #### Key function of GSEA -----------------------------------------------
            
            geneList<-DE_result_sel_with_ENTREZ_sel_ordered$log2FoldChange
            
            names(geneList)<-DE_result_sel_with_ENTREZ_sel_ordered$ENTREZID
            
            
            if(DEBUG ==1){
              
              cat("geneList_0\n")
              str(geneList)
              cat("\n")
            }
            
            FLAG_overlap<-length(which(names(geneList)%in%selected_collection_df$gene))
            # 
            # if(DEBUG ==1){
            
            cat("FLAG_overlap_0\n")
            str(FLAG_overlap)
            cat("\n")
            # }
            
            
            if(FLAG_overlap > 0){
              
              res = GSEA(geneList, 
                         TERM2GENE=selected_collection_df,                            
                         maxGSSize=maxGSSize_spec,
                         minGSSize=minGSSize_spec,
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH")
              
              
              
              
              if(DEBUG ==1){
                
                cat("res_0\n")
                # str(res)
                cat("\n")
              }
              
              
              res_df <- res@result
              
              if(dim(res_df)[1] >0){
                
                if(DEBUG ==1){
                  
                  cat("res_df_0\n")
                  str(res_df)
                  cat("\n")
                }
                
                res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust))
                
                if(DEBUG ==1){
                  
                  cat("res_df_3\n")
                  str(res_df)
                  cat("\n")
                }
                
                
                FLAG_ZERO<-length(which(unique(res_df$core_enrichment) == ""))
                
                
                cat("FLAG_ZERO_0\n")
                str(FLAG_ZERO)
                cat("\n")
                
                
                if(FLAG_ZERO == 0){
                  
                  res_df_long<-unique(as.data.frame(cSplit(res_df,sep = '/', direction = "long",
                                                           splitCols = "core_enrichment"),stringsAsFactors=F))
                  
                  res_df_long$core_enrichment<-as.character(res_df_long$core_enrichment)
                  
                  if(DEBUG ==1){
                    str(res_df_long)
                    cat("\n")
                  }
                  
                  res_df_long$core_enrichment <- mapIds(org.Hs.eg.db, keys=res_df_long$core_enrichment, keytype="ENTREZID",
                                                        column="SYMBOL", multiVals=multiVals)
                  
                  if(DEBUG ==1){
                    str(res_df_long)
                    cat("\n")
                  }
                  
                  res_df_long.dt<-data.table(res_df_long, key=colnames(res_df_long)[-which(colnames(res_df_long) == 'core_enrichment')])
                  
                  
                  res_df_long_collapsed<-as.data.frame(res_df_long.dt[,.(core_enrichment=paste(core_enrichment, collapse='/')), by=key(res_df_long.dt)], stringsAsFactors=F)
                  
                  if(DEBUG ==1){
                    str(res_df_long_collapsed)
                    cat("\n")
                  }
                  
                  
                  
                  contrast_df<-rbind(res_df_long_collapsed, contrast_df)     
                  
                }else{
                  
                  contrast_df<-rbind(res_df, contrast_df)          
                  
                  
                  
                }# FLAG_ZERO == 0
                
                ##########
                
                
                
                List_bg_files[[selected_collection]]<-res
                
                names(List_bg_files)<-gsub("\\.entrez_selected\\.rds$","",gsub(paste(out_path,'/',sep=""),"",names(List_bg_files)))
                
                
                
              }#dim(res_df)[1] >0
              
            }else{
              
              # do nothing
              
            }#FLAG_overlap > 0
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
            
            contrast_list[[k]]<-List_bg_files
            
            names(contrast_list)[k]<-contrast_sel
            
            
            
            
          }# dim(contrast_df)[1] >0    
          
          
          
        }# dim(DE_result_sel_with_ENTREZ_sel)[1] > 0
        
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
        
        GSEA_df<-rbind(identity_df,GSEA_df)
        
        
        Leading_edge_list[[identity_sel]]<-contrast_list
        
        
        
        
      }#dim(identity_df)[1] >0
      
      
    }#dim(DE_result_sel)[1] >0
    
    
   
  }# i in 1:length(array_identities)
  
  
  if(dim(GSEA_df)[1] >0){
    
    cat("GSEA_df_0\n")
    str(GSEA_df)
    cat("\n")
    
    GSEA_df_SIG<-GSEA_df[which(GSEA_df$p.adjust <= pval_threshold),]

    cat("GSEA_df_SIG_0\n")
    str(GSEA_df_SIG)
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(GSEA_df_SIG$identity))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(GSEA_df_SIG$identity)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(GSEA_df_SIG$contrast))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(GSEA_df_SIG$contrast)))))
    cat("\n")
    
    GSEA_df_SIG$contrast<-factor(GSEA_df_SIG$contrast,
                                 levels=rev(levels(DE_results_NO_NA$contrast)),
                                 ordered=T)
    
    GSEA_df_SIG$identity<-factor(GSEA_df_SIG$identity,
                                 levels=rev(levels(DE_results_NO_NA$identity)),
                                 ordered=T)
    
    GSEA_df_SIG$minuslog10padj<--log10(GSEA_df_SIG$p.adjust)
    
    cat("GSEA_df_SIG_1\n")
    str(GSEA_df_SIG)
    cat("\n")
    cat(sprintf(as.character(names(summary(GSEA_df_SIG$identity)))))
    cat("\n")
    cat(sprintf(as.character(summary(GSEA_df_SIG$identity))))
    cat("\n")
    cat(sprintf(as.character(names(summary(GSEA_df_SIG$contrast)))))
    cat("\n")
    cat(sprintf(as.character(summary(GSEA_df_SIG$contrast))))
    cat("\n")
    
    ##### SAVE RESULTS ----------------------------------
    
    
    setwd(out)
    
    
    write.table(GSEA_df_SIG, file=paste("GSEA_results_significant_",Diff_sel,".tsv", sep=''), sep="\t", quote=F, row.names = F)
    
    saveRDS(Leading_edge_list, file=paste("GSEA_complete_results_",Diff_sel,".rds", sep=''))
    
    saveRDS(GSEA_df_SIG, file=paste("GSEA_results_significant_",Diff_sel,".rds",sep=''))
    
    
    
  }#dim(GSEA_df)[1] >0
   
}

Leading_edge_printer = function(option_list)
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
  
  #### READ GSEA_complete_results.rds ----

  
  List_GSEA<-readRDS(file=opt$List_GSEA)
  
  cat("List_GSEA_0\n")
  # str(List_GSEA)
  cat("\n")
  
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
  
  
  #### LOOP of identities  -----
  
  
  array_identities<-names(List_GSEA)
  
  
  cat("array_identities_0\n")
  cat(str(array_identities))
  cat("\n")
  
  DEBUG<-0
  

  
  array_contrasts<-levels(DE_results_NO_NA$contrast)
  
  cat("array_contrasts\n")
  str(array_contrasts)
  cat("\n")
  
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
    
    List_GSEA_identity_sel<-List_GSEA[[identity_sel]]
    
    
    indexes_contrasts<-names(List_GSEA_identity_sel)
    
    
    cat("indexes_contrasts\n")
    str(indexes_contrasts)
    cat("\n")
    
    FLAG_NA<-sum(!is.na(indexes_contrasts) == TRUE)
    
    if(FLAG_NA >0){
      
      indexes_contrasts_NO_NA<-indexes_contrasts[!is.na(indexes_contrasts)]
      
    }else{
      
      indexes_contrasts_NO_NA<-indexes_contrasts
      
      
    }#FLAG_NA >0
    
    
    FLAG_EMPTY<-length(which(indexes_contrasts == ""))
    
    if(FLAG_EMPTY >0){
      
      indexes_contrasts_NO_NA_NO_EMPTY<-indexes_contrasts_NO_NA[which(indexes_contrasts_NO_NA != "")]
      
    }else{
      
      indexes_contrasts_NO_NA_NO_EMPTY<-indexes_contrasts_NO_NA
      
      
    }#FLAG_EMPTY >0
    
    if(length(indexes_contrasts_NO_NA_NO_EMPTY) >0){
      for(k in 1:length(indexes_contrasts_NO_NA_NO_EMPTY)){
        
        indexes_contrasts_NO_NA_NO_EMPTY_sel<-indexes_contrasts_NO_NA_NO_EMPTY[k]
        contrasts_sel<-array_contrasts[as.numeric(indexes_contrasts_NO_NA_NO_EMPTY_sel)]
        
        
        cat("----------------->\t")
        cat(sprintf(as.character(indexes_contrasts_NO_NA_NO_EMPTY_sel)))
        cat("\t")
        cat(sprintf(as.character(contrasts_sel)))
        cat("\n")
        
        
        
        path_contrast_sel<-paste(path_identity_sel,contrasts_sel,'/',sep='')
        
        if (file.exists(path_contrast_sel)){
          
          
          
        }else{
          
          dir.create(path_contrast_sel)
        }
        
        Leading_edge_plots_dir<-paste0(path_contrast_sel,'GSEA_Leading_edge_plots','/')
        
        if(file.exists(Leading_edge_plots_dir)){
          
          unlink(Leading_edge_plots_dir, recursive =T)
          
          dir.create(Leading_edge_plots_dir)
        }else{
          
          dir.create(Leading_edge_plots_dir)
        }
        
        
        cat("Leading_edge_plots_dir_0\n")
        cat(sprintf(as.character(Leading_edge_plots_dir)))
        cat("\n")
        
        
        List_GSEA_identity_sel_contrast_sel<-List_GSEA_identity_sel[[indexes_contrasts_NO_NA_NO_EMPTY_sel]]
        
        
        
        array_collections<-names(List_GSEA_identity_sel_contrast_sel)
        
        cat("array_collections_0\n")
        str(array_collections)
        cat("\n")
        
        for(h in 1:length(array_collections)){
          
          array_collections_sel<-array_collections[h]
          
          
          cat("--->\t")
          cat(sprintf(as.character(array_collections_sel)))
          cat("\n")
          
          List_GSEA_identity_sel_contrast_sel_collection_sel<-List_GSEA_identity_sel_contrast_sel[[array_collections_sel]]
          
          array_gene_sets<-List_GSEA_identity_sel_contrast_sel_collection_sel$Description
          
          cat("array_gene_sets_0\n")
          str(array_gene_sets)
          cat("\n")
          
          for(l in 1:length(array_gene_sets)){
            
            array_gene_sets_sel<-array_gene_sets[l]
            
            
            cat(">\t")
            cat(sprintf(as.character(array_gene_sets_sel)))
            cat("\n")
            
            
            
            
            Leading_edge_plot<-gseaplot(List_GSEA_identity_sel_contrast_sel_collection_sel, geneSetID = l, title = List_GSEA_identity_sel_contrast_sel_collection_sel@result$Description[l])
            
            
            
            setwd(Leading_edge_plots_dir)
            
            
            filename<-gsub("\\s+","_",array_gene_sets_sel)
            
            
            FLAG_nchar<-nchar(filename)
            
            if(FLAG_nchar >= 40){
              
              
              svgname<-paste(paste(unlist(strsplit(filename, split=''))[c(1:40)],collapse=''),'.svg', sep='')
              
              
              
            }else{
              
              
              svgname<-paste(filename,'.svg', sep='')
              
              
            }#FLAG_nchar >= 30
            
            
            
            makesvg = TRUE
            
            if (makesvg == TRUE)
            {
              ggsave(svgname, plot= Leading_edge_plot,
                     device="svg")
            }
            
            
            
            
          }#l in 1:length(array_gene_sets
        }# h in 1:length(array_collections)
        
      }# k in 1:length(indexes_contrasts_NO_NA_NO_EMPTY)
      
    }#length(indexes_contrasts_NO_NA_NO_EMPTY) >0
  }# i in 1:length(array_identities)
  
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
  
  

  GSEA_result<-readRDS(file=opt$GSEA_result)
  
  cat("GSEA_result_0\n")
  str(GSEA_result)
  cat("\n")
  
  array_identities<-levels(GSEA_result$identity)
  
  
  cat("array_identities_0\n")
  cat(str(array_identities))
  cat("\n")
  
  array_contrasts<-levels(GSEA_result$contrast)
  
  
  cat("array_contrasts_0\n")
  cat(str(array_contrasts))
  cat("\n")
  
  GSEA_result$contrast<-factor(GSEA_result$contrast,
                               levels = rev(array_contrasts),
                               ordered = TRUE)
  
  
  cat("GSEA_result_1\n")
  str(GSEA_result)
  cat("\n")
  
  # Maximize pvalue by ID, identity, contrast -----------------------------------
  
  
  GSEA_result.dt<-data.table(GSEA_result, key=c("ID","identity","contrast"))
  
  GSEA_result_MAX<-as.data.frame(GSEA_result.dt[,.SD[which.max(minuslog10padj)],by=key(GSEA_result.dt)], stringsAsFactors=F)
  
  cat("GSEA_result_MAX_0\n")
  str(GSEA_result_MAX)
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
  
  GSEA_result_MAX$PATH_class<-NA
  
  GSEA_result_MAX$PATH_class[grep(TF_terms_to_match,GSEA_result_MAX$ID)]<-'TF_targets'
  
  GSEA_result_MAX$PATH_class[-grep(TF_terms_to_match,GSEA_result_MAX$ID)]<-'other'
  
  GSEA_result_MAX$PATH_class<-factor(GSEA_result_MAX$PATH_class, levels=c('TF_targets','other'), ordered=T)
  
  
  cat("GSEA_result_MAX_1\n")
  str(GSEA_result_MAX)
  cat("\n")
  
  ## Annotate genes in pathways ------------------------------------------
  
  ind.dep<-which(colnames(GSEA_result_MAX)%in%c('core_enrichment'))
  
  
  GSEA_result_MAX.dt<-data.table(GSEA_result_MAX, key=colnames(GSEA_result_MAX)[-ind.dep])
  
  
  
  GSEA_result_MAX_with_count<-as.data.frame(GSEA_result_MAX.dt[,.(core_enrichment = core_enrichment,
                                                                  Count=length(unlist(strsplit(core_enrichment, split='/')))), by=key(GSEA_result_MAX.dt)], stringAsFactors=F)
  

  cat("GSEA_result_MAX_with_count_0\n")
  str(GSEA_result_MAX_with_count)
  cat("\n")
  
  
  
  ## Annotate genes in pathways ------------------------------------------
  
  indx<-c(which(colnames(GSEA_result_MAX_with_count)%in%c('ID','identity','core_enrichment','PATH_class')))
  
  GSEA_result_MAX_with_count_sub<-unique(GSEA_result_MAX_with_count[,indx])
  
  cat("GSEA_result_MAX_with_count_sub_0\n")
  str(GSEA_result_MAX_with_count_sub)
  cat("\n")
  
  GSEA_result_MAX_with_count_sub_long<-unique(as.data.frame(cSplit(GSEA_result_MAX_with_count_sub,sep = '/', direction = "long",
                                                        splitCols = "core_enrichment") , stringsAsFactors=F))
  
  colnames(GSEA_result_MAX_with_count_sub_long)[which(colnames(GSEA_result_MAX_with_count_sub_long) == 'core_enrichment')]<-'gene'
  
  cat("GSEA_result_MAX_with_count_sub_long_0\n")
  str(GSEA_result_MAX_with_count_sub_long)
  cat("\n")
  
  GSEA_result_MAX_with_count_sub_long.dt<-data.table(GSEA_result_MAX_with_count_sub_long, key=c('gene','identity','PATH_class'))
  
  gene_annotation<-as.data.frame(GSEA_result_MAX_with_count_sub_long.dt[,.(string_ID=paste(unique(ID), collapse='|')), by=key(GSEA_result_MAX_with_count_sub_long.dt)], stringsAsFactors=F)
  
  cat("gene_annotation_0\n")
  str(gene_annotation)
  cat("\n")
  
  gene_annotation_wide<-as.data.frame(pivot_wider(gene_annotation, id_cols=c('gene','identity'), names_from=PATH_class, values_from=string_ID), stringsAsFactors=F)
  
  cat("gene_annotation_wide_0\n")
  str(gene_annotation_wide)
  cat("\n")
  
  setwd(out)
  
  write.table(gene_annotation_wide, file=paste("genes_GSEA_annotated_",Diff_sel,".tsv",sep=''), sep="\t", quote=F, row.names = F)
  
  ### Lolliplot -----------------------
  
  GSEA_result_MAX_with_count$PATH_class<-factor(as.character(GSEA_result_MAX_with_count$PATH_class), levels=rev(c('TF_targets','other')), ordered=T)
  
  
  cat("GSEA_result_MAX_with_count_REMEMBER\n")
  str(GSEA_result_MAX_with_count)
  cat("\n")
  
  
 
  

  GSEA_result_MAX_with_count_Thresholded<-GSEA_result_MAX_with_count[which(GSEA_result_MAX_with_count$Count >= Threshold_number_of_genes),]
  
  
  cat("GSEA_result_MAX_with_count_Thresholded_REMEMBER\n")
  str(GSEA_result_MAX_with_count_Thresholded)
  cat("\n")
  
  
  
  #### LOOP of identities  -----
  
  
  DEBUG<-0
 
  
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
    
    
    GSEA_result_MAX_with_count_Thresholded_sel<-GSEA_result_MAX_with_count_Thresholded[which(GSEA_result_MAX_with_count_Thresholded$identity == identity_sel),]
    
    if(DEBUG == 1){
      
      cat("GSEA_result_MAX_with_count_Thresholded_sel_0\n")
      cat(str(GSEA_result_MAX_with_count_Thresholded_sel))
      cat("\n")
    }
    
    if(dim(GSEA_result_MAX_with_count_Thresholded_sel)[1] >0){
      
      GSEA_result_MAX_with_count_Thresholded_sel<-GSEA_result_MAX_with_count_Thresholded_sel[order(GSEA_result_MAX_with_count_Thresholded_sel$PATH_class),]
      
      levels_ID<-unique(as.character(GSEA_result_MAX_with_count_Thresholded_sel$ID))
      
      GSEA_result_MAX_with_count_Thresholded_sel$DUMMY<-factor(GSEA_result_MAX_with_count_Thresholded_sel$ID, levels=levels_ID, ordered=T)
      
      if(DEBUG ==1){
        
        cat("GSEA_result_MAX_with_count_Thresholded_sel_0\n")
        str(GSEA_result_MAX_with_count_Thresholded_sel)
        cat("\n")
      }
      
      
      breaks_gene_sets<-as.numeric(GSEA_result_MAX_with_count_Thresholded_sel$DUMMY)
      labels_gene_sets<-as.character(gsub("\\..+$","",GSEA_result_MAX_with_count_Thresholded_sel$DUMMY))
      
      
      Gene_set_lolliplot<-ggplot(data=GSEA_result_MAX_with_count_Thresholded_sel, 
                                 aes(y=as.numeric(DUMMY),
                                     x=minuslog10padj)) +
        geom_segment(data=GSEA_result_MAX_with_count_Thresholded_sel[which(GSEA_result_MAX_with_count_Thresholded_sel$minuslog10padj > 0), ],
                     aes(y=as.numeric(DUMMY),
                         yend=as.numeric(DUMMY),
                         x=0,
                         xend=minuslog10padj,
                         color=NES),
                     size=0.8)+
        geom_point(data=GSEA_result_MAX_with_count_Thresholded_sel[which(GSEA_result_MAX_with_count_Thresholded_sel$minuslog10padj > 0), ],
                   aes(color=NES),size=5, stroke=1, shape=21,fill="white")+
        geom_text(data=GSEA_result_MAX_with_count_Thresholded_sel[which(GSEA_result_MAX_with_count_Thresholded_sel$minuslog10padj > 0), ],
                  aes(x=minuslog10padj, y=as.numeric(DUMMY), label=Count),color="black",size=2, family="sans",fontface="bold")+
        geom_segment(data=GSEA_result_MAX_with_count_Thresholded_sel[which(GSEA_result_MAX_with_count_Thresholded_sel$minuslog10padj < 0), ],
                     aes(y=as.numeric(DUMMY),
                         yend=as.numeric(DUMMY),
                         x=0,
                         xend=minuslog10padj,
                         color=NES),
                     size=0.8)+
        geom_point(data=GSEA_result_MAX_with_count_Thresholded_sel[which(GSEA_result_MAX_with_count_Thresholded_sel$minuslog10padj < 0), ],
                   aes(color=NES), stroke=1, shape=21, fill="white")+
        geom_text(data=GSEA_result_MAX_with_count_Thresholded_sel[which(GSEA_result_MAX_with_count_Thresholded_sel$minuslog10padj < 0), ],
                  aes(x=minuslog10padj, y=as.numeric(DUMMY), label=Count),color="black",size=2, family="sans",fontface="bold")+
        scale_color_gradient2(name=paste("Normalized","Enrichment","Score", sep="\n"),
                              low = "blue", high = "red",mid="white",midpoint=0,
                              na.value = NA)
      
      
      Gene_set_lolliplot <-Gene_set_lolliplot+
        theme_cowplot(font_size = 2,
                      font_family = "sans")+
        facet_grid(. ~ identity + contrast, scales='free_x', space='free_x', switch="y", drop=TRUE)+
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
              axis.title.x=element_text(size=12,color="black", family="sans"),
              axis.text.y=element_text(size=10,color="black", family="sans", face='bold'),
              axis.text.x=element_text(size=10,color="black", family="sans"))+
        theme(legend.title = element_text(size=12),
              legend.text = element_text(size=10),
              legend.key.size = unit(0.5, 'cm'), #change legend key size
              legend.key.height = unit(0.5, 'cm'), #change legend key height
              legend.key.width = unit(0.5, 'cm'), #change legend key width
              legend.position="right")
      
      
      setwd(path_identity_sel)
      
      svgname<-paste(paste("Lolliplot",'GSEA',sep='_'),".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= Gene_set_lolliplot,
               device="svg", width=13)
      }
      
      
    }# dim(GSEA_result_MAX_with_count_Thresholded_sel)[1] >0
    
  }#i in 1:length(array_identities)
 
}

Lolliplot_and_gene_annotation_subset = function(option_list)
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
  
  
  
  GSEA_result<-readRDS(file=opt$GSEA_result)
  
  cat("GSEA_result_0\n")
  str(GSEA_result)
  cat("\n")
  
  array_identities<-levels(GSEA_result$identity)
  
  
  cat("array_identities_0\n")
  cat(str(array_identities))
  cat("\n")
  
  array_contrasts<-levels(GSEA_result$contrast)
  
  
  cat("array_contrasts_0\n")
  cat(str(array_contrasts))
  cat("\n")
  
  GSEA_result$contrast<-factor(GSEA_result$contrast,
                               levels = rev(array_contrasts),
                               ordered = TRUE)
  
  
  cat("GSEA_result_1\n")
  str(GSEA_result)
  cat("\n")
  
  # Maximize pvalue by ID, identity, contrast -----------------------------------
  
  
  GSEA_result.dt<-data.table(GSEA_result, key=c("ID","identity","contrast"))
  
  GSEA_result_MAX<-as.data.frame(GSEA_result.dt[,.SD[which.max(minuslog10padj)],by=key(GSEA_result.dt)], stringsAsFactors=F)
  
  cat("GSEA_result_MAX_0\n")
  str(GSEA_result_MAX)
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
  
  GSEA_result_MAX$PATH_class<-NA
  
  GSEA_result_MAX$PATH_class[grep(TF_terms_to_match,GSEA_result_MAX$ID)]<-'TF_targets'
  
  GSEA_result_MAX$PATH_class[-grep(TF_terms_to_match,GSEA_result_MAX$ID)]<-'other'
  
  GSEA_result_MAX$PATH_class<-factor(GSEA_result_MAX$PATH_class, levels=c('TF_targets','other'), ordered=T)
  
  
  cat("GSEA_result_MAX_1\n")
  str(GSEA_result_MAX)
  cat("\n")
  
  ## Annotate genes in pathways ------------------------------------------
  
  ind.dep<-which(colnames(GSEA_result_MAX)%in%c('core_enrichment'))
  
  
  GSEA_result_MAX.dt<-data.table(GSEA_result_MAX, key=colnames(GSEA_result_MAX)[-ind.dep])
  
  
  
  GSEA_result_MAX_with_count<-as.data.frame(GSEA_result_MAX.dt[,.(core_enrichment = core_enrichment,
                                                                  Count=length(unlist(strsplit(core_enrichment, split='/')))), by=key(GSEA_result_MAX.dt)], stringAsFactors=F)
  
  
  cat("GSEA_result_MAX_with_count_0\n")
  str(GSEA_result_MAX_with_count)
  cat("\n")
  
  
  
 
  
  ### Lolliplot -----------------------
  
  GSEA_result_MAX_with_count$PATH_class<-factor(as.character(GSEA_result_MAX_with_count$PATH_class), levels=rev(c('TF_targets','other')), ordered=T)
  
  
  cat("GSEA_result_MAX_with_count_REMEMBER\n")
  str(GSEA_result_MAX_with_count)
  cat("\n")
  
  
  
  
  
  GSEA_result_MAX_with_count_Thresholded<-GSEA_result_MAX_with_count[which(GSEA_result_MAX_with_count$Count >= Threshold_number_of_genes),]
  
  
  cat("GSEA_result_MAX_with_count_Thresholded_REMEMBER\n")
  str(GSEA_result_MAX_with_count_Thresholded)
  cat("\n")
  
  ID_selected<-paste(c('DESCARTES_MAIN_FETAL_MEGAKARYOCYTES','Cell cycle','G2M','HP_INCREASED_MEAN_PLATELET_VOLUME'), collapse="|")
  
  cat("ID_selected_0\n")
  str(ID_selected)
  cat("\n")
  
  
  GSEA_result_MAX_with_count_Thresholded_subset<-GSEA_result_MAX_with_count_Thresholded[grep(ID_selected, GSEA_result_MAX_with_count_Thresholded$ID),]
  
  cat("GSEA_result_MAX_with_count_Thresholded_subset_0\n")
  str(GSEA_result_MAX_with_count_Thresholded_subset)
  cat("\n")
  
  #### LOOP of identities  -----
  
  
  DEBUG<-0
  
  
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
    
    
    GSEA_result_MAX_with_count_Thresholded_subset_sel<-GSEA_result_MAX_with_count_Thresholded_subset[which(GSEA_result_MAX_with_count_Thresholded_subset$identity == identity_sel),]
    
    if(DEBUG == 1){
      
      cat("GSEA_result_MAX_with_count_Thresholded_subset_sel_0\n")
      cat(str(GSEA_result_MAX_with_count_Thresholded_subset_sel))
      cat("\n")
    }
    
    if(dim(GSEA_result_MAX_with_count_Thresholded_subset_sel)[1] >0){
      
      GSEA_result_MAX_with_count_Thresholded_subset_sel<-GSEA_result_MAX_with_count_Thresholded_subset_sel[order(GSEA_result_MAX_with_count_Thresholded_subset_sel$PATH_class),]
      
      levels_ID<-unique(as.character(GSEA_result_MAX_with_count_Thresholded_subset_sel$ID))
      
      GSEA_result_MAX_with_count_Thresholded_subset_sel$DUMMY<-factor(GSEA_result_MAX_with_count_Thresholded_subset_sel$ID, levels=levels_ID, ordered=T)
      
      if(DEBUG ==1){
        
        cat("GSEA_result_MAX_with_count_Thresholded_subset_sel_0\n")
        str(GSEA_result_MAX_with_count_Thresholded_subset_sel)
        cat("\n")
      }
      
      
      breaks_gene_sets<-as.numeric(GSEA_result_MAX_with_count_Thresholded_subset_sel$DUMMY)
      labels_gene_sets<-as.character(gsub("\\..+$","",GSEA_result_MAX_with_count_Thresholded_subset_sel$DUMMY))
      
      
      Gene_set_lolliplot<-ggplot(data=GSEA_result_MAX_with_count_Thresholded_subset_sel, 
                                 aes(y=as.numeric(DUMMY),
                                     x=minuslog10padj)) +
        geom_segment(data=GSEA_result_MAX_with_count_Thresholded_subset_sel[which(GSEA_result_MAX_with_count_Thresholded_subset_sel$minuslog10padj > 0), ],
                     aes(y=as.numeric(DUMMY),
                         yend=as.numeric(DUMMY),
                         x=0,
                         xend=minuslog10padj,
                         color=NES),
                     size=0.8)+
        geom_point(data=GSEA_result_MAX_with_count_Thresholded_subset_sel[which(GSEA_result_MAX_with_count_Thresholded_subset_sel$minuslog10padj > 0), ],
                   aes(color=NES),size=5, stroke=1, shape=21,fill="white")+
        geom_text(data=GSEA_result_MAX_with_count_Thresholded_subset_sel[which(GSEA_result_MAX_with_count_Thresholded_subset_sel$minuslog10padj > 0), ],
                  aes(x=minuslog10padj, y=as.numeric(DUMMY), label=Count),color="black",size=2, family="sans",fontface="bold")+
        geom_segment(data=GSEA_result_MAX_with_count_Thresholded_subset_sel[which(GSEA_result_MAX_with_count_Thresholded_subset_sel$minuslog10padj < 0), ],
                     aes(y=as.numeric(DUMMY),
                         yend=as.numeric(DUMMY),
                         x=0,
                         xend=minuslog10padj,
                         color=NES),
                     size=0.8)+
        geom_point(data=GSEA_result_MAX_with_count_Thresholded_subset_sel[which(GSEA_result_MAX_with_count_Thresholded_subset_sel$minuslog10padj < 0), ],
                   aes(color=NES), stroke=1, shape=21, fill="white")+
        geom_text(data=GSEA_result_MAX_with_count_Thresholded_subset_sel[which(GSEA_result_MAX_with_count_Thresholded_subset_sel$minuslog10padj < 0), ],
                  aes(x=minuslog10padj, y=as.numeric(DUMMY), label=Count),color="black",size=2, family="sans",fontface="bold")+
        scale_color_gradient2(name=paste("Normalized","Enrichment","Score", sep="\n"),
                              low = "blue", high = "red",mid="white",midpoint=0,
                              na.value = NA)
      
      
      Gene_set_lolliplot <-Gene_set_lolliplot+
        theme_cowplot(font_size = 2,
                      font_family = "sans")+
        facet_grid(contrast ~ identity , scales='free_x', space='free_x', switch="y", drop=TRUE)+
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
              axis.title.x=element_text(size=12,color="black", family="sans"),
              axis.text.y=element_text(size=10,color="black", family="sans", face='bold'),
              axis.text.x=element_text(size=10,color="black", family="sans"))+
        theme(legend.title = element_text(size=12),
              legend.text = element_text(size=10),
              legend.key.size = unit(0.5, 'cm'), #change legend key size
              legend.key.height = unit(0.5, 'cm'), #change legend key height
              legend.key.width = unit(0.5, 'cm'), #change legend key width
              legend.position="right")
      
      
      setwd(path_identity_sel)
      
      svgname<-paste(paste("Lolliplot",'GSEA','selected',sep='_'),".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= Gene_set_lolliplot,
               device="svg")
      }
      
      
    }# dim(GSEA_result_MAX_with_count_Thresholded_subset_sel)[1] >0
    
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
    make_option(c("--List_GSEA"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GSEA_result"), type="character", default=NULL, 
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
  
  GSEA_function(opt)
  Leading_edge_printer(opt)
  Lolliplot_and_gene_annotation(opt)
  Lolliplot_and_gene_annotation_subset(opt)
  
  
  
}


###########################################################################

system.time( main() )