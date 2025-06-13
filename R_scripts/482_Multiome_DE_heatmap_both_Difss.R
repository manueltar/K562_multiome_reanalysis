.libPaths()
assign(".lib.loc", "/home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis/lib/R/library", envir = environment(.libPaths))
.libPaths()
# # sessionInfo()


# .libPaths()
# .libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis/lib/R/library"))
# .libPaths()
# sessionInfo()

library("optparse")
suppressMessages(library(dplyr)) 
suppressMessages(library(ggplot2)) 
suppressMessages(library(Matrix)) 
suppressMessages(library(data.table)) 
suppressMessages(library(ggpubr)) 
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
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
suppressMessages(library("tibble")) 
suppressMessages(library(future))

opt = NULL

options(warn = 1)


multiVals <- function(x) paste(x,collapse=";")

log_info_simple <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(timestamp, "INFO:", message, "\n")
}

heatmap_function = function(option_list)
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
  

  #### READ and transform DE_results ----
  
  
  DE_results<-readRDS(file=opt$DE_results)
  
  cat("DE_result_0\n")
  cat(str(DE_results))
  cat("\n")
  
  #### READ and transform normalised_counts ----
  
  
  normalised_counts<-readRDS(file=opt$normalised_counts)
  
  cat("normalised_counts_0\n")
  cat(str(normalised_counts))
  cat("\n")
  
  #### LOOP of identities  -----
  
  
  array_identities<-levels(DE_results$identity)
  
  
  cat("array_identities_0\n")
  cat(str(array_identities))
  cat("\n")
  
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
  
   
   DE_result_sel<-droplevels(DE_results[which(DE_results$identity == identity_sel),])
   
   if(DEBUG == 1){
     
     cat("DE_result_sel_0\n")
     cat(str(DE_result_sel))
     cat("\n")
   }
   
   
   
   array_contrasts<-levels(DE_result_sel$contrast)
   
   for(k in 1:length(array_contrasts)){
     
     
     contrast_sel<-array_contrasts[k]
     
     cat("------------------------->\t")
     cat(sprintf(as.character(contrast_sel)))
     cat("\n")
     
     
     DE_result_sel_contrast_sel<-DE_result_sel[which(DE_result_sel$contrast == contrast_sel),]
     
     
     if(DEBUG == 1){
       
       cat("DE_result_sel_contrast_sel_0\n")
       cat(str(DE_result_sel_contrast_sel))
       cat("\n")
     }
     
     SIG_genes<-DE_result_sel_contrast_sel[which(DE_result_sel_contrast_sel$minuslog10padj >= 1.3),]
     
     
     if(DEBUG == 1){
       
       cat("SIG_genes_0\n")
       cat(str(SIG_genes))
       cat("\n")
     }
     
     
     if(dim(SIG_genes)[1] >0){
       
       
       
       
       path_contrast_sel<-paste(path_identity_sel,contrast_sel,'/',sep='')
       
       if (file.exists(path_contrast_sel)){
         
       }else{
         
         dir.create(path_contrast_sel)
       }
       
       
       REP<-normalised_counts[which(normalised_counts$gene%in%SIG_genes$gene &
                                      normalised_counts$identity == identity_sel),]
       
       if(DEBUG == 1){
         
         cat("REP_0\n")
         cat(str(REP))
         cat("\n")
       }
       
       
       REP_wide<-as.data.frame(pivot_wider(REP, id_cols=c('gene'), 
                                           names_from=c('identity',"clone_line","time_point"), 
                                           values_from='count',
                                           names_sep='|'), stringsAsFactors=F)
       
       if(DEBUG == 1){
         
         cat("REP_wide_0\n")
         cat(str(REP_wide))
         cat("\n")
       }
       
       GeneEXP_matrix<-as.matrix(REP_wide[,-which(colnames(REP_wide)%in%c('gene'))])
       
       row.names(GeneEXP_matrix)<-REP_wide$gene
       
       if(DEBUG == 1){
         
         cat("GeneEXP_matrix_0\n")
         cat(str(GeneEXP_matrix))
         cat("\n")
       }           
       
       annotation_col<- data.frame(matrix(vector(), length(colnames(GeneEXP_matrix)), 3,
                                          dimnames=list(c(),
                                                        c("identity","time_point","clone_line"))),stringsAsFactors=F)
       
       row.names(annotation_col)<-colnames(GeneEXP_matrix)
       

       
       if(DEBUG == 1){
         cat("annotation_col_0\n")
         cat(str(annotation_col))
         cat("\n")
         cat(str(row.names(annotation_col)))
         cat("\n")
         
       }
       
       
       annotation_col$identity<-gsub("\\|.+$","", row.names(annotation_col))
       annotation_col$clone_line<-gsub("^[^\\|]+\\|","", row.names(annotation_col))
       annotation_col$clone_line<-gsub("\\|.+$","", annotation_col$clone_line)
       
       annotation_col$time_point<-gsub("^[^\\|]+\\|[^\\|]+\\|","", row.names(annotation_col))
       
       if(DEBUG == 1){
         cat("annotation_col_1\n")
         cat(str(annotation_col))
         cat("\n")
       }
       
       annotation_col$Genotype<-NA
       
       annotation_col$Genotype[which(annotation_col$clone_line%in%c('wt_1','wt_2','wt_3'))]<-'wt'
       annotation_col$Genotype[which(annotation_col$clone_line%in%c('rs139141690_HET_1'))]<-'rs139141690_HET'
       annotation_col$Genotype[which(annotation_col$clone_line%in%c('rs139141690_1','rs139141690_2','rs139141690_3'))]<-'rs139141690'
       annotation_col$Genotype[which(annotation_col$clone_line%in%c('Del_16bp_1'))]<-'Del_16bp'
       annotation_col$Genotype[which(annotation_col$clone_line%in%c('Del_80bp_1','Del_80bp_2','Del_80bp_3'))]<-'Del_80bp'
       
       if(DEBUG == 1){
         cat("annotation_col_PRE\n")
         cat(str(annotation_col))
         cat("\n")
         names(summary(as.factor(annotation_col$Genotype)))
         cat("\n")
         names(summary(as.factor(annotation_col$clone_line)))
         cat("\n")
         names(summary(as.factor(annotation_col$identity)))
         cat("\n")
       }
       
       annotation_col$Genotype<-factor(annotation_col$Genotype,
                                       levels=c('wt','rs139141690_HET','rs139141690','Del_16bp','Del_80bp'),
                                       ordered=T)
       
       annotation_col$clone_line<-factor(annotation_col$clone_line,
                                         levels=levels(normalised_counts$clone_line),
                                         ordered=T)
       
       annotation_col$identity<-factor(annotation_col$identity,
                                       levels=array_identities,
                                       ordered=T)
       
       annotation_col$time_point<-factor(annotation_col$time_point,
                                         levels=levels(normalised_counts$time_point),
                                         ordered=T)
       
       if(DEBUG == 1){
         cat("annotation_col_POST\n")
         cat(str(annotation_col))
         cat("\n")
         names(summary(as.factor(annotation_col$Genotype)))
         cat("\n")
         names(summary(as.factor(annotation_col$clone_line)))
         cat("\n")
         names(summary(as.factor(annotation_col$identity)))
         cat("\n")
       }
       
       
       vector_colors_clone_line<- c(brewer.pal(9, "Greens")[c(5,6,7)],
                                                  brewer.pal(9, "YlOrRd")[c(2)],
                                                  brewer.pal(9, "Reds")[c(5,6,7)],
                                                  brewer.pal(9, "Purples")[c(7)],
                                                  brewer.pal(9, "Blues")[c(4,5,6)])

        names(vector_colors_clone_line)<-levels(annotation_col$clone_line)
        
        vector_colors_Genotype<-c(brewer.pal(9, "Greens")[c(5)],
                                  brewer.pal(9, "YlOrRd")[c(2)],
                                    brewer.pal(9, "Reds")[c(5)],
                                    brewer.pal(9, "Purples")[c(7)],
                                    brewer.pal(9, "Blues")[c(4)])
        
        names(vector_colors_Genotype)<-levels(annotation_col$Genotype)
        
        vector_colors_time_point<-c(brewer.pal(9, "Greys")[c(1)],
                                    brewer.pal(9, "Greys")[c(3)],
                                    brewer.pal(9, "Greys")[c(5)],
                                    brewer.pal(9, "Greys")[c(7)])
        
        names(vector_colors_time_point)<-levels(annotation_col$time_point)
        
        
        
        
        vector_colors_seurat_clusters<-c(brewer.pal(9, "YlOrRd")[c(7)],
                                         brewer.pal(9, "Blues")[c(5)],
                                         brewer.pal(9, "RdPu")[c(6)],
                                         brewer.pal(9, "Blues")[c(4)],
                                         brewer.pal(9, "Greens")[c(6)],
                                         brewer.pal(9, "YlOrRd")[c(6)],
                                         brewer.pal(9, "Blues")[c(3)],
                                         brewer.pal(9, "RdPu")[c(5)],
                                         brewer.pal(9, "Greens")[c(5)],
                                         brewer.pal(9, "Greens")[c(4)],
                                         brewer.pal(9, "Blues")[c(2)],
                                         brewer.pal(9, "RdPu")[c(4)],
                                         brewer.pal(9, "RdPu")[c(3)])
        
        values<-vector_colors_seurat_clusters
        

        # harcoded!!!!
        
        names<-array_identities
        
        values<-values[1:length(names)]
        
        
        vector_colors_identity<-setNames(values, names)
        
        if(DEBUG == 1){
          cat("vector_colors_identity_0\n")
          cat(str(vector_colors_identity))
          cat("\n")
        }
        
        vector_colors_identity<-vector_colors_identity[which(names(vector_colors_identity)%in%levels(annotation_col$identity))]

                                  
        
        if(DEBUG == 1){
          cat("vector_colors_identity_1\n")
          cat(str(vector_colors_identity))
          cat("\n")
        }
        
        
        ## ChatGPT check missing identitites
        
        missing_identities <- setdiff(levels(annotation_col$identity), names(vector_colors_identity))
        
        if (length(missing_identities) > 0) {
          cat("Missing colors for the following identities:\n")
          print(missing_identities)
          stop("Please add these identities to the color vector.")
        }
        
       
       ann_colors <- list( clone_line = vector_colors_clone_line,
                           Genotype = vector_colors_Genotype,
                           time_point = vector_colors_time_point,
                           identity =vector_colors_identity)
      
       
       if(DEBUG == 1){
         cat("ann_colors_0\n")
         cat(str(ann_colors))
         cat("\n")
       }
       
       
       FLAG_log_pval<-NA
       
       
         if(dim(SIG_genes)[1] <= 50 & dim(SIG_genes)[1] > 1){
           
           
           heatmap<-pheatmap(GeneEXP_matrix, display_numbers = FALSE,
                             show_colnames=FALSE,
                             angle_col = "0",
                             clustering_method="ward.D2",
                             fontsize_row = 8, 
                             fontsize_col = 8,
                             breaks=seq(-2,2,length.out=101),
                             color=colorRampPalette(c("blue","white","red"))(100),
                             scale="row",
                             cluster_cols=FALSE,
                             border_color='black',
                             treeheight_row=70, treeheight_col=70, cutree_cols=7,
                             annotation_col = annotation_col,
                             annotation_colors = ann_colors)
           
           setwd(path_contrast_sel)
           
           svgname<-paste("Heatmap_DE_genes_",Diff_sel,".svg",sep='')
           
           ggsave(svgname,plot=heatmap, device ='svg')
           
           FLAG_log_pval<-1
         
           
         }else{
           
           if(dim(SIG_genes)[1] > 50){
             
             heatmap<-pheatmap(GeneEXP_matrix, display_numbers = FALSE,
                               show_colnames=FALSE,
                               show_rownames=FALSE,
                               angle_col = "0",
                               clustering_method="ward.D2",
                               fontsize_row = 8, 
                               fontsize_col = 8,
                               breaks=seq(-2,2,length.out=101),
                               color=colorRampPalette(c("blue","white","red"))(100),
                               scale="row",
                               cluster_cols=FALSE,
                               border_color='black',
                               treeheight_row=70, treeheight_col=70, cutree_cols=7,
                               annotation_col = annotation_col,
                               annotation_colors = ann_colors)
             
             setwd(path_contrast_sel)
             
             svgname<-paste("Heatmap_DE_genes_",Diff_sel,".svg",sep='')
             
             ggsave(svgname,plot=heatmap, device ='svg')
             
             FLAG_log_pval<-1
             
           }else{
             
             cat("no heatmap\n")
             
             FLAG_log_pval<-0
             
           }# dim(SIG_genes)[1] > 50
           
         }#dim(SIG_genes)[1] <= 50
         
         
         
        
      
       
       
       
       ##### logpval plot ----
       
       if(FLAG_log_pval == 1){
         selected_genes_after_heatmap_clustering<-heatmap$tree_row$labels[heatmap$tree_row$order]
         
         if(DEBUG == 1){
           cat("selected_genes_after_heatmap_clustering_0\n")
           cat(str(selected_genes_after_heatmap_clustering))
           cat("\n")
         }
         
         logpval_df<-DE_result_sel[which(DE_result_sel$gene%in%SIG_genes$gene),]
         
         
         if(DEBUG == 1){
           cat("logpval_df_0\n")
           cat(str(logpval_df))
           cat("\n")
         }
         
         logpval_df$gene<-factor(logpval_df$gene, levels=rev(selected_genes_after_heatmap_clustering), ordered=T)
         
         
         
         logpval_df$SIG<-NA
         
         logpval_df$SIG[which(logpval_df$minuslog10padj >= 1.3)]<-'YES'
         logpval_df$SIG[which(logpval_df$minuslog10padj < 1.3)]<-'NO'
         
         
         logpval_df$SIG<-factor(logpval_df$SIG, levels=c('NO','YES'), ordered=T)
         
         if(DEBUG == 1){
           cat("logpval_df_1\n")
           cat(str(logpval_df))
           cat("\n")
         }
         
         
         vector_fill<-c(brewer.pal(9, "YlOrRd")[c(2)],
                        brewer.pal(9, "Reds")[c(5)],
                        brewer.pal(9, "Purples")[c(7)],
                        brewer.pal(9, "Blues")[c(4)],
                        "steelblue")
         

         
         logpval_dotplot<-ggplot(data=logpval_df,
                                 aes(y=gene,
                                     x=contrast))+
           geom_point(aes(size=minuslog10padj, 
                          color=SIG,
                          fill=log2FoldChange),
                      stroke=1, shape=21)+
           scale_size(range = c(0,6), name='-log10pval')+
           scale_y_discrete(name=NULL)+
           scale_x_discrete(name=NULL)+
           scale_fill_gradient2(
             low = "blue", 
             mid = "white", 
             high = "red", 
             midpoint = 0)+
           scale_color_manual(name='p < 0.05',values=c('gray','black'))+
           theme_classic()+
           theme(axis.title.y=element_blank(),
                 axis.title.x=element_blank(),
                 axis.text.y=element_text(size=8, color="black", family="sans"),
                 axis.text.x=element_text(size=4, color="black", family="sans"),
                 axis.line.x = element_line(size = 0.4),
                 axis.ticks.x = element_blank(),
                 axis.ticks.y = element_line(size = 0.4),
                 axis.line.y = element_line(size = 0.4))+
           theme(legend.title = element_text(size=12),
                 legend.text = element_text(size=8),
                 legend.key.size = unit(0.5, 'cm'), #change legend key size
                 legend.key.height = unit(0.5, 'cm'), #change legend key height
                 legend.key.width = unit(0.5, 'cm'), #change legend key width
                 legend.position="right")+
           ggeasy::easy_center_title()
         
         setwd(path_contrast_sel)
         
         svgname<-paste("logpval_dotplot_DE_genes_",Diff_sel,".svg",sep='')
         
         ggsave(svgname,plot=logpval_dotplot, device ='svg')
         
       }#FLAG_log_pval == 1
       
     
       
       
     }#dim(SIG_genes)[1] >0
   }#k in 1:length(array_contrasts
}# i in 1:length(array_identities)
  
  
  
  
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
    make_option(c("--DE_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--normalised_counts"), type="character", default=NULL, 
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
  
  heatmap_function(opt)
  
  
  
  
}


###########################################################################

system.time( main() )