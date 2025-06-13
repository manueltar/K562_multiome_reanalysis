.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis/lib/R/library"))
.libPaths()
# sessionInfo()


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
library("ggrepel")
suppressMessages(library(future))


library("optparse")


opt = NULL

options(warn = 1)


multiVals <- function(x) paste(x,collapse=";")

log_info_simple <- function(message) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  cat(timestamp, "INFO:", message, "\n")
}

volcano_plots_function = function(option_list)
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
  
  cat("DE_results_0\n")
  cat(str(DE_results))
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
  
   
   DE_results_sel<-droplevels(DE_results[which(DE_results$identity == identity_sel),])
   
   if(DEBUG == 1){
     
     cat("DE_results_sel_0\n")
     cat(str(DE_results_sel))
     cat("\n")
   }
   
   
   
   array_contrasts<-levels(DE_results_sel$contrast)
   
   for(k in 1:length(array_contrasts)){
     
     
     contrast_sel<-array_contrasts[k]
     
     cat("------------------------->\t")
     cat(sprintf(as.character(contrast_sel)))
     cat("\n")
     
     
     DE_results_sel_contrast_sel<-DE_results_sel[which(DE_results_sel$contrast == contrast_sel),]
     
     
     if(DEBUG == 1){
       
       cat("DE_results_sel_contrast_sel_0\n")
       cat(str(DE_results_sel_contrast_sel))
       cat("\n")
     }
     
     SIG_genes<-DE_results_sel_contrast_sel[which(DE_results_sel_contrast_sel$minuslog10padj >= 1.3),]
     
     
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
       
         if(dim(SIG_genes)[1] <= 50){
           
           volcano_plot<-ggplot(data=DE_results_sel_contrast_sel,
                                aes(x=log2FoldChange,
                                    y=minuslog10padj)) +
             geom_vline(xintercept=c(-0.1,0.1), color="gray", linetype='dashed',linewidth=1)+
             geom_hline(yintercept=c(1.3), color="gray", linetype='dashed',linewidth=1)+
             geom_point(data=DE_results_sel_contrast_sel[which(DE_results_sel_contrast_sel$minuslog10padj < 1.3 & DE_results_sel_contrast_sel$abslogfc < 0.1),],
                        color="gray",fill="gray", stroke=0.2, shape=21, size=2)+
             geom_point(data=DE_results_sel_contrast_sel[which(DE_results_sel_contrast_sel$minuslog10padj >= 1.3 & DE_results_sel_contrast_sel$log2FoldChange <= -0.1),],
                        color="black",fill="blue", stroke=0.2, shape=21, size=4)+
             geom_point(data=DE_results_sel_contrast_sel[which(DE_results_sel_contrast_sel$minuslog10padj >= 1.3 & DE_results_sel_contrast_sel$log2FoldChange >= 0.1),],
                        color="black",fill="red", stroke=0.2, shape=21, size=4)+
             scale_y_continuous(name="-log10 pval")+
             theme_classic()+
             ggtitle(paste(paste(paste(contrast_sel, sep=''), paste(identity_sel, sep=''),sep='  '),paste(dim(SIG_genes)[1], "DE genes", sep=" "), sep="\n"))+
             theme(plot.title=element_text(color="black", family="sans", size=14),
                   axis.title.y=element_text(color="black", family="sans", size=12),
                   axis.title.x=element_text(color="black", family="sans", size=12),
                   axis.text.y=element_text(color="black", family="sans", size=8),
                   axis.text.x=element_text(color="black", family="sans", size=8))+
             theme(legend.title = element_text(family="sans"),
                   legend.text = element_text(family="sans"),
                   legend.key.size = unit(0.35, 'cm'), #change legend key size
                   legend.key.height = unit(0.35, 'cm'), #change legend key height
                   legend.key.width = unit(0.35, 'cm'), #change legend key width
                   legend.position="hidden")+
             ggeasy::easy_center_title()
           
           volcano_plot <-volcano_plot+
             geom_text_repel(data=DE_results_sel_contrast_sel[which(DE_results_sel_contrast_sel$minuslog10padj >= 1.3 & DE_results_sel_contrast_sel$log2FoldChange <= -0.1),],
                             aes(x=log2FoldChange,
                                 y=minuslog10padj,
                                 label=gene),
                             family="sans",
                             fontface='italic',
                             color='blue',
                             segment.size  = 0.25,
                             segment.color = "blue",
                             force=25,
                             size=4,
                             box.padding = 1,
                             max.overlaps = Inf)
           
           volcano_plot <-volcano_plot+
             geom_text_repel(data=DE_results_sel_contrast_sel[which(DE_results_sel_contrast_sel$minuslog10padj >= 1.3 & DE_results_sel_contrast_sel$log2FoldChange >= 0.1),],
                             aes(x=log2FoldChange,
                                 y=minuslog10padj,
                                 label=gene),
                             family="sans",
                             fontface='italic',
                             color='red',
                             segment.size  = 0.25,
                             segment.color = "red",
                             force=25,
                             size=4,
                             box.padding = 1,
                             max.overlaps = Inf)
           
         }else{
           
           REP_genes<-SIG_genes[order(SIG_genes$abslogfc, decreasing =T),]
           
           if(DEBUG == 1){
             
             cat("REP_genes_0\n")
             cat(str(REP_genes))
             cat("\n")
           }
           
           
           top25_genes<-REP_genes[c(1:25),]
           
           if(DEBUG == 1){
             
             cat("top25_genes_0\n")
             cat(str(top25_genes))
             cat("\n")
           }
           
           
           
           volcano_plot<-ggplot(data=DE_results_sel_contrast_sel,
                                aes(x=log2FoldChange,
                                    y=minuslog10padj)) +
             geom_vline(xintercept=c(-0.1,0.1), color="gray", linetype='dashed',linewidth=1)+
             geom_hline(yintercept=c(1.3), color="gray", linetype='dashed',linewidth=1)+
             geom_point(data=DE_results_sel_contrast_sel[which(DE_results_sel_contrast_sel$minuslog10padj < 1.3 & DE_results_sel_contrast_sel$abslogfc < 0.1),],
                        color="gray",fill="gray", stroke=0.2, shape=21, size=2)+
             geom_point(data=DE_results_sel_contrast_sel[which(DE_results_sel_contrast_sel$minuslog10padj >= 1.3 & DE_results_sel_contrast_sel$log2FoldChange <= -0.1),],
                        color="black",fill="blue", stroke=0.2, shape=21, size=4)+
             geom_point(data=DE_results_sel_contrast_sel[which(DE_results_sel_contrast_sel$minuslog10padj >= 1.3 & DE_results_sel_contrast_sel$log2FoldChange >= 0.1),],
                        color="black",fill="red", stroke=0.2, shape=21, size=4)+
             scale_y_continuous(name="-log10 pval")+
             theme_classic()+
             ggtitle(paste(paste(paste(contrast_sel, sep=''), paste(identity_sel, sep=''),sep='  '),paste(dim(SIG_genes)[1], "DE genes", sep=" "), sep="\n"))+
             theme(plot.title=element_text(color="black", family="sans", size=14),
                   axis.title.y=element_text(color="black", family="sans", size=12),
                   axis.title.x=element_text(color="black", family="sans", size=12),
                   axis.text.y=element_text(color="black", family="sans", size=8),
                   axis.text.x=element_text(color="black", family="sans", size=8))+
             theme(legend.title = element_text(family="sans"),
                   legend.text = element_text(family="sans"),
                   legend.key.size = unit(0.35, 'cm'), #change legend key size
                   legend.key.height = unit(0.35, 'cm'), #change legend key height
                   legend.key.width = unit(0.35, 'cm'), #change legend key width
                   legend.position="hidden")+
             ggeasy::easy_center_title()
           
           volcano_plot <-volcano_plot+
             geom_text_repel(data=top25_genes[which(top25_genes$log2FoldChange <= -0.1),],
                             aes(x=log2FoldChange,
                                 y=minuslog10padj,
                                 label=gene),
                             family="sans",
                             fontface='italic',
                             color='blue',
                             segment.size  = 0.25,
                             segment.color = "blue",
                             force=25,
                             size=4,
                             box.padding = 1,
                             max.overlaps = Inf)
           
           volcano_plot <-volcano_plot+
             geom_text_repel(data=top25_genes[which(top25_genes$log2FoldChange >= 0.1),],
                             aes(x=log2FoldChange,
                                 y=minuslog10padj,
                                 label=gene),
                             family="sans",
                             fontface='italic',
                             color='red',
                             segment.size  = 0.25,
                             segment.color = "red",
                             force=25,
                             size=4,
                             box.padding = 1,
                             max.overlaps = Inf)
           
           
         }#dim(SIG_genes)[1] <= 50
         
         
         
        
       setwd(path_contrast_sel)
       
       svgname<-paste("volcano_plot_",Diff_sel,".svg",sep='')
       
       ggsave(svgname,plot=volcano_plot, device ='svg', height=5, width=5)
         

       
       
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
  
  volcano_plots_function(opt)
  
  
  
  
}


###########################################################################

system.time( main() )