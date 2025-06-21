
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ActivePathways", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dorothea", lib="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("decoupleR", lib="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rWikiPathways", lib="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

create_gmt = function(option_list)
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
  
  #### new_sets_Manual_curation ----
  
  new_sets<-read.table(opt$new_sets, sep=",", header=T)
  
  
  cat("new_sets_0\n")
  cat(str(new_sets))
  cat("\n")
  cat(str(unique(new_sets$GeneSet)))
  cat("\n")
  
  DEBUG<-1
  
  
  adapted_df<-data.frame()
  
  for(i in 1:dim(new_sets)[2]){
    
    new_sets_sel<-new_sets[,i]
    
    if(DEBUG == 1){
      
      cat("new_sets_sel_0\n")
      cat(str(new_sets_sel))
      cat("\n")
      
    }#DEBUG == 1
    
    FLAG_empty<-which(new_sets_sel == "")
    
    if(length(FLAG_empty) >0){
      
      new_sets_sel<-new_sets_sel[-which(new_sets_sel == "")]
      
    }else{
      
      new_sets_sel<-new_sets_sel
      
    }#length(FLAG_empty) >0
    
    if(DEBUG == 1){
      
      cat("new_sets_sel_1\n")
      cat(str(new_sets_sel))
      cat("\n")
      
    }#DEBUG == 1
    
    new_sets_sel_collapsed<-paste(new_sets_sel, collapse=",")
    
    if(DEBUG == 1){
      
      cat("new_sets_sel_collapsed_1\n")
      cat(str(new_sets_sel_collapsed))
      cat("\n")
      
    }#DEBUG == 1
    
    GeneSet<-colnames(new_sets)[i]
    
    if(DEBUG == 1){
      
      cat("GeneSet\n")
      cat(sprintf(as.character((GeneSet))))
      cat("\n")
      
    }#DEBUG == 1
    
    a_df<-as.data.frame(cbind(GeneSet,new_sets_sel_collapsed), stringsAsFactors=F)
    
    colnames(a_df)<-c("GeneSet","Genes")
    
    if(DEBUG == 1){
      
      cat("a_df_0\n")
      cat(str(a_df))
      cat("\n")
      
    }#DEBUG == 1
    
    adapted_df<-rbind(a_df,adapted_df)
    

    
   
    
  }#i in 1:dim(new_sets)[1]
  
  cat("adapted_df_0\n")
  cat(str(adapted_df))
  cat("\n")
  
  new_sets_long<-unique(as.data.frame(cSplit(adapted_df,sep = ',', direction = "long",
                                                       splitCols = "Genes"),stringsAsFactors=F))
  
  cat("new_sets_long_0\n")
  cat(str(new_sets_long))
  cat("\n")
  cat(str(unique(new_sets_long$GeneSet)))
  cat("\n")
  
  new_sets_long$ENTREZID<-mapIds(org.Hs.eg.db, keys=new_sets_long$Genes, keytype="SYMBOL",column="ENTREZID")
  new_sets_long$ensembl_gene_id<-mapIds(org.Hs.eg.db, keys=new_sets_long$Genes, keytype="SYMBOL",column="ENSEMBL")
  
  cat("new_sets_long_1\n")
  cat(str(new_sets_long))
  cat("\n")
  
  new_sets_long_NO_NA<-new_sets_long[!is.na(new_sets_long$ENTREZID),]
  
  cat("new_sets_long_NO_NA_0\n")
  cat(str(new_sets_long_NO_NA))
  cat("\n")
  
  #### Prepare the file for gmt export ----
  
  for_gmt_new_sets_long_NO_NA<-unique(new_sets_long_NO_NA[,c(which(colnames(new_sets_long_NO_NA) == 'GeneSet'),which(colnames(new_sets_long_NO_NA) == 'ENTREZID'))])
  
  colnames(for_gmt_new_sets_long_NO_NA)[which(colnames(for_gmt_new_sets_long_NO_NA) == 'GeneSet')]<-'id'
  colnames(for_gmt_new_sets_long_NO_NA)[which(colnames(for_gmt_new_sets_long_NO_NA) == 'ENTREZID')]<-'gene'
  
  cat("for_gmt_new_sets_long_NO_NA_0\n")
  cat(str(for_gmt_new_sets_long_NO_NA))
  cat("\n")
  
  
  for_gmt_new_sets_long_NO_NA$name<-for_gmt_new_sets_long_NO_NA$id
  
  
  
  cat("for_gmt_new_sets_long_NO_NA_1\n")
  cat(str(for_gmt_new_sets_long_NO_NA))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(for_gmt_new_sets_long_NO_NA$name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(for_gmt_new_sets_long_NO_NA$name)))))
  cat("\n")
  
  #### Reorder
  
  indx.reorder<-c(which(colnames(for_gmt_new_sets_long_NO_NA) == 'id'),which(colnames(for_gmt_new_sets_long_NO_NA) == 'name'),which(colnames(for_gmt_new_sets_long_NO_NA) == 'gene'))
  
  
  for_gmt_new_sets_long_NO_NA_reordered<-for_gmt_new_sets_long_NO_NA[,indx.reorder]
  
  cat("for_gmt_new_sets_long_NO_NA_reordered_0\n")
  cat(str(for_gmt_new_sets_long_NO_NA_reordered))
  cat("\n")
  
  
  #### Table_of_gene_sets_Manual_curation ----
  
  Table_of_gene_sets<-read.table(opt$Table_of_gene_sets, sep="\t", header=T)
  

  cat("Table_of_gene_sets_0\n")
  cat(str(Table_of_gene_sets))
  cat("\n")
  cat(str(unique(Table_of_gene_sets$GeneSet)))
  cat("\n")
  
  Table_of_gene_sets_long<-unique(as.data.frame(cSplit(Table_of_gene_sets,sep = ',', direction = "long",
                                               splitCols = "Genes"),stringsAsFactors=F))
  
  cat("Table_of_gene_sets_long_0\n")
  cat(str(Table_of_gene_sets_long))
  cat("\n")
  cat(str(unique(Table_of_gene_sets_long$GeneSet)))
  cat("\n")
  
  Table_of_gene_sets_long$ENTREZID<-mapIds(org.Hs.eg.db, keys=Table_of_gene_sets_long$Genes, keytype="SYMBOL",column="ENTREZID")
  Table_of_gene_sets_long$ensembl_gene_id<-mapIds(org.Hs.eg.db, keys=Table_of_gene_sets_long$Genes, keytype="SYMBOL",column="ENSEMBL")
  
  cat("Table_of_gene_sets_long_1\n")
  cat(str(Table_of_gene_sets_long))
  cat("\n")
  
  Table_of_gene_sets_long_NO_NA<-Table_of_gene_sets_long[!is.na(Table_of_gene_sets_long$ENTREZID),]
  
  cat("Table_of_gene_sets_long_NO_NA_0\n")
  cat(str(Table_of_gene_sets_long_NO_NA))
  cat("\n")
  
  #### Prepare the file for gmt export ----
  
  for_gmt_Table_of_gene_sets_long_NO_NA<-unique(Table_of_gene_sets_long_NO_NA[,c(which(colnames(Table_of_gene_sets_long_NO_NA) == 'GeneSet'),which(colnames(Table_of_gene_sets_long_NO_NA) == 'ENTREZID'))])
  
  colnames(for_gmt_Table_of_gene_sets_long_NO_NA)[which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'GeneSet')]<-'id'
  colnames(for_gmt_Table_of_gene_sets_long_NO_NA)[which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'ENTREZID')]<-'gene'
  
  cat("for_gmt_Table_of_gene_sets_long_NO_NA_0\n")
  cat(str(for_gmt_Table_of_gene_sets_long_NO_NA))
  cat("\n")
  
  
  for_gmt_Table_of_gene_sets_long_NO_NA$name<-for_gmt_Table_of_gene_sets_long_NO_NA$id
  
  
  
  cat("for_gmt_Table_of_gene_sets_long_NO_NA_1\n")
  cat(str(for_gmt_Table_of_gene_sets_long_NO_NA))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(for_gmt_Table_of_gene_sets_long_NO_NA$name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(for_gmt_Table_of_gene_sets_long_NO_NA$name)))))
  cat("\n")
  
  #### Reorder
  
  indx.reorder<-c(which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'id'),which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'name'),which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'gene'))
  
  
  for_gmt_Table_of_gene_sets_long_NO_NA_reordered<-for_gmt_Table_of_gene_sets_long_NO_NA[,indx.reorder]
  
  cat("for_gmt_Table_of_gene_sets_long_NO_NA_reordered_0\n")
  cat(str(for_gmt_Table_of_gene_sets_long_NO_NA_reordered))
  cat("\n")
  
  #### miRNA_table_Manual_curation ----
  
  miRNA_table<-read.table(opt$miRNA_table, sep=",", header=T)
  
  
  cat("miRNA_table_0\n")
  cat(str(miRNA_table))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(miRNA_table$Evidence))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(miRNA_table$Evidence)))))
  cat("\n")
  
  miRNA_table$ENTREZID<-mapIds(org.Hs.eg.db, keys=miRNA_table$Target, keytype="SYMBOL",column="ENTREZID")
  miRNA_table$ensembl_gene_id<-mapIds(org.Hs.eg.db, keys=miRNA_table$Target, keytype="SYMBOL",column="ENSEMBL")
  
  cat("miRNA_table_1\n")
  cat(str(miRNA_table))
  cat("\n")
  
  miRNA_table_NO_NA<-miRNA_table[!is.na(miRNA_table$ENTREZID),]
  
  cat("miRNA_table_NO_NA_0\n")
  cat(str(miRNA_table_NO_NA))
  cat("\n")
  
  indx.experimental<-which(miRNA_table_NO_NA$Evidence%in%c('experimental (all)',
                                                           'experimental (all) + predicted (intersection) + predicted (union)',
                                                           'experimental (all) + predicted (union)'))
  
  
  
  
  
  miRNA_table_experimental<- data.frame(matrix(vector(), length(indx.experimental), 2,
                             dimnames=list(c(),
                                           c("GeneSet","ENTREZID"))),stringsAsFactors=F)
  
  
  
  miRNA_table_experimental$ENTREZID<-miRNA_table_NO_NA$ENTREZID[indx.experimental]
  miRNA_table_experimental$GeneSet<-rep("mir5739_experimental_targets",length(indx.experimental))                                    
  

  indx.intersection<-c(indx.experimental, which(miRNA_table_NO_NA$Evidence%in%c('predicted (intersection) + predicted (union)')))
  
  
  
  miRNA_table_intersection<- data.frame(matrix(vector(), length(indx.intersection), 2,
                                               dimnames=list(c(),
                                                             c("GeneSet","ENTREZID"))),stringsAsFactors=F)
  
  
  
  miRNA_table_intersection$ENTREZID<-miRNA_table_NO_NA$ENTREZID[indx.intersection]
  miRNA_table_intersection$GeneSet<-rep("mir5739_experimental_and_intersection_targets",length(indx.intersection))
  

  
  indx.union<-c(indx.experimental,indx.intersection, which(miRNA_table_NO_NA$Evidence%in%c('predicted (union)')))
  
  
  
  miRNA_table_union<- data.frame(matrix(vector(), length(indx.union), 2,
                                               dimnames=list(c(),
                                                             c("GeneSet","ENTREZID"))),stringsAsFactors=F)
  
  
  
  miRNA_table_union$ENTREZID<-miRNA_table_NO_NA$ENTREZID[indx.union]
  miRNA_table_union$GeneSet<-rep("mir5739_experimental_intersection_and_union_targets",length(indx.union))            
  
  
  #### Prepare the file for gmt export ----
  
  for_gmt_mir5739_targets<-rbind(miRNA_table_experimental,miRNA_table_intersection,miRNA_table_union)
  
  colnames(for_gmt_mir5739_targets)[which(colnames(for_gmt_mir5739_targets) == 'GeneSet')]<-'id'
  colnames(for_gmt_mir5739_targets)[which(colnames(for_gmt_mir5739_targets) == 'ENTREZID')]<-'gene'
  
  cat("for_gmt_mir5739_targets_0\n")
  cat(str(for_gmt_mir5739_targets))
  cat("\n")
  
  
  for_gmt_mir5739_targets$name<-for_gmt_mir5739_targets$id
  
  
  
  cat("for_gmt_mir5739_targets_1\n")
  cat(str(for_gmt_mir5739_targets))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(for_gmt_mir5739_targets$name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(for_gmt_mir5739_targets$name)))))
  cat("\n")
  
  
  #### Reorder
  
  indx.reorder<-c(which(colnames(for_gmt_mir5739_targets) == 'id'),which(colnames(for_gmt_mir5739_targets) == 'name'),which(colnames(for_gmt_mir5739_targets) == 'gene'))
  
  
  for_gmt_mir5739_targets_reordered<-for_gmt_mir5739_targets[,indx.reorder]
  
  cat("for_gmt_mir5739_targets_reordered_0\n")
  cat(str(for_gmt_mir5739_targets_reordered))
  cat("\n")
  
  
  
  #### CUX1_TARGET_GENES ----
  
  c3_collection<-readRDS(opt$c3_collection)
  
  
  cat("c3_collection_0\n")
  cat(str(c3_collection))
  cat("\n")
  
  
  CUX1_TARGET_GENES_df<-c3_collection[which(c3_collection$term == 'CUX1_TARGET_GENES'),]
  
  
  cat("CUX1_TARGET_GENES_df_0\n")
  cat(str(CUX1_TARGET_GENES_df))
  cat("\n")
  
  #### Prepare the file for gmt export ----
  
  for_gmt_CUX1_TARGET_GENES_df<-unique(CUX1_TARGET_GENES_df[,c(which(colnames(CUX1_TARGET_GENES_df) == 'term'),which(colnames(CUX1_TARGET_GENES_df) == 'gene'))])
  
  colnames(for_gmt_CUX1_TARGET_GENES_df)[which(colnames(for_gmt_CUX1_TARGET_GENES_df) == 'term')]<-'id'

  cat("for_gmt_CUX1_TARGET_GENES_df_0\n")
  cat(str(for_gmt_CUX1_TARGET_GENES_df))
  cat("\n")
  
  
  for_gmt_CUX1_TARGET_GENES_df$name<-for_gmt_CUX1_TARGET_GENES_df$id
  
  
  
  cat("for_gmt_CUX1_TARGET_GENES_df_1\n")
  cat(str(for_gmt_CUX1_TARGET_GENES_df))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(for_gmt_CUX1_TARGET_GENES_df$name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(for_gmt_CUX1_TARGET_GENES_df$name)))))
  cat("\n")
  
  #### Reorder
  
  indx.reorder<-c(which(colnames(for_gmt_CUX1_TARGET_GENES_df) == 'id'),which(colnames(for_gmt_CUX1_TARGET_GENES_df) == 'name'),which(colnames(for_gmt_CUX1_TARGET_GENES_df) == 'gene'))
  
  
  for_gmt_CUX1_TARGET_GENES_df_reordered<-for_gmt_CUX1_TARGET_GENES_df[,indx.reorder]
  
  cat("for_gmt_CUX1_TARGET_GENES_df_reordered_0\n")
  cat(str(for_gmt_CUX1_TARGET_GENES_df_reordered))
  cat("\n")
  #############################################################
  
  #### RUNX1_TARGET_GENES ----
  
  Dorothea_ABC<-readRDS(opt$Dorothea_ABC)
  
  
  cat("Dorothea_ABC_0\n")
  cat(str(Dorothea_ABC))
  cat("\n")
  
  
  RUNX1_TARGET_GENES_df<-Dorothea_ABC[which(Dorothea_ABC$term == 'Dorothea_ABC_RUNX1_targets'),]
  
  
  cat("RUNX1_TARGET_GENES_df_0\n")
  cat(str(RUNX1_TARGET_GENES_df))
  cat("\n")
  
  #### RUNX1_TARGET_GENES ----
  
  Dorothea_ABCD<-readRDS(opt$Dorothea_ABCD)
  
  
  cat("Dorothea_ABCD_0\n")
  cat(str(Dorothea_ABCD))
  cat("\n")
  
  
  CUX1_TARGET_GENES_df<-Dorothea_ABCD[which(Dorothea_ABCD$term == 'Dorothea_ABCD_CUX1_targets'),]
  
  
  cat("CUX1_TARGET_GENES_df_0\n")
  cat(str(CUX1_TARGET_GENES_df))
  cat("\n")
  
  
  ### merge and Csplitlong---
  
  doro<-rbind(RUNX1_TARGET_GENES_df,
              CUX1_TARGET_GENES_df)
  
  cat("doro_0\n")
  cat(str(doro))
  cat("\n")
  
  
  doro.dt<-data.table(doro, key="gene")
  
  
  doro_collapsed<-as.data.frame(doro.dt[,.(term_string=paste(unique(term), collapse=" & ")), by=key(doro.dt)], stringsAsFactors=F)
  
  cat("doro_collapsed_0\n")
  cat(str(doro_collapsed))
  cat("\n")
  
  doro_collapsed$term_string[which(doro_collapsed$term_string == 'Dorothea_ABCD_CUX1_targets')]<-'Exclusive_Dorothea_ABCD_CUX1_targets'
  doro_collapsed$term_string[which(doro_collapsed$term_string == 'Dorothea_ABC_RUNX1_targets')]<-'Exclusive_Dorothea_ABC_RUNX1_targets'
  
  
  colnames(doro_collapsed)[which(colnames(doro_collapsed) == 'term_string')]<-'term'
                                        
  cat("doro_collapsed_1\n")
  cat(str(doro_collapsed))
  cat("\n")                                  
  
  
  
  #### Prepare the file for gmt export ----
  
  for_gmt_doro_collapsed<-unique(doro_collapsed[,c(which(colnames(doro_collapsed) == 'term'),which(colnames(doro_collapsed) == 'gene'))])
  
  colnames(for_gmt_doro_collapsed)[which(colnames(for_gmt_doro_collapsed) == 'term')]<-'id'
  
  cat("for_gmt_doro_collapsed_0\n")
  cat(str(for_gmt_doro_collapsed))
  cat("\n")
  
  
  for_gmt_doro_collapsed$name<-for_gmt_doro_collapsed$id
  
  
  
  cat("for_gmt_doro_collapsed_1\n")
  cat(str(for_gmt_doro_collapsed))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(for_gmt_doro_collapsed$name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(for_gmt_doro_collapsed$name)))))
  cat("\n")
  
  #### Reorder
  
  indx.reorder<-c(which(colnames(for_gmt_doro_collapsed) == 'id'),which(colnames(for_gmt_doro_collapsed) == 'name'),which(colnames(for_gmt_doro_collapsed) == 'gene'))
  
  
  for_gmt_doro_collapsed_reordered<-for_gmt_doro_collapsed[,indx.reorder]
  
  cat("for_gmt_doro_collapsed_reordered_0\n")
  cat(str(for_gmt_doro_collapsed_reordered))
  cat("\n")
  #############################################################
  
  #### SIMBA_GRN_targets ----
  
  SIMBA_GRN_targets<-read.table(opt$SIMBA_GRN_targets, sep="\t", header=T)
  
  
  cat("SIMBA_GRN_targets_0\n")
  cat(str(SIMBA_GRN_targets))
  cat("\n")
  cat(str(unique(SIMBA_GRN_targets$Motif_Group)))
  cat("\n")
 
  
  SIMBA_GRN_targets$ENTREZID<-mapIds(org.Hs.eg.db, keys=SIMBA_GRN_targets$Gene_Name, keytype="SYMBOL",column="ENTREZID")
  SIMBA_GRN_targets$ensembl_gene_id<-mapIds(org.Hs.eg.db, keys=SIMBA_GRN_targets$Gene_Name, keytype="SYMBOL",column="ENSEMBL")
  
  cat("SIMBA_GRN_targets_1\n")
  cat(str(SIMBA_GRN_targets))
  cat("\n")
  
  SIMBA_GRN_targets_NO_NA<-SIMBA_GRN_targets[!is.na(SIMBA_GRN_targets$ENTREZID),]
  
  cat("SIMBA_GRN_targets_NO_NA_0\n")
  cat(str(SIMBA_GRN_targets_NO_NA))
  cat("\n")
  
  #### Prepare the file for gmt export ----
  
  for_gmt_SIMBA_GRN_targets_NO_NA<-unique(SIMBA_GRN_targets_NO_NA[,c(which(colnames(SIMBA_GRN_targets_NO_NA) == 'Motif_Group'),which(colnames(SIMBA_GRN_targets_NO_NA) == 'ENTREZID'))])
  
  colnames(for_gmt_SIMBA_GRN_targets_NO_NA)[which(colnames(for_gmt_SIMBA_GRN_targets_NO_NA) == 'Motif_Group')]<-'id'
  colnames(for_gmt_SIMBA_GRN_targets_NO_NA)[which(colnames(for_gmt_SIMBA_GRN_targets_NO_NA) == 'ENTREZID')]<-'gene'
  
  cat("for_gmt_SIMBA_GRN_targets_NO_NA_0\n")
  cat(str(for_gmt_SIMBA_GRN_targets_NO_NA))
  cat("\n")
  
  
  for_gmt_SIMBA_GRN_targets_NO_NA$name<-for_gmt_SIMBA_GRN_targets_NO_NA$id
  
  
  
  cat("for_gmt_SIMBA_GRN_targets_NO_NA_1\n")
  cat(str(for_gmt_SIMBA_GRN_targets_NO_NA))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(for_gmt_SIMBA_GRN_targets_NO_NA$name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(for_gmt_SIMBA_GRN_targets_NO_NA$name)))))
  cat("\n")
  
  #### Reorder
  
  indx.reorder<-c(which(colnames(for_gmt_SIMBA_GRN_targets_NO_NA) == 'id'),which(colnames(for_gmt_SIMBA_GRN_targets_NO_NA) == 'name'),which(colnames(for_gmt_SIMBA_GRN_targets_NO_NA) == 'gene'))
  
  
  for_gmt_SIMBA_GRN_targets_NO_NA_reordered<-for_gmt_SIMBA_GRN_targets_NO_NA[,indx.reorder]
  
  cat("for_gmt_SIMBA_GRN_targets_NO_NA_reordered_0\n")
  cat(str(for_gmt_SIMBA_GRN_targets_NO_NA_reordered))
  cat("\n")
  
 #### put together everything ----
  
  DEF<-rbind(for_gmt_new_sets_long_NO_NA_reordered,
             for_gmt_Table_of_gene_sets_long_NO_NA_reordered,
             for_gmt_mir5739_targets_reordered,
             for_gmt_CUX1_TARGET_GENES_df_reordered,
             for_gmt_doro_collapsed_reordered,
             for_gmt_SIMBA_GRN_targets_NO_NA_reordered)
  
  cat("DEF_0\n")
  cat(str(DEF))
  cat("\n")
 
 ### save as gmt ----
  
  setwd(out)
  
  writeGMT(DEF, paste("Custom_Soranzo",'_Hs.entrez.gmt',sep=''))
  
  
  
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
    make_option(c("--Table_of_gene_sets"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--new_sets"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--miRNA_table"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--c3_collection"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Dorothea_ABC"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Dorothea_ABCD"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--SIMBA_GRN_targets"), type="character", default=NULL, 
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
  
  create_gmt(opt)

  
}


###########################################################################

system.time( main() )