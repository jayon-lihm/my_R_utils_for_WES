### tree utilities ###
source("utilities.R")

library("rjson")
library("data.tree")
library("DiagrammeR")
library(ggtree)
library("htmltools")
library("ggrepel")
library(paletteer)
library(rsvg)
library(scales)

## Read tree structure from json file
read_tree_from_json <- function(tree_file){
  tree <- fromJSON(file=tree_file)
  return(tree)
}

## convert it to Nodes and draw tree to HTML
phylo_tree_nodes <- function(tree, tree_html_file=""){
  phylo_tree <- as.Node(
    tree$trees[[1]]$topology,
    mode = "explicit",
    nameName = "clone_id",
    childrenName = "children",
    nodeName = NULL
  )
  if(tree_html_file != ""){
    cat("Drawing tree to HTML file:", tree_html_file, "\n")
    htmltools::save_html(html = plot(phylo_tree), file = tree_html_file)
  }
  return(phylo_tree)
}

## convert top tree into dataframe                                                                                                                                 
tree_to_data_frame <- function(phylo_tree){
  tree_df <- ToDataFrameNetwork(phylo_tree)
  tree_level <- phylo_tree$Get("level") - 1
  tree_nmuts <-  phylo_tree$Get(function(node) length(node$clone_mutations))
  tree_info <- merge(data.frame(clone=names(tree_level), tree_level), 
                     data.frame(clone=names(tree_nmuts), n_muts=tree_nmuts),
                     by="clone")
  tree_df <- merge(tree_df, tree_info, by.x="to", by.y="clone", all.x=TRUE)
  tree_df <- tree_df[order(tree_df$tree_level, tree_df$from, tree_df$to),]
  return(tree_df)
}

find_upper_nodes <- function(node, tree_df){
  ## "tree" should be for one patient
  ### account tree structure
  ### 0 -> 1 -> 2 -> 3
  ## child=3, upper_clones=0, 1, 2  
  
  child=node
  parent_list <- NULL
  flag=TRUE
  while (flag){
    parent=tree_df$from[tree_df$to==child]
    if (parent==0){
      parent_list <- c(parent, parent_list)
      flag=FALSE
      break
    } else {
      parent_list <- c(parent, parent_list)
      child <- parent
    }
  }
  return(parent_list)
}

tree_summary <- function(tree_file, sample_name, tree_html_file=""){
  ## Read tree structure from json file                                                                                                                            
  tree <- read_tree_from_json(tree_file)
  
  ## convert it to Nodes and draw tree to HTML                                                                                                                     
  tree_nodes <- phylo_tree_nodes(tree, tree_html_file)
  
  ## convert it to dataframe                                                                                                                                       
  tree_df <- tree_to_data_frame(tree_nodes)
  
  ## annotate upper nodes                                                                                                                                          
  for (i in 1:nrow(tree_df)){
    upper_nodes <- find_upper_nodes(tree_df$to[i], tree_df)
    tree_df$upper_nodes[i] <- paste(upper_nodes, collapse=",")
  }
  
  return(tree_df)
}

clone_assignmnet <- function(mut_df, tree_nodes){
  ## extract mutation IDs per clone                                                                                                                               
  clone_mut_ids <- tree_nodes$Get(function(node){ paste(node$clone_mutations, collapse=",") })
  
  ## sort mutation names                                                                                                                                          
  clone_mut_ids <- sapply(clone_mut_ids, function(x) paste(sort(strsplit(x, ",")[[1]]), collapse=","))
  
  clone_mut_ids_list <- lapply(clone_mut_ids, function(x) strsplit(x, ",")[[1]] )
  
  ## find clone name
  mut_df$clone_name <- sapply(mut_df$ID, function(x, y_list=clone_mut_ids_list){
    c_index <- which(unlist(lapply(y_list, function(y, mut=x){ mut %in% y })))
    if ( length(c_index)==0){
      clone_name <- "."
    } else {
      clone_name <- names(y_list)[c_index]
    }
    return(clone_name)
  })
  
  return(mut_df)
}

mutation_clone_summary <- function(sample_name, tree_df, tree_file, mutation_file, nonsyn_mut_types,
                                   file_type="vcf"){
  ## get phylogeny tree info                                                                                                                                       
  tree <- read_tree_from_json(tree_file)
  tree_nodes <- phylo_tree_nodes(tree)
  
  ## read mutation files
  if (file_type=="vcf"){
    mut_df <- read_vcf_file(mutation_file)
    mut_df <- parse_vcf_df(mut_df, sample_name=sample_name)
  }
  
  if (file_type=="final_mut_table"){
    mut_df <- read_final_mut_table(mutation_file)
  }
  
  mut_annot <- mutation_annotations(mut_df)
  mut_df$sample_name <- sample_name
  
  ## map mutID and clone                                                                                                                                           
  mut_annot <- clone_assignmnet(mut_annot, tree_nodes)
  
  ## compute mutation counts
  colnames(mut_annot)[colnames(mut_annot)=="mut_type"] <- "Mutation_Type"
  mut_clone_summary <- mut_annot %>%
    group_by(clone_name, sample_name) %>%
    summarize(                                                                                                                                         
      n_missense=sum( grepl("missense", Mutation_Type) ),
      n_frameshift=sum( grepl("frameshift", Mutation_Type) ),
      n_nonsyn=sum( grepl( paste(nonsyn_mut_types, collapse="|"), Mutation_Type)))
  
  return(tree_df %>% left_join(mut_clone_summary, by=c("to"="clone_name")))
}

###### neoantigen #######                                                                                                                                         
# read_neoag_files <- function(cfit_prefix, sample_name_list){
#   neo_df <- NULL
#   for (s in sample_name_list){
#     file_name <- Sys.glob(paste(cfit_prefix, "/Neoantigen_qualities_*/nf_", s, "_full.txt", sep=""))
#     neo_s <- read.table(file_name,
#                         header=T, as.is=T, sep="\t")
#     neo_s$sample_name <- s
#     neo_df <- rbind(neo_df, neo_s)
#   }
#   neo_df$mutation <- sapply(neo_df$neoantigen, function(x) paste(strsplit(x, "_")[[1]][1:4], collapse="_"), 
#                             USE.NAMES = FALSE)
#   return(neo_df)
# }

neoantigen_clone_summary <- function(tree_df, cfit_prefix, sample_name_list){
  ## neoantigen burden (<500nM), clone fitness, best quality, kD, A, R, C                                                                                         
  neo_df <- read_neoag_files(cfit_prefix, sample_name_list)
  
  x <- neo_df %>%
    group_by(sample_name, clone_number) %>%
    summarize(clone_fitness=unique(clone_fitness),
              CCF=unique(X), exc_CCF=unique(Y),
              strong_binders=sum(kDmt<500),
              all_binders=n(),
              min_kD=min(kDmt),
              max_A=max(A),
              max_R=max(R),
              max_D=max(D),
              max_quality=max(quality),
              max_AR=max(A*R)) %>%
    mutate(clone_number.c=as.character(clone_number))
  
  tree_df <- tree_df %>%
    left_join(x, by=c("to"="clone_number.c", "sample_name"="sample_name")) %>%
    arrange(as.numeric(to))
  
  return(tree_df)
}

get_CCF_mutations_from_clone <- function(tree_node){
  clone_info <- tree_node$Get(function(node) paste(node$X, ## CCF
                                                   node$x, ## exclusive CCF
                                                   paste(node$clone_mutations, collapse=","), ## mutations
                                                   sep=";") )
  #print(clone_info)
  return(clone_info)
}

tabularize_clone_info <- function(clone_info){
  temp <- strsplit(clone_info, ";")
  #print(temp)
  clone_names <- names(temp)
  clone_df <- NULL
  for (i in 1:length(temp)){
    clone_i <- clone_names[i]
    if (clone_i != "0"){
      clone_CCF <- as.numeric(temp[[i]][1])
      clone_excCCF <- as.numeric(temp[[i]][2])
      clone_muts <- strsplit(temp[[i]][3], ",")[[1]]
      
      clone_df <- rbind(clone_df, data.frame(mutations=clone_muts, clone_name=clone_i, 
                                             CCF=clone_CCF, excCCF=clone_excCCF))
    }
  }
  return(clone_df)
}

retrieve_clone_info <- function(tree, rank=1){
  ## tree -> read_tree_from_json
  ## get sample list
  num_samples <- length(tree$time_points[[1]]$samples)
  patient_clone_df <- NULL
  for (i in 1:num_samples){
    sname <- tree$time_points[[1]]$samples[[i]]$id
    ## best tree
    best_tree <- tree$time_points[[1]]$samples[[i]]$sample_trees[[rank]]
    
    best_phylo_tree <- as.Node(
      best_tree$topology,
      mode = "explicit",
      nameName = "clone_id",
      childrenName = "children",
      nodeName = NULL)
    
    clone_info <- get_CCF_mutations_from_clone(best_phylo_tree)
    clone_df <- tabularize_clone_info(clone_info)
    clone_df$sample_name <- sname
    patient_clone_df <- rbind(patient_clone_df, clone_df)
  }
  return(patient_clone_df)
}

retrieve_clone_info_for_short_json <- function(tree, tree_rank=1){
  ## tree -> read_tree_from_json
  ## get sample list
  num_samples <- length(tree)
  patient_clone_df <- NULL
  for (i in 1:num_samples){
    sname <- tree[[i]]$id
    ## best tree
    best_tree <- tree[[i]]$sample_trees[[tree_rank]]
    
    best_phylo_tree <- as.Node(
      best_tree$topology,
      mode = "explicit",
      nameName = "clone_id",
      childrenName = "children",
      nodeName = NULL)
    
    clone_info <- get_CCF_mutations_from_clone(best_phylo_tree)
    clone_df <- tabularize_clone_info(clone_info)
    clone_df$sample_name <- sname
    patient_clone_df <- rbind(patient_clone_df, clone_df)
  }
  return(patient_clone_df)
}

ColorScaleFunction <- function(Range, high = "red", low = "blue", ValueToLookup) {
  negative_colors <- seq_gradient_pal(low, "white")
  positive_colors <- seq_gradient_pal("white", high)
  if ( ValueToLookup == 0){
    color <- "white"
  } else if (ValueToLookup > 0){
    color <- positive_colors( ValueToLookup / Range[2] )
  } else{
    color <- negative_colors( 1 - abs(ValueToLookup) / abs(Range[1]) )
  }
  return(color)
}
