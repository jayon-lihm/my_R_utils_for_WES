## 2024/2/21 

library(readxl)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyr)
library("ggpubr")
library(ComplexHeatmap) ## for oncoprint
#library("ensembldb") ## for oncoprint
#library(EnsDb.Hsapiens.v86) ## for oncoprint hg38
#library(EnsDb.Hsapiens.v75) ## for oncoprint hg19
library(ggbeeswarm)

## final mutation information table header
mut_header <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "Patient_Name", "Normal_Name", "Tumor_Name",
                "Gene_Name", "Gene_ID", "Feature_Type", "Feature_ID", "Transcript_BioType", "Mutation_Type", "Putative_Impact",
                "DNA_Mutation", "Amino_Acid_Mutation", "CDS_position.CDS_length", "Protein_position.Protein_length", "Callers",
                "NORMAL_Tot_Cov", "NORMAL_Ref_Cov", "NORMAL_Alt_Cov", "NORMAL_VAF.pct", "TUMOR_Tot_Cov", "TUMOR_Ref_Cov", "TUMOR_Alt_Cov", "TUMOR_VAF.pct")

text_vcf_header <- "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR"
vcf_header_columns <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")

b37_chrom_names <- c(1:22, "X", "Y")
human_canonical_chrom_names <- c(1:22, "X", "Y")

axis_rotation <- theme(axis.text.x = element_text(angle=34, hjust=1))

nonsyn_types <- c("missense", "frameshift", "disruptive", "stop_gained", "stop_lost", "start_lost", "splice")

get_first_line <- function(file_name){
  if ( endsWith(file_name, ".gz") ){
    cmd <- paste("gzcat", file_name, " | grep CHROM")
  } else {
    cmd <- paste("grep CHROM", file_name)
  }
  first_line <- system(cmd, intern = TRUE)
  return(first_line)
}

get_vcf_header <- function(vcf_file){
  first_line <- get_first_line(vcf_file)
  header <- strsplit(first_line, "\t")[[1]]
  header[1] <- gsub("#", "", header[1])
  return(header)
}

read_vcf_file <- function(vcf_file){
  if ( endsWith(vcf_file, ".gz") ){
    vcf_df <- read.table( gzfile(vcf_file), header=F, as.is=T, sep="\t")
  } else {
    vcf_df <- read.table(vcf_file, header=F, as.is=T, sep="\t")
  }
  colnames(vcf_df) <- get_vcf_header(vcf_file)
  return(vcf_df)
}

read_final_mut_table <- function(final_mut_file, sample_name=NULL){
  mut_df <- read.table(final_mut_file, header=F, as.is=T, sep="\t", col.names = mut_header)

  ## remove MT chroms                                                                                                                                              
  if (any(mut_df$CHROM %in% c("M", "MT", "chrM", "chrMT"))){
    mut_df <- mut_df[-which(mut_df$CHROM %in% c("M", "MT", "chrM", "chrMT")),]
  }
  
  if (!is.null(sample_name)){
    mut_df$sample_name <- sample_name
  }
  return(mut_df)
}

extract_ann_field <- function(info_column){
  info_list <- strsplit(info_column, ";")[[1]]
  ann_field <- info_list[which(startsWith(info_list, "ANN=") |
                               substring(info_list, 1, 2) %in% c("A|", "C|", "G|", "T|") )]
  return(select_ann_field(ann_field) )
}

select_ann_field <- function(ann_field){
  ann_list <- strsplit(gsub("ANN=", "", ann_field), ",")[[1]]
  mutType_list <- sapply(ann_list, "get_mutation_type", USE.NAMES = FALSE)
  if (any(grepl("missense", mutType_list))){
    selected_ann_i <- grep("missense", mutType_list)
  } else {
    selected_ann_i <- 1
  }
  return(ann_list[selected_ann_i])
}


get_mutation_type <- function(ann){
  return(strsplit(ann, "\\|")[[1]][2])
}

get_gene_name <- function(ann){
  return(strsplit(ann, "\\|")[[1]][4])
}

get_gene_ensembl_ID <- function(ann){
  return(strsplit(ann, "\\|")[[1]][5])
}

get_mutation_region <- function(ann){
  return(strsplit(ann, "\\|")[[1]][6])
}

get_transcript_ensembl_ID <- function(ann){
  return(strsplit(ann, "\\|")[[1]][7])
}

get_amino_acid_change <- function(ann){
  return(strsplit(ann, "\\|")[[1]][11])
}

get_total_cov <- function(x){
  total_cov <- strsplit(x, ":")[[1]][1]
  if ( total_cov == "."){
    total_cov <- 0
  } else {
    total_cov <- as.numeric(total_cov)
  }
  return( total_cov )
}

get_ref_cov <- function(x){
  ref_cov <- strsplit(x, ":")[[1]][2]
  if ( ref_cov == "."){
    ref_cov <- 0
  } else {
    ref_cov <- as.numeric(ref_cov)
  }
  return( ref_cov )
}

get_AD_field.Dragen <- function(x){
  return(unlist(strsplit(x, ":"))[3])
}

get_ref_cov.Dragen <- function(x){
  AD <- get_AD_field.Dragen(x)
  ref <- strsplit(AD, ",")[[1]][1]
  return(as.numeric(ref))
}

get_alt_cov.Dragen <- function(x){
  AD <- get_AD_field.Dragen(x)
  alt <- strsplit(AD, ",")[[1]][2]
  return(as.numeric(alt))
}

parse_vcf_df <- function(vcf_df, sample_name=NULL, vcf_format="Greenbaum"){
  ## parse info field
  x <- sapply(vcf_df$INFO, "extract_ann_field", USE.NAMES = FALSE)
  
  mut_type <- sapply(x, "get_mutation_type", USE.NAMES = FALSE)
  gene_name <- sapply(x, "get_gene_name", USE.NAMES = FALSE)
  gene_id <- sapply(x, "get_gene_ensembl_ID", USE.NAMES = FALSE)
  mut_region <- sapply(x, "get_mutation_region", USE.NAMES = FALSE)
  transcript_id <- sapply(x, "get_transcript_ensembl_ID", USE.NAMES = FALSE)
  AA_change <- sapply(x, "get_amino_acid_change", USE.NAMES = FALSE)
  
  ## parse normal and tumor coverage
  if ( toupper(vcf_format) == "GREENBAUM"){
    normal_total <- sapply(vcf_df$NORMAL, "get_total_cov", USE.NAMES = FALSE)
    normal_ref <- sapply(vcf_df$NORMAL, "get_ref_cov", USE.NAMES = FALSE)
    normal_alt <- normal_total - normal_ref
    tumor_total <- sapply(vcf_df$TUMOR, "get_total_cov", USE.NAMES = FALSE)
    tumor_ref <- sapply(vcf_df$TUMOR, "get_ref_cov", USE.NAMES = FALSE)
    tumor_alt <- tumor_total - tumor_ref
  }
  
  if ( toupper(vcf_format) == "DRAGEN"){
    normal_ref <- sapply(vcf_df$NORMAL, "get_ref_cov.Dragen", USE.NAMES = FALSE)
    normal_alt <- sapply(vcf_df$NORMAL, "get_alt_cov.Dragen", USE.NAMES = FALSE)
    normal_total <- normal_ref + normal_alt
  
    tumor_ref <- sapply(vcf_df$TUMOR, "get_ref_cov.Dragen", USE.NAMES = FALSE)
    tumor_alt <- sapply(vcf_df$TUMOR, "get_alt_cov.Dragen", USE.NAMES = FALSE)
    tumor_total <- tumor_ref + tumor_alt
  }
  
  new_df <- cbind(vcf_df[, c("CHROM", "POS", "ID", "REF", "ALT")], 
                  data.frame(mut_type, gene_name, gene_id, mut_region, transcript_id, AA_change, 
                           normal_total, normal_ref, normal_alt, tumor_total, tumor_ref, tumor_alt))
  
  if (!is.null(sample_name)){
    new_df$sample_name <- sample_name
  }
  return(new_df)
}

mutation_annotations <- function(mut_df){
  ## trim sample specific columns (such as VAF, coverage...)                                                                                                      
  x <- mut_df %>% dplyr::select(-contains("NORMAL"), -contains("TUMOR"), 
                                -contains("normal"), -contains("tumor")) %>%
    distinct()
  return(x)
}

tag_filtering_status <- function(vcf_df, n_normal=1, ignore_normal=FALSE, roundTo=NULL){
  ## vcf_df: after "parse_vcf_df" function
  low_conditions <- c("t_tot_thres"=10, "t_alt_thres"=5, "t_vaf_thres"=0.02,
                      "n_tot_thres"=7*n_normal, "n_vaf_thres"=0.01)
  strict_conditions <- c("t_tot_thres"=10, "t_alt_thres"=9, "t_vaf_thres"=0.04,
                         "n_tot_thres"=7*n_normal, "n_vaf_thres"=0.01)
  
  if (ignore_normal){
    vcf_df2 <- vcf_df %>% 
      mutate(tumor_alt = tumor_total - tumor_ref,
             tumor_vaf = tumor_alt / tumor_total)
    if ( !is.null(roundTo)){
      vcf_df2$tumor_vaf <- round(vcf_df2$tumor_vaf, roundTo)
    }
  } else {
    vcf_df2 <- vcf_df %>% 
      mutate(normal_alt = normal_total - normal_ref,
             normal_vaf = normal_alt / normal_total,
             tumor_alt = tumor_total - tumor_ref,
             tumor_vaf = tumor_alt / tumor_total) 
    if ( !is.null(roundTo)){
      vcf_df2$tumor_vaf <- round(vcf_df2$tumor_vaf, roundTo)
      vcf_df2$normal_vaf <- round(vcf_df2$normal_vaf, roundTo)
    }
  }
  
  if (ignore_normal){
    x <- vcf_df2 %>%
      mutate(filtering = case_when(
          tumor_total >= strict_conditions["t_tot_thres"] & 
          tumor_alt >= strict_conditions["t_alt_thres"] & 
          tumor_vaf >= strict_conditions["t_vaf_thres"] ~ "passed(str)",
        
          tumor_total >= low_conditions["t_tot_thres"] & 
          tumor_alt >= low_conditions["t_alt_thres"] & 
          tumor_vaf >= low_conditions["t_vaf_thres"] ~ "passed(low)",
        
        tumor_vaf == 0 ~ "no_alt_allele",
        TRUE ~ "low_qual")) %>%
      pull(filtering)
  } else {
    x <- vcf_df2 %>% 
      mutate(filtering = case_when(
          normal_total >= strict_conditions["n_tot_thres"] & 
          normal_vaf <= strict_conditions["n_vaf_thres"] &
          normal_vaf >= 0 &
          tumor_total >= strict_conditions["t_tot_thres"] & 
          tumor_alt >= strict_conditions["t_alt_thres"] & 
          tumor_vaf >= strict_conditions["t_vaf_thres"] ~ "passed(str)",
        
        normal_total >= low_conditions["n_tot_thres"] & 
          normal_vaf <= low_conditions["n_vaf_thres"] &
          normal_vaf >= 0 &
          tumor_total >= low_conditions["t_tot_thres"] & 
          tumor_alt >= low_conditions["t_alt_thres"] & 
          tumor_vaf >= low_conditions["t_vaf_thres"] ~ "passed(low)",
        
        tumor_vaf == 0 ~ "no_alt_allele",
        TRUE ~ "low_qual")) %>%
      pull(filtering)
  }
  return(x)
}

filter_nonsyn_mut_types <- function(vcf_df, nonsyn_types){
  a <- vcf_df %>%
    dplyr::filter( grepl(paste(nonsyn_types, collapse="|"), mut_type ))
  return(a)
}

filter_genes <- function(vcf_df, gene_list, exclude_intergenic=TRUE){
  ## by default: exclude intergenic regions
  
  gene_list2 <- paste("\\<", gene_list, "\\>", sep="")
  if (exclude_intergenic){
    ftd <- vcf_df %>%
      dplyr::filter( !grepl("intergenic", mut_type)) %>%
      dplyr::filter( grepl(paste(gene_list2, collapse="|"), gene_name ))
  } else {
    ftd <- vcf_df %>%
      filter( grepl(paste(gene_list2, collapse="|"), gene_name))
  }
  return(ftd)
}

summarize_mutation_counts <- function(vcf_df, nonsyn_types=NULL){
  ## "vcf_df" from parse_vcf_df function
  
  mut_count_summary <- vcf_df %>%
    group_by(sample_name) %>%
    summarize(total_muts=n())
  
  if ( !is.null(nonsyn_types) ) {
    nonsyn_counts <- vcf_df %>%
      filter( grepl(paste(nonsyn_types, collapse="|"), mut_type )) %>%
      group_by(sample_name) %>%
      summarize(nonsyn_muts=n())
    mut_count_summary <- merge(mut_count_summary, nonsyn_counts, by="sample_name")
  }
  
  return(mut_count_summary)
}

summarize_gene_set_alterations <- function(vcf_df, gene_list, gene_set_name=NULL, nonsyn_types=NULL){
  ## This function is to check how many genes and how many mutations are altered
  ## in a given gene set
  
  sample_list <- unique(vcf_df$sample_name)
  
  ## get nonsyn muts
  if (!is.null(nonsyn_types)) {
    selected_vcf_df <- filter_nonsyn_mut_types(vcf_df, nonsyn_types)
  } else {
    selected_vcf_df <- vcf_df
  }
  
  ## subset genes
  gene_vcf_df <- filter_genes(selected_vcf_df, gene_list)
  
  ## count mutations per gene (row: sample, column: gene)
  mut_counts <- gene_vcf_df %>% 
    group_by(sample_name, gene_name) %>%
    summarize(n_nonsyn_muts=n()) %>%
    spread(gene_name, n_nonsyn_muts) %>%
    replace(is.na(.), 0)

  ## fill in empty genes
  genes_zero_counts <- setdiff(gene_list, colnames(mut_counts))
  mut_counts[, genes_zero_counts] <- 0
  
  ## fill in empty samples and 
  empty_samples <- 
    data.frame(sample_name=setdiff(sample_list, mut_counts$sample_name))
  empty_samples[, colnames(mut_counts)[-1]] <- 0
  
  mut_counts <- rbind(mut_counts, empty_samples)
  
  ## count how many mutations across the gene set
  num_muts_altered <- apply( mut_counts[, gene_list], 1, function(x){
    sum(x)
  })
  ## count number of genes that are affected
  num_genes_altered <- apply( mut_counts[, gene_list], 1, function(x){
    sum(x>0)
  })
  
  if (!is.null(gene_set_name)){
    mut_counts[, paste("num_muts_altered", gene_set_name, sep=".")] <- num_muts_altered
    mut_counts[, paste("num_genes_altered", gene_set_name, sep=".")] <- num_genes_altered
  } else {
    mut_counts$num_muts_altered <- num_muts_altered
    mut_counts$num_genes_altered <- num_genes_altered
  }
  return(mut_counts)
}

compute_TMB <- function(vcf_df, TMB_mut_type=NULL){
  ## Compute TMB from a parsed VCF dataframe
  ## If TMB_mut_type is given, consider those mutations for TMB computation
  y <- NULL
  if (!is.null(TMB_mut_type)) {
    y <- vcf_df %>%
      filter( grepl(paste(TMB_mut_type, collapse="|"), mut_type) ) %>%
      group_by(sample_name) %>%
      summarize(n_muts=n())
  } else {
    y <- vcf_df %>%
      group_by(sample_name) %>%
      summarize(n_muts=n())
  }
  return(y)
}

create_vcf_df <- function(vcf_list, file_pattern){
  all_vcf_df <- NULL
  for (v in vcf_list){
    raw_vcf_df <- read_vcf_file(v)
    vcf_df <- parse_vcf_df(raw_vcf_df, sample_name=get_sample_name_from_vcf_file(v, file_pattern=file_pattern))
    all_vcf_df <- rbind(all_vcf_df, vcf_df)
  }
  
  ## compute VAF
  all_vcf_df$normal_vaf <- (all_vcf_df$normal_total - all_vcf_df$normal_ref) / all_vcf_df$normal_total
  all_vcf_df$tumor_vaf <- (all_vcf_df$tumor_total - all_vcf_df$tumor_ref) / all_vcf_df$tumor_total
  
  return(all_vcf_df)
}


get_sample_name_from_vcf_file <- function(vcf_file_name, file_pattern){
  return(gsub(file_pattern, "", basename(vcf_file_name)))
}

write_header <- function(vcf_file_name, header){
  writeLines(paste(header, sep="\t"), vcf_file_name)
}

write_vcf_from_vcf_df <- function(vcf_df, file_name_suffix="filtered", output_dir="./",
                                  create_sample_dir=FALSE){
  sample_list <- unique(vcf_df$sample_name)
  select_columns <- c("CHROM", "POS", "ID", "REF", "ALT", "normal_total", "normal_ref", "tumor_total", "tumor_ref")
  
  for (s in sample_list){
    ## If create_sample_dir=TRUE, create sample_name directory and output 
    if (create_sample_dir){
      new_output_dir <- paste(output_dir, "/", s, sep="")
    } else {
      new_output_dir <- output_dir
    }
    
    if (! dir.exists(new_output_dir)){
      dir.create(new_output_dir, recursive=TRUE)
    }
    
    output_file <- paste(new_output_dir, "/", s, "_", file_name_suffix, ".vcf", sep="")
    
    this_sample_vcf_df <- vcf_df[which(vcf_df$sample_name == s), select_columns]
    ## make VCF folumns
    this_sample_vcf_df$QUAL <- "."
    this_sample_vcf_df$FILTER <- "PASS"
    this_sample_vcf_df$INFO <- "INFO"
    this_sample_vcf_df$FORMAT <- "DP:AP"
    this_sample_vcf_df$NORMAL <- paste(this_sample_vcf_df$normal_total, this_sample_vcf_df$normal_ref, sep=":")
    this_sample_vcf_df$TUMOR <- paste(this_sample_vcf_df$tumor_total, this_sample_vcf_df$tumor_ref, sep=":")
    
    write_header(output_file, text_vcf_header)
    write.table(this_sample_vcf_df[, vcf_header_columns], output_file,
                col.names=F, row.names=F, quote=F, sep="\t",
                append = TRUE)
  }
}

##########################################
########### ONCOPRINT ########

recoding_mutation_type <- function(mutation_type){
  recoded_mutation_type <- NULL
  ## Order of selection: 
  ## frameshift > missense > stop_lost
  ## > stop_gain, start_lost 
  ## > disruptive_inframe > splice
  
  if (grepl("frameshift", mutation_type)){
    recoded_mutation_type <- "frameshift"
  } else if ( grepl("missense", mutation_type)){
    recoded_mutation_type <- "missense"
  } else if ( grepl("stop_lost", mutation_type)){
    recoded_mutation_type <- "stop_lost"
  } else if ( grepl("stop_gain", mutation_type)){
    recoded_mutation_type <- "stop_gain"
  } else {
    recoded_mutation_type <- "other"
  }
  
}

collect_muts_by_gene <- function(vcf_df, recoding=TRUE){
  ## This function collects muts by gene
  ## so that each sample can have one row per gene
  ## Multiple mutations are collapsed by ";"
  ## Order of selection: 
  ## frameshift > missense > stop_lost
  ## > stop_gain, start_lost 
  ## > disruptive_inframe > splice
  
  if (recoding){
    vcf_df$mut_type <- sapply(vcf_df$mut_type, "recoding_mutation_type", USE.NAMES = FALSE)
  }
  
  out <- vcf_df %>%
    dplyr::select(gene_name, mut_type, sample_name) %>%
    group_by(sample_name, gene_name) %>%
    summarize(mut_types=paste( unique(mut_type), collapse=";"))
  
  return(out)
}


long_to_wide_by_gene <- function(vcf_df, nonsyn_types=NULL, recoding=TRUE){
  ## column: gene
  ## row: sample
  ## value: mut types (separated by ";")
  
  ## convert format
  gene_simple_mut_df <- collect_muts_by_gene(vcf_df, recoding=recoding)
  
  ## spread it to matrix
  out <- gene_simple_mut_df %>% 
    dplyr::select(gene_name, mut_types, sample_name) %>%
    spread(gene_name, mut_types) %>%
    replace(is.na(.), " ")
  return(out)
}


convert_to_oncoprint_mat <- function(vcf_df, gene_list, nonsyn_types=NULL, recoding=TRUE){
  ## [Expected output]
  ## row: gene name
  ## column: sample name
  ## 1) filter for nonsynonymous mutations
  ## 2) select genes
  ## 3) convert to matrix
  ## 4) fill in empty samples (samples without any mutation under given conditions)
  sample_list <- unique(vcf_df$sample_name)
  
  ## filter for mutation types
  if ( ! is.null(nonsyn_types) ){
    selected_vcf_df <- filter_nonsyn_mut_types(vcf_df, nonsyn_types)
  } else {
    selected_vcf_df <- vcf_df
  }
  
  ## filter for genes
  gene_vcf_df <- filter_genes(selected_vcf_df, gene_list)
  
  ## spread to matrix
  gene_mut_mat <- long_to_wide_by_gene(gene_vcf_df, nonsyn_types = nonsyn_types,
                                       recoding=recoding)
  sample_row_names <- gene_mut_mat$sample_name
  gene_mut_mat <- as.matrix(gene_mut_mat[, -1])
  rownames(gene_mut_mat) <- sample_row_names
  
  ## add empty samples
  samples_with_no_muts <- setdiff(sample_list, sample_row_names)
  empty_mat <- matrix(" ", nrow=length(samples_with_no_muts), ncol=ncol(gene_mut_mat))
  rownames(empty_mat) <- samples_with_no_muts
  colnames(empty_mat) <- colnames(gene_mut_mat) ## gene
  gene_mut_mat <- rbind(gene_mut_mat, empty_mat)
  
  return(t(gene_mut_mat))
}

## https://www.cbioportal.org/oncoprinter

assign_oncoprinter_alteration_type <- function(mut_type){
  ##alteration_list <- c("MISSENSE", "INFRAME", "TRUNC", "PROMOTER", "SPLICE", "OTHER")
  
  if ( grepl( paste("frameshift", "stop_gained", "start_lost", "stop_lost", sep="|"), mut_type) ){
    alt_type <- "TRUNC"
  } else if (grepl ("missense", mut_type)){
    alt_type <- "MISSENSE"
  } else if (grepl("splice", mut_type)){
    alt_type <- "SPLICE"
  } else if (grepl("inframe", mut_type)) {
    alt_type <- "INFRAME"
  } else {
    alt_type <- "OTHER"
  }
   
  return(alt_type)
  
}


#if (ref_genome %in% c("hg19", "hg37", "b37")){
#  EnsDb <- EnsDb.Hsapiens.v75
#}

#if (ref_genome %in% c("hg38")){
#  EnsDb <- EnsDb.Hsapiens.v86
#}


# get_splice_protein_location <- function(chrom, position, transcript_id){
#   #https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/coordinate-mapping.html#Mapping_between_genome,_transcript_and_protein_coordinates
#   
#   # vcf_df[grep("splice", vcf_df$mut_type)[2],]
#   chrom <- "12"
#   position <- 49433216
#   edbx <- filter(EnsDb, filter = ~ seq_name == chrom)
#   gnm <- GRanges(chrom, IRanges(start = c(position-2),
#                               width = c(4)) )
#   
#   ## BRCA2, 13,	32893213, G>C, 	c.68-1G>C, p.X23_splice, ENST00000380152
#   chrom <- "13"
#   position <- 32893213
#   edbx <- filter(EnsDb, filter = ~ seq_name == chrom)
#   edbx <- filter(EnsDb, filter = ~ seq_name == chrom)
#   
#   
#   gnm <- GRanges(chrom, IRanges(start = c(position-10),
#                                 width = c(20)) )
#   gnm_prt <- genomeToProtein(gnm, EnsDb)
#   
#   transcript_id <- "ENST00000380152"
#   prts <- proteins(edbx, filter = TxIdFilter(transcript_id),
#                   return.type = "AAStringSet")
# 
# }



get_oncoprinter_genomic_alterations <- function(vcf_df){
  ## Expected output
  ## Sample Gene Alteration Type
  
  alteration_list <- c("MISSENSE", "INFRAME", "TRUNC", "PROMOTER", "SPLICE", "OTHER")
  vcf_df$Alteration <- vcf_df$AA_change
  vcf_df$Type <- NA

  for (i in 1:nrow(vcf_df)){
    vcf_df$Type[i] <- assign_oncoprinter_alteration_type(vcf_df$mut_type[i])
    if (vcf_df$Type[i] == "SPLICE"){
      vcf_df$Alteration[i] <- "splice"
    }
  }
  return(vcf_df[, c("sample_name", "gene_name", "Alteration", "Type")])
}

create_oncoprinter_genomic_input <- function(vcf_df, gene_list, nonsyn_types=NULL){
  ## This is for web oncoprinter data prep
  
  ## filter for mutation types
  if ( ! is.null(nonsyn_types) ){
    selected_vcf_df <- filter_nonsyn_mut_types(vcf_df, nonsyn_types)
  } else {
    selected_vcf_df <- vcf_df
  }
  
  ## filter for genes
  gene_vcf_df <- filter_genes(selected_vcf_df, gene_list)
  
  ## Oncoprinter input table
  onco_alt_table <- get_oncoprinter_genomic_alterations(gene_vcf_df)
  
  ## original sample list
  ori_sample_list <- unique(vcf_df$sample_name)
  
  ## if some samples are missing, fill in with empty entries
  if ( !all(ori_sample_list %in% onco_alt_table$sample_name)){
    fill_in_sample <- setdiff(ori_sample_list, onco_alt_table$sample_name)
    onco_alt_table <- rbind(onco_alt_table, data.frame(sample_name=fill_in_sample, gene_name="", Alteration="", Type=""))
  }
  
  return(onco_alt_table)
}

get_cmoID_patient <- function(x){
  y <- sapply(x, function(x){
    return(paste(strsplit(x, "_")[[1]][2:3], collapse="_"))
  }, USE.NAMES = FALSE)
  return(y)
}

get_cmoID_status <- function(x){
  y <- sapply(x, function(x){
    status <- strsplit(x, "_")[[1]][4]
    if (substr(status, 1, 1) == "N"){
      status2 <- "Normal"      
    } else if (substr(status, 1, 1)== "P"){
      status2 <- "Primary"
    } else if (substr(status, 1, 1) == "M"){
      status2 <- "Mets"
    }
    return(status2)
  }, USE.NAMES = FALSE)
  return(y)
}

insert_new_mut2vcfDf <- function(new_id, read_replacement_df){
  read_replacement_df[which(read_replacement_df$mutation_ID == new_id), ]
}

replace_reads <- function(vcf_df, read_replacement_df, replace_priority="higher_vaf"){
  vcf_df$in_IMPACT_table <- FALSE
  vcf_df$read_replaced <- FALSE
  for (i in 1:nrow(read_replacement_df)){
    ## first check if mutation in read_replacement_df is in vcf_df
    if ( read_replacement_df$mutation_ID[i] %in% vcf_df$ID |
         read_replacement_df$mutation_ID[i] %in% vcf_df$ID2 ){
      
      ## extract info
      replacement_mutid = read_replacement_df$mutation_ID[i]
      vcf_df_index <- unique(which(vcf_df$ID == read_replacement_df$mutation_ID[i] | vcf_df$ID2 == read_replacement_df$mutation_ID[i] ))
      
      vcf_df$in_IMPACT_table[vcf_df_index] <- TRUE
      
      ## update only if vaf is higher
      if (replace_priority == "higher_vaf"){
        if ( read_replacement_df$t_vaf[i] > vcf_df$tumor_vaf[vcf_df_index] ){
          vcf_df$tumor_alt[vcf_df_index] <- read_replacement_df$t_alt_count[i]
          vcf_df$tumor_ref[vcf_df_index] <- read_replacement_df$t_ref_count[i]
          vcf_df$tumor_total[vcf_df_index] <- vcf_df$tumor_alt[vcf_df_index] + vcf_df$tumor_ref[vcf_df_index]
          vcf_df$read_replaced[vcf_df_index] <- TRUE
        }
      }
    } else {
      ## if mutation in read_replacement_df is NOT in "vcf_df", then add a new entry to "vcf_df"
      ## new entries
      new_vcf_df = as.data.frame(matrix(NA, nrow=1, ncol=ncol(vcf_df)))
      colnames(new_vcf_df) <- colnames(vcf_df)
      new_vcf_df$ID <- replacement_mutid
      new_vcf_df$tumor_ref <- read_replacement_df$t_ref_count[i]
      new_vcf_df$tumor_total <- read_replacement_df$t_alt_count[i] + read_replacement_df$t_ref_count[i]
      new_vcf_df$read_replaced <- TRUE
      ## add
      vcf_df <- rbind(vcf_df, new_vcf_df)
    }
  }
  
  return(vcf_df)
}
  
  
vcfDf_to_ssm <- function(vcf_df, ssm_file, sex="Female", read_replacement_df=NULL, replace_priority="higher_vaf"){
  ## This function converts vcf_df (vcf file parsed) to SSM file for phylogeny tree reconstruction.
  ## for phylowgs
  
  ## ssm_file: output
  ## no return value
  header <- c("id", "gene", "a", "d", "mu_r", "mu_v")
  ssm_df <- data.frame(id=paste("s", 1:nrow(vcf_df) - 1, sep=""),
             gene=vcf_df$ID,
             a = vcf_df$tumor_ref,
             d = vcf_df$tumor_total,
             mu_r = 0.999,
             mu_v = 0.499
             )
  
  ## correct for female male
  if (sex == "Male" | "Y" %in% vcf_df$CHROM){
    ssm_df$mu_v[which(ssm_df$CHROM=="X")] = 0.001
    ssm_df$mu_v[which(ssm_df$CHROM=="Y")] = 0.001
  }
  
  ## replace with IMPACT (or any other)
  if ( !is.null(read_replacement_df) & nrow(read_replacement_df) > 0){
    cat("read counts are updated by given read_replacement_df...\n")
    cat(nrow(read_replacement_df), "mutations are in read_replacement_df...\n")
    vcf_df <- replace_reads(vcf_df, read_replacement_df, replace_priority=replace_priority)
  } else {
    cat("No mutations are in read_replacement_df...\n")
  }
  
  cat("Writing ", ssm_file, "\n")
  write.table(ssm_df, ssm_file, col.names=T, row.names=F, quote=F, sep="\t")
  return(vcf_df)
}


vcfDf_to_orchardSSM <- function(vcf_df, ssm_file, sex="Female", 
                                read_replacement_df=NULL, replace_priority="higher_vaf"){
  ## This function converts vcf_df (vcf file parsed) to SSM file
  ## for pairtree, or orchard
  
  ## ssm_file: output
  ## no return value
  header <- c("id", "name", "var_reads", "total_reads", "var_read_prob")
  ssm_df <- data.frame(id=paste("s", 1:nrow(vcf_df) - 1, sep=""),
                       name=vcf_df$ID,
                       var_reads = vcf_df$tumor_alt,
                       total_reads = vcf_df$tumor_total,
                       var_read_prob = 0.5
  )
  
  ## correct for  male one copy 
  if (sex == "Male" | "Y" %in% vcf_df$CHROM){
    ssm_df$var_read_prob[which( startsWith(vcf_df$CHROM, "X") | startsWith(vcf_df$CHROM, "chrX"))] = 1
    ssm_df$var_read_prob[which( startsWith(vcf_df$CHROM, "Y") | startsWith(vcf_df$CHROM, "chrY"))] = 1
  }
  
  ## replace with IMPACT (or any other)
  if ( !is.null(read_replacement_df) & nrow(read_replacement_df) > 0){
    cat("read counts are updated by given read_replacement_df...\n")
    cat(nrow(read_replacement_df), "mutations are in read_replacement_df...\n")
    vcf_df <- replace_reads(vcf_df, read_replacement_df, replace_priority=replace_priority)
  } else {
    cat("No mutations are in read_replacement_df...\n")
  }
  
  cat("Writing ", ssm_file, "\n")
  write.table(ssm_df, ssm_file, col.names=T, row.names=F, quote=F, sep="\t")
  return(vcf_df)
}

