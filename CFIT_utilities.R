## This script is to make functions for processing CFIT files

find_nf_file <- function(cfit_prefix, sample_name, mid_directory="Neoantigen_qualities__OS"){
  nf_file_path <- paste(cfit_prefix, "/", mid_directory, "/nf_*.txt", sep="")
  nf_file_list <- Sys.glob(nf_file_path)
  return(nf_file_list[grep(sample_name, nf_file_list)])
}

read_nf_files <- function(cfit_prefix, sample_name_list, mid_directory="Neoantigen_qualities__OS"){
  neo_df <- NULL
  for (s in sample_name_list){
    nf_file_name <- find_nf_file(cfit_prefix, s, mid_directory)
    if ( is.null(nf_file_name) | length(nf_file_name)==0 ){
      cat("no matched nf file:", s, "\n")
    } else {
      for (f in nf_file_name){
        cat("reading", f, "...\n")
        neo_s <- read.table(f,
                            header=T, as.is=T, sep="\t")
        neo_s$sample_name <- s
        neo_df <- rbind(neo_df, neo_s)
      }
    }
  }
  neo_df$mutation <- sapply(neo_df$neoantigen, function(x) paste(strsplit(x, "_")[[1]][1:4], collapse="_") )
  return(neo_df)
}

find_fit_statistics_file <- function(cfit_prefix, sample_name){
  file_path <- paste(cfit_prefix, "/sample_statistics/fitness_sample_statistics_*.txt", sep="")
  file_list <- Sys.glob(file_path)
  return(file_list[grep(sample_name, file_list)])
}

read_fitness_statistics <- function(cfit_prefix, sample_name_list){
  fit_stat_df <- NULL
  for (s in sample_name_list){
    file_name <- find_fit_statistics_file(cfit_prefix, s)
    if ( is.null(file_name) | length(file_name)==0 ){
      cat("no matched fitness statistics file:", s, "\n")
    } else {
      for (f in file_name){
        cat("reading", f, "...\n")
        sample_df <- read.table(f,
                            header=T, as.is=T, sep="\t")
        sample_df$sample_name <- s
        fit_stat_df <- rbind(fit_stat_df, sample_df)
      }
    }
  }
  return(fit_stat_df)
}

find_CCF_file <- function(cfit_prefix, sample_name){
  file_path <- paste(cfit_prefix, "/*/*/Mutations/CCF_*.txt", sep="")
  file_list <- Sys.glob(file_path)
  return(file_list[grep(sample_name, file_list)])
}

read_CCF_files <- function(cfit_prefix, sample_name_list){
  ccf_df <- NULL
  for (s in sample_name_list){
    file_name <- find_CCF_file(cfit_prefix, s)
    if ( is.null(file_name) | length(file_name)==0 ){
      cat("no matched CCF file:", s, "\n")
    } else {
      for (f in file_name){
        cat("reading", f, "...\n")
        sample_df <- read.table(f,
                                header=T, as.is=T, sep="\t")
        sample_df$sample_name <- s
        ccf_df <- rbind(ccf_df, sample_df)
      }
    }
  }
  return(ccf_df)
}


find_mut_statistics_file <- function(cfit_prefix, sample_name){
  file_path <- paste(cfit_prefix, "/sample_statistics/sample_statistics_*.txt", sep="")
  file_list <- Sys.glob(file_path)
  return(file_list[grep(sample_name, file_list)])
}

read_mutation_statistics <- function(cfit_prefix, sample_name_list){
  mut_stat_df <- NULL
  for (s in sample_name_list){
    file_name <- find_mut_statistics_file(cfit_prefix, s)
    if ( is.null(file_name) | length(file_name)==0 ){
      cat("no matched fitness statistics file:", s, "\n")
    } else {
      for (f in file_name){
        cat("reading", f, "...\n")
        sample_df <- read.table(f,
                                header=T, as.is=T, sep="\t")
        sample_df$sample_name <- s
        mut_stat_df <- rbind(mut_stat_df, sample_df)
      }
    }
  }
  return(mut_stat_df)
}