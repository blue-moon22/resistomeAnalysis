#' Combines mapping data with metadata and CARD data
#'
#' @param filenames A character vector of filenames
#' @param output_file A string of the location of combined file
#' @param metadata_file A string of the location of the metadata file
#' @param card_dir A string of the CARD metadata directory that includes aro_index and aro_categories_index
#' @param subsampled A boolean stating whether reads were subsampled or not before mapping to CARD
#'
#' @examples
#'
#' @export
#'
#' @import dplyr purrr
combineBedtools <- function(filenames, output_file, metadata_file, card_dir, subsampled = FALSE) {

  # Function to read bedtools file
  readFile <- function(filename){
    if (file.info(filename)$size > 0){
      data <- read.delim(filename, header = F, stringsAsFactors = F)
      data$filename <- filename
      return(data)
    }
  }
  df_bedtools <- map_df(filenames, function(x) readFile(x))

  # Remove no mapping
  df_bedtools <- df_bedtools[df_bedtools$V7 != 0,]

  # Extract ARO for mapping
  df_bedtools$ARO <- gsub("\\|.*", "", gsub(".*ARO:", "", df_bedtools$V1))

  # Add IDs
  ID <- sapply(df_bedtools$filename, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) %>%
    sapply(function(x) gsub(".bam.bedtools.coverage.txt", "", x))
  if(subsampled){
    seq_num <- sapply(ID, function(x) strsplit(x, "-")[[1]][2]) %>% as.numeric()
    ID <- sapply(ID, function(x) strsplit(x, "-")[[1]][1])
    df_bedtools$seq_num <- seq_num
  }
  df_bedtools$ID <- ID

  # Read metadata
  metadata <- read.csv(metadata_file, stringsAsFactors = F)
  aro_index <- read.delim(paste0(card_dir, "aro_index.csv"), stringsAsFactors = F)
  aro_categories_index <- read.delim(paste0(card_dir, "aro_categories_index.csv"), stringsAsFactors = F)
  aro_index <- left_join(aro_index, aro_categories_index, by = "Protein.Accession")
  efflux_mech <- read.csv(paste0(card_dir, "efflux_mechanism.csv"), stringsAsFactors = FALSE)
  aro_index <- left_join(aro_index, efflux_mech, by = c("ARO.Accession", "ARO.Name"))

  # Combine raw and metadata
  df_bedtools <- left_join(df_bedtools, metadata, by = "ID")

  # Combine with aro_index
  ARO <- sapply(aro_index$ARO.Accession, function(x) gsub("ARO:", "", x))
  aro_index$ARO <- ARO
  df_bedtools <- left_join(df_bedtools, aro_index, by = "ARO")

  # Fill in the gaps
  if (sum(df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)) > 0) {
    df_bedtools$AMR.Gene.Family[df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)] <- "chloramphenicol acetyltransferase (CAT)"
    df_bedtools$Drug.Class[df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)] <- "phenicol antibiotic"
    df_bedtools$ARO.Name[df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)] <- "cat"
    df_bedtools$Resistance.Mechanism[df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)] <- "antibiotic inactivation"
  }

  # Write csv
  write.csv(df_bedtools, output_file, row.names = F)
}

#' Gets the best BLAST hit from overlapping ORFs
#'
#' @param blast_result A dataframe of filtered blast results
#' @param qseqid_id A character vector of the query sequence IDs to get the best BLAST hit
#' @examples
#'
#' @export
getBestBlastHit <- function(blast_result, qseqid_id){
  blast_result_tmp <- blast_result[blast_result$qseqid == qseqid_id,]
  remove_ids <- c()
  for(i in 1:(nrow(blast_result_tmp)-1)){

    for(j in (i+1):nrow(blast_result_tmp)){

      if(blast_result_tmp$qend[i] >= blast_result_tmp$qstart[j] & blast_result_tmp$qend[j] >= blast_result_tmp$qstart[i]){
        overlap_prop <- (blast_result_tmp$qend[i] - blast_result_tmp$qstart[j]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else if(blast_result_tmp$qend[j] >= blast_result_tmp$qstart[i] & blast_result_tmp$qstart[j] <= blast_result_tmp$qend[i]){
        overlap_prop <- (blast_result_tmp$qend[j] - blast_result_tmp$qstart[i]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else if(blast_result_tmp$qend[i] >= blast_result_tmp$qend[j] & blast_result_tmp$qstart[i] <= blast_result_tmp$qstart[j]){
        overlap_prop <- (blast_result_tmp$qend[i] - blast_result_tmp$qstart[i]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else if(blast_result_tmp$qend[j] >= blast_result_tmp$qend[i] & blast_result_tmp$qstart[j] <= blast_result_tmp$qstart[i]){
        overlap_prop <- (blast_result_tmp$qend[j] - blast_result_tmp$qstart[j]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else {
        overlap_prop <- 0
      }

      if(overlap_prop > 0.2 & blast_result_tmp[i, "evalue"] != blast_result_tmp[j, "evalue"]){
        ind <- which.max(blast_result_tmp$evalue[c(i,j)])
        remove_id <- which(blast_result$qseqid == blast_result_tmp$qseqid[c(i,j)[ind]] &
                             blast_result$sseqid == blast_result_tmp$sseqid[c(i,j)[ind]] &
                             blast_result$qstart == blast_result_tmp$qstart[c(i,j)[ind]] &
                             blast_result$qend == blast_result_tmp$qend[c(i,j)[ind]])
        remove_ids <- unique(c(remove_ids, remove_id))

      } else if (overlap_prop > 0.2 & blast_result_tmp[i, "evalue"] == blast_result_tmp[j, "evalue"] & blast_result_tmp[i, "pident"] != blast_result_tmp[j, "pident"]){
        ind <- which.min(blast_result_tmp$pident[c(i,j)])
        remove_id <- which(blast_result$qseqid == blast_result_tmp$qseqid[c(i,j)[ind]] &
                             blast_result$sseqid == blast_result_tmp$sseqid[c(i,j)[ind]] &
                             blast_result$qstart == blast_result_tmp$qstart[c(i,j)[ind]] &
                             blast_result$qend == blast_result_tmp$qend[c(i,j)[ind]])
        remove_ids <- unique(c(remove_ids, remove_id))
      }
    }
  }
  return(remove_ids)
}

#' Filters blast output file
#'
#' @param filename A string of the file name
#' @param column_names A character vector of the headers of the blast file
#' @param e_value_filter A numeric e-value cut-off
#' @param identity_filter A numeric identity_filter cut-off
#' @param min_length_filter A numeric minimum contig length cut-off
#'
#' @examples
#'
#' @export
#'
#' @import dplyr
filterBlastOutput <- function(filename, column_names, e_value_filter = NA, identity_filter = NA, min_length_filter = NA){
  df_tmp <- read.delim(file=filename, sep = "\t", stringsAsFactors = F, header = FALSE)

  names(df_tmp) <- column_names

  df_tmp$qseqid_mod <- gsub("_length.*", "", df_tmp$qseqid)

  if(!is.na(e_value_filter)){
    df_tmp <- df_tmp[df_tmp$evalue <= e_value_filter,]
  }
  if(!is.na(identity_filter)){
    df_tmp <- df_tmp[df_tmp$pident >= identity_filter,]
  }
  if(!is.na(min_length_filter)){
    df_tmp <- df_tmp[df_tmp$length >= min_length_filter,]
  }

  # Filter
  df_tmp <- df_tmp %>% group_by(qseqid, qstart, qend) %>%
    filter(evalue == min(evalue)) %>%
    filter(pident == max(pident))

  df_tmp <- df_tmp %>% group_by(qseqid, qstart) %>%
    filter(evalue == min(evalue)) %>%
    filter(pident == max(pident))

  df_tmp <- df_tmp %>% group_by(qseqid, qend) %>%
    filter(evalue == min(evalue)) %>%
    filter(pident == max(pident))

  # Extract ORFs with best hit
  df_tmp_dup_qseqids <- unique(df_tmp$qseqid[duplicated(df_tmp$qseqid)])
  remove_ids <- unlist(lapply(df_tmp_dup_qseqids, function(x) getBestBlastHit(df_tmp, x)))
  if(length(remove_ids) > 0){
    df_tmp <- df_tmp[-remove_ids,]
  }

  if(nrow(df_tmp) > 0){
    df_tmp$filename <- filename
    return(df_tmp)
  }
}

#' Combines CARD blast output files with metadata and CARD data
#'
#' @param filenames A character vector of bedtools filenames
#' @param card_dir A string of the CARD metadata directory that includes aro_index and aro_categories_index
#' @param metadata_file A string of the location of the metadata file
#' @param output_file A string of the location of the output
#' @examples
#'
#' @export
#'
#' @import dplyr purrr
combineCARDBlast <- function(filenames, card_dir, metadata_file, output_file){
  column_headers <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  card_comb <- map_df(.x = blast_card_filenames, .f = function(.x) filterBlastOutput(.x, column_names = column_headers, e_value_filter = 1e-50, identity_filter = 80))

  # CARD metadata
  aro_index <- read.delim(paste0(card_dir, "aro_index.csv"), stringsAsFactors = F)
  aro_categories_index <- read.delim(paste0(card_dir, "aro_categories_index.csv"), stringsAsFactors = F)
  card_metadata <- full_join(aro_index, aro_categories_index, by = "Protein.Accession")

  # Join blast hits with metadata
  card_comb$ARO.Accession <- paste0("ARO:", gsub("_.*", "", gsub(".*:", "", card_comb$sseqid, perl = TRUE)))
  card_metadata <- card_metadata[!duplicated(card_metadata$ARO.Accession),]
  card_comb <- left_join(card_comb, card_metadata, by = "ARO.Accession")

  # Add IDs
  ID <- sapply(card_comb$filename, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) %>%
    sapply(function(x) gsub("_card.out", "", x))
  card_comb$ID <- ID

  # Add and join metadata
  metadata <- read.csv(metadata_file, stringsAsFactors = F)
  card_comb <- left_join(card_comb, metadata, by = "ID")

  # Write file
  write.csv(card_comb, output_file, row.names = FALSE)
}

#' Combines PATRIC blast output files with metadata and PATRIC data
#'
#' @param filenames A character vector of bedtools filenames
#' @param patric_headers A dataframe of the faa headers in the patric database
#' @param metadata_file A string of the location of the metadata file
#' @param output_file A string of the location of the output
#' @examples
#'
#' @export
#'
#' @import dplyr purrr
combinePATRICBlast <- function(filenames, patric_headers, metadata_file, output_file){
  # Read patric metadata headers
  patric_headers$sseqid <- gsub("\\|", "_", gsub(">", "", gsub("[ ].*", "", patric_headers$V1)))
  patric_headers$putative <- !grepl("hypothetical protein", patric_headers$V1)
  patric_headers$phage <- gsub("\\ \\|.*", "", gsub(".*\\[", "", patric_headers$V1))

  # Combine patric blast results
  column_headers <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  blastp_comb <- map_df(.x = filenames, .f = function(.x) filterBlastOutput(.x, column_names = column_headers, e_value_filter = 1e-50, identity_filter = 80))
  blastp_comb <- left_join(blastp_comb, patric_headers, by = "sseqid")

  # Select putative proteins only
  blastp_comb <- blastp_comb[blastp_comb$putative,]

  # Add IDs
  ID <- sapply(blastp_comb$filename, function(x) strsplit(x, "/")[[1]][length(strsplit(x, "/")[[1]])]) %>%
    sapply(function(x) gsub(".out", "", x))
  blastp_comb$ID <- ID

  # Add and join metadata
  metadata <- read.csv(metadata_file, stringsAsFactors = F)
  blastp_comb <- left_join(blastp_comb, metadata, by = "ID")

  # Remove duplicate qseqids with same qstart, qend and phage
  blastp_comb <- blastp_comb[!duplicated(paste0(blastp_comb$qseqid, blastp_comb$qstart, blastp_comb$qend, blastp_comb$phage)),]

  # Write file
  write.csv(blastp_comb, output_file, row.names = FALSE)
}

#' Combines metaphlan2 output data
#'
#' @param filenames A character vector of metaphlan2 output filenames
#' @param output_file A string of the location of combined metaphlan2 file
#'
#' @examples
#'
#' @export
#'
#' @import dplyr
combineMetaphlanSamples <- function(filenames, output_file){
  metaphlan <- NULL
  for(i in 1:length(filenames)){
    df <- read.delim(filenames[i], stringsAsFactors = FALSE)
    sample <- gsub("_profile.txt", "", strsplit(filenames[i], "/")[[1]][4])
    names(df) <- c("Taxa", sample)
    if(is.null(metaphlan)){
      metaphlan <- df
    } else {
      metaphlan <- full_join(metaphlan, df, by = "Taxa")
    }
  }
  metaphlan[is.na(metaphlan)] <- 0
  row.names(metaphlan) <- metaphlan$Taxa
  metaphlan <- metaphlan[,-1]
  metaphlan <- metaphlan[rowSums(metaphlan == 0) != ncol(metaphlan),]

  write.csv(metaphlan, output_file)
}



