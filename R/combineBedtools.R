#' Combine mapping data and metadata for all samples and write to a csv file
#'
#' @param filenames A character vector of filenames of mapping data from Bedtools
#' @param output_csv A string of the file path of the output csv file
#' @param metadata_csv A string of the file path of the metadata file
#' @param card_dir A string of the file directory of CARD metadata
#' @param subsampled Logical: if TRUE, the filenames will be split into IDs and the sequence number
#'
#' @return None
#'
#' @export
#'
#' @importFrom purrr map_df
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#' @importFrom graphics legend par plot
#' @importFrom stats cor.test cutree dist hclust model.matrix p.adjust quantile wilcox.test
#' @importFrom utils read.csv read.delim write.csv
#' @importFrom tools texi2pdf
combineBedtools <- function(filenames, output_csv, metadata_csv, card_dir, subsampled = FALSE){

  # Function to read bedtools file
  readFile <- function(filename){
    data <- read.delim(filename, header = FALSE, stringsAsFactors = FALSE)
    data$filename <- as.character(filename)
    return(data)
  }
  df_map_bedtools <- map_df(filenames, function(x) readFile(x))

  # Remove no mapping
  df_map_bedtools <- df_map_bedtools[df_map_bedtools$V7 != 0,]

  # Extract ARO for mapping and RGI
  ARO_mapping <- sapply(df_map_bedtools$V1, function(x) strsplit(x, "\\|")[[1]][5]) %>%
    sapply(function(x) strsplit(x, ":")[[1]][2])
  df_map_bedtools$ARO <- ARO_mapping

  # Add ID
  ID <- sapply(df_map_bedtools$filename, function(x) strsplit(x, "/")[[1]][4]) %>%
    sapply(function(x) gsub(".bam.bedtools.coverage.txt", "", x))

  # Add the number of sequences if subsampled reads
  if(subsampled){
    seq_num <- sapply(ID, function(x) strsplit(x, "-")[[1]][2]) %>% as.numeric()
    ID <- sapply(ID, function(x) strsplit(x, "-")[[1]][1])
    df_map_bedtools$sub_seqnum <- seq_num
  }

  df_map_bedtools$ID <- as.character(ID)

  # Read metadata
  metadata <- read.csv(metadata_csv, stringsAsFactors = FALSE)
  aro_index <- read.delim(paste0(card_dir, "aro_index.csv"), stringsAsFactors = FALSE)
  aro_categories_index <- read.delim(paste0(card_dir, "aro_categories_index.csv"), stringsAsFactors = FALSE)
  aro_index <- left_join(aro_index, aro_categories_index, by = "Protein.Accession")

  # Combine raw and metadata
  df_map_bedtools <- left_join(df_map_bedtools, metadata, by = "ID")

  # Combine with aro_index
  ARO <- sapply(aro_index$ARO.Accession, function(x) gsub("ARO:", "", x))
  aro_index$ARO <- ARO
  df_map_bedtools <- left_join(df_map_bedtools, aro_index, by = "ARO")

  # Remove rows with no metadata (can happen)
  df_map_bedtools <- df_map_bedtools[!is.na(df_map_bedtools$sample_type),]

  # Fill in the gaps
  df_map_bedtools$AMR.Gene.Family[df_map_bedtools$ARO == "3002670" & is.na(df_map_bedtools$Resistance.Mechanism)] <- "chloramphenicol acetyltransferase (CAT)"
  df_map_bedtools$Drug.Class[df_map_bedtools$ARO == "3002670" & is.na(df_map_bedtools$Resistance.Mechanism)] <- "phenicol antibiotic"
  df_map_bedtools$Resistance.Mechanism[df_map_bedtools$ARO == "3002670" & is.na(df_map_bedtools$Resistance.Mechanism)] <- "antibiotic inactivation"

  # Write csv
  write.csv(df_map_bedtools, output_csv, row.names = F)
}
