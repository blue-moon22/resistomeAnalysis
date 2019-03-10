#' Combines bedtools mapping data with metadata and CARD data
#'
#' @param filenames A character vector of bedtools filenames
#' @param output_file A string of the location of combined bedtools and metadata file
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
    data <- read.delim(filename, header = F, stringsAsFactors = F)
    data$filename <- filename
    return(data)
  }
  df_bedtools <- map_df(filenames, function(x) readFile(x))

  # Remove no mapping
  df_bedtools <- df_bedtools[df_bedtools$V7 != 0,]

  # Extract ARO for mapping and RGI
  df_bedtools$ARO <- gsub("\\|.*", "", gsub(".*ARO:", "", df_bedtools$V1))

  # Add IDs
  ID <- sapply(df_bedtools$filename, function(x) strsplit(x, "/")[[1]][4]) %>%
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

  # Combine raw and metadata
  df_bedtools <- left_join(df_bedtools, metadata, by = "ID")

  # Combine with aro_index
  ARO <- sapply(aro_index$ARO.Accession, function(x) gsub("ARO:", "", x))
  aro_index$ARO <- ARO
  df_bedtools <- left_join(df_bedtools, aro_index, by = "ARO")

  # Fill in the gaps
  if (sum(df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)) > 0) {
    df_bedtools$AMR.Gene.Family[df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)] <- "chloramphenicol acetyltransferase (CAT)"
  }
  if (sum(df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)) > 0) {
    df_bedtools$Drug.Class[df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)] <- "phenicol antibiotic"
  }
  if (sum(df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)) > 0) {
    df_bedtools$Resistance.Mechanism[df_bedtools$ARO == "3002670" & is.na(df_bedtools$Resistance.Mechanism)] <- "antibiotic inactivation"
  }

  # Write csv
  write.csv(df_bedtools, output_file, row.names = F)
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



