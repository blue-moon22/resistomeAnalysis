#' Read a csv file of combined mapping data and metadata
#'
#' @param filename A character vector of the csv
#' @param without_US_duplicates Logical: If TRUE, includes the US longitudinal data
#' @param coverage_threshold A numerical value were include only ARGs that have a coverage over the coverage threshold
#'
#' @return A dataframe of processed, combined mapping data and metadata
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_dup <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = FALSE)
#'
#' @export
#'
#' @import dplyr
readMappingData <- function(filename, without_US_duplicates = TRUE, coverage_threshold = 0.9){
  # Read mapping data
  df_map <- read.csv(filename, stringsAsFactors = FALSE)

  # Filter over by coverage
  df_map <- df_map[df_map$V7 > coverage_threshold,]

  # Create alternative drug class colum where more than 3 ARG classes are multidrug
  countCharOccurrences <- function(char, s) {
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))
  }
  drug_class <- as.character(df_map$Drug.Class)
  drug_class[as.logical(sapply(as.character(df_map$Drug.Class), function(x) countCharOccurrences(";", x)) > 2)] <- "multidrug"
  df_map$Drug.Class_mod <- drug_class

  # Create alternative sample type of oral and stool
  df_map$sample_type_mod <- ifelse(df_map$sample_type == "stool", "stool", "oral")

  # Remove hmp duplicates
  if(without_US_duplicates){
    df_map <- df_map[!(df_map$Visit_Number > 1 & df_map$Country == "US"),]
  }

  # Rename some of the ARGs in V1
  df_map$V1 <- gsub("ARO:3004039\\|Escherichia", "ARO:3004039\\|Escherichia_coli_emrE", df_map$V1)
  df_map$V1 <- gsub("ARO:3004043\\|Escherichia", "ARO:3004043\\|Escherichia_coli_acrA", df_map$V1)
  df_map$V1 <- gsub("ARO:3001328\\|Escherichia", "ARO:3001328\\|Escherichia_coli_mdfA", df_map$V1)
  df_map$V1 <- gsub("ARO:3004290\\|Escherichia", "ARO:3004290\\|Escherichia_coli_ampC_beta-lactamase", df_map$V1)
  df_map$V1 <- gsub("ARO:3003730\\|Bifidobacteria", "ARO:3003730\\|Bifidobacterium_ileS_conferring_resistance_to_mupirocin", df_map$V1)
  df_map$V1 <- gsub("ARO:3004041\\|Klebsiella", "ARO:3004041\\|Klebsiella_pneumoniae_acrA", df_map$V1)
  df_map$V1 <- gsub("ARO:3004122\\|Klebsiella", "ARO:3004122\\|Klebsiella_pneumoniae_OmpK37", df_map$V1)
  df_map$V1 <- gsub("ARO:3004454\\|Campylobacter", "ARO:3004454\\|Campylobacter_coli_chloramphenicol_acetyltransferase", df_map$V1)
  df_map$V1 <- gsub("ARO:3004042\\|Enterobacter", "ARO:3004042\\|Enterobacter_cloacae_acrA", df_map$V1)

  # Get total reads per sample and ARG richness
  total_reads <- df_map %>%
    group_by(ID, Cohort, sample_type) %>%
    summarise(total_reads_per_sample = sum(V4), arg_richness = n_distinct(ARO)) %>%
    ungroup()
  df_map <- left_join(df_map, total_reads, by = c("ID", "Cohort", "sample_type"))

  # Calculate RPKM
  rpkm <- (df_map$V4 * 10^3 * 10*6) / (df_map$total_reads_per_sample * df_map$V3)
  df_map$rpkm <- rpkm
  df_map$rpkm_log <- log10(df_map$rpkm)

  return(df_map)
}
