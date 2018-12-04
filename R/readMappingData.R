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

  # Add country (Philippines okay)
  df_map$Country <- df_map$Cohort
  df_map$Country[df_map$Country == "HMP1"] <- "US"
  df_map$Country[df_map$Country %in% c("Twin", "Guys")] <- "UK"
  df_map$Country[df_map$Country == "Chinese"] <- "China"

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

  # Remove ARGs with N/A drug class (not used in CARD tool)
  df_map <- df_map[df_map$Drug.Class != "N/A",]

  if(without_US_duplicates){
    # Remove hmp duplicates
    df_map <- df_map[order(df_map$ID),]
    df_map_unique <- df_map[!duplicated(df_map$ID),]
    nondup_ids <- df_map_unique$ID[!duplicated(df_map_unique[,names(df_map_unique) %in% c("Sample.name", "sample_type")])]
    df_map <- df_map[df_map$ID %in% nondup_ids,]
  }

  # Rename some of the ARGs in V1
  new_V1 <- df_map$V1
  new_V1[new_V1 == "gb|Z11877.1|+|0-1244|ARO:3004039|Escherichia"] <- "gb|Z11877.1|+|0-1244|ARO:3004039|Escherichia_coli_emrE"
  new_V1[new_V1 == "gb|NC_000913.3|-|484425-485619|ARO:3004043|Escherichia"] <- "gb|NC_000913.3|-|484425-485619|ARO:3004043|Escherichia_coli_acrA"
  new_V1[new_V1 == "gb|JQ394987|+|0-1233|ARO:3001328|Escherichia"] <- "gb|JQ394987|+|0-1233|ARO:3001328|Escherichia_coli_mdfA"
  new_V1[new_V1 == "gb|NC_000913.3|-|4377810-4378944|ARO:3004290|Escherichia"] <- "gb|NC_000913.3|-|4377810-4378944|ARO:3004290|Escherichia_coli_ampC_beta-lactamase"
  new_V1[new_V1 == "gb|NC_014638|-|1610636-1613960|ARO:3003730|Bifidobacteria"] <- "gb|NC_014638|-|1610636-1613960|ARO:3003730|Bifidobacterium_ileS_conferring_resistance_to_mupirocin"
  new_V1[new_V1 == "gb|AJ318073.1|+|793-1990|ARO:3004041|Klebsiella"] <- "gb|NC_014638|-|1610636-1613960|ARO:3003730|Klebsiella_pneumoniae_acrA"
  new_V1[new_V1 == "gb|AJ011502|+|0-1463|ARO:3004122|Klebsiella"] <- "gb|NC_014638|-|1610636-1613960|ARO:3003730|Klebsiella_pneumoniae_OmpK37"
  df_map$V1 <- new_V1

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
