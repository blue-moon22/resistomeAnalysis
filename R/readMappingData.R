#' Read a csv file of combined mapping data and metadata
#'
#' @param filename A character vector of the csv
#' @param without_US_duplicates Logical: If TRUE, includes the US longitudinal data
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
readMappingData <- function(filename, without_US_duplicates = TRUE){

  # Read mapping data
  df_map <- read.csv(filename, stringsAsFactors = FALSE)

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

  # Get total reads per sample and ARG richness
  df_map$Drug.Class.Efflux <- df_map$Drug.Class
  df_map$Drug.Class.Efflux[df_map$Resistance.Mechanism == "antibiotic efflux"] <- paste(df_map$Drug.Class.Efflux[df_map$Resistance.Mechanism == "antibiotic efflux"], "(efflux pump)")
  df_map$ARO.Name[!is.na(df_map$Multi.Efflux.Pump)] <- paste0(df_map$ARO.Name[!is.na(df_map$Multi.Efflux.Pump)], " [", df_map$Multi.Efflux.Role[!is.na(df_map$Multi.Efflux.Pump)], " ", df_map$Multi.Efflux.Pump[!is.na(df_map$Multi.Efflux.Pump)], "]")

  df_map$Prot[!is.na(df_map$Multi.Efflux.Pump)] <- df_map$Multi.Efflux.Pump[!is.na(df_map$Multi.Efflux.Pump)]
  total_reads <- df_map %>%
    group_by(ID, Cohort, sample_type) %>%
    summarise(total_reads_per_sample = sum(V4), arg_richness = n_distinct(ARO.Name)) %>%
    ungroup()
  df_map <- left_join(df_map, total_reads, by = c("ID", "Cohort", "sample_type"))

  # Calculate RPKM
  rpkm <- (df_map$V4 * 10^3 * 10^6) / (df_map$total_reads_per_sample * df_map$V3)
  df_map$rpkm <- rpkm
  df_map$rpkm_log <- log10(df_map$rpkm)

  return(df_map)
}
