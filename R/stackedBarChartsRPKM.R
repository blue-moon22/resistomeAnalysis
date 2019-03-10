#' Get relative abundance of total rpkm per Location, sample type and drug class
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#'
#' @return A dataframe of relative abundances of drug classes for every Location and sample type
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_rel <- getRelativeAbundance(df_map)
#'
#' @export
#'
#' @import dplyr
getRelativeAbundance <- function(df_map){

  # Get relative rpkm abundance
  df_map_rel <- df_map %>%
    group_by(Location, Health, sample_type) %>%
    summarise(sum_rpkm = sum(rpkm)) %>%
    ungroup() %>%
    right_join(df_map) %>%
    mutate(rel_rpkm = rpkm/sum_rpkm) %>%
    group_by(Location, sample_type, Drug.Class) %>%
    summarise(sum_scale_rpkm = sum(rel_rpkm))

  # Extract most abundant ARG classes and label others as Other
  df_map_rel_top <- df_map_rel %>%
    group_by(Drug.Class) %>%
    summarise(total_sum_rel_rpkm = sum(sum_scale_rpkm)) %>%
    arrange(desc(total_sum_rel_rpkm)) %>%
    filter(total_sum_rel_rpkm > 0.5)
  Drug.Class.Alt <- df_map_rel$Drug.Class
  Drug.Class.Alt[!Drug.Class.Alt %in% df_map_rel_top$Drug.Class] <- "Other"
  df_map_rel$Drug.Class.Alt <- Drug.Class.Alt
  return(df_map_rel)
}

#' Plot stacked bar charts of the relative abundance of total rpkm per Location, sample type and drug class
#'
#' @param df_map_abundance A dataframe of relative abundances created from \code{\link{getRelativeAbundance}}
#' @param cols A character vector, named by ARG class
#'
#' @return A list of ggplot object where each graph represents a Location
#'
#' @examples
#' library(RColorBrewer)
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_rel <- getRelativeAbundance(df_map)
#' cols <- brewer.pal(length(unique(df_map_rel$Drug.Class.Alt)), "Set3")
#' names(cols) <- unique(df_map_rel$Drug.Class.Alt)
#' g <- plotARGClassAbundance(df_map_rel[df_map_rel$Location == "HMP1",], cols)
#'
#' @export
#'
#' @import ggplot2
plotARGClassAbundance <- function(df_map_abundance, cols){

  g <- ggplot(df_map_abundance, aes(sample_type, sum_scale_rpkm, fill = factor(Drug.Class.Alt))) +
    geom_bar(stat = "identity") +
    xlab("") + ylab("Relative Abundance") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 11)) +
    scale_fill_manual("ARG Class", values = cols) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))

  return(g)
}


#' Get relative abundance of total rpkm per individual, sample type and drug class
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#'
#' @return A dataframe of relative abundances of drug classes for every individual and sample type
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_rel_ind <- getRelativeAbundanceIndividuals(df_map)
#'
#' @export
#'
#' @import dplyr
getRelativeAbundanceIndividuals <- function(df_map){
  # Get top abundant classes
  df_map_rel_top <- df_map %>%
    group_by(Location, Health, sample_type) %>%
    summarise(sum_rpkm = sum(rpkm)) %>%
    ungroup() %>%
    right_join(df_map) %>%
    mutate(rel_rpkm = rpkm/sum_rpkm) %>%
    group_by(Location, sample_type, Drug.Class) %>%
    summarise(sum_scale_rpkm = sum(rel_rpkm)) %>%
    group_by(Drug.Class) %>%
    summarise(total_sum_rel_rpkm = sum(sum_scale_rpkm)) %>%
    arrange(desc(total_sum_rel_rpkm)) %>%
    filter(total_sum_rel_rpkm > 0.5)

  # Get individual relative abundance
  df_map_rel_ind <- df_map %>%
    group_by(ID, Location, Health, sample_type) %>%
    summarise(sum_rpkm = sum(rpkm)) %>%
    ungroup() %>%
    right_join(df_map) %>%
    mutate(rel_rpkm = rpkm/sum_rpkm) %>%
    group_by(ID, Location, sample_type, Drug.Class) %>%
    summarise(sum_scale_rpkm = sum(rel_rpkm))
  Drug.Class.Alt <- df_map_rel_ind$Drug.Class
  Drug.Class.Alt[!Drug.Class.Alt %in% df_map_rel_top$Drug.Class] <- "Other"
  df_map_rel_ind$Drug.Class.Alt <- Drug.Class.Alt

  return(df_map_rel_ind)

}

#' Plot stacked bar charts of the relative abundance of total rpkm per individual, sample type and drug class
#'
#' @param df_map_abundance A dataframe of relative abundances created from \code{\link{getRelativeAbundanceIndividuals}}
#' @param cols A character vector, named by ARG class
#'
#' @return A list of ggplot object where each graph represents a Location
#'
#' @examples
#' library(RColorBrewer)
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_rel_ind <- getRelativeAbundanceIndividuals(df_map)
#' cols <- brewer.pal(length(unique(df_map_rel_ind$Drug.Class.Alt)), "Set3")
#' names(cols) <- unique(df_map_rel_ind$Drug.Class.Alt)
#' g <- plotIndividualAbundance(df_map_rel_ind[df_map_rel_ind$Location == "HMP1",], cols)
#'
#' @export
#'
#' @import ggplot2
plotIndividualAbundance <- function(df_map_abundance, cols){

  g <- ggplot(df_map_abundance, aes(ID, sum_scale_rpkm, fill = factor(Drug.Class.Alt))) +
    geom_bar(stat = "identity") +
    xlab("Samples") + ylab("Relative Abundance") +
    scale_fill_manual(values = cols) +
    guides(fill=guide_legend(title="ARG Class")) +
    facet_grid(~ sample_type, scales = "free", space = "free") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank())
  return(g)
}

