#' Creates a dataframe of paired sample types with a column named group
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param pair_list A list of vectors of length two containing paired sample types
#'
#' @return A dataframe of paired sample types with a column named group i.e. saliva vs. stool
#'
#' @examples
#' df_map_subsampled <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/subsampled_argrich_merged.csv", without_US_duplicates = TRUE)
#' pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"), c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"), c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
#' df_map_subsampled_pairs <- createPairedData(df_map_subsampled, pair_list)
#'
#' @export
createPairedData <- function(df_map, pair_list){

  # Create groups e.g. stool vs dental
  df_map_pairs <- data.frame()
  for(i in 1:length(pair_list)){
    samples_one <- unique(df_map$Sample.name[df_map$sample_type == pair_list[[i]][1]])
    samples_two <- unique(df_map$Sample.name[df_map$sample_type == pair_list[[i]][2]])
    df_map_pair <- df_map[(df_map$Sample.name %in% Reduce(intersect,list(samples_one, samples_two))) & (df_map$sample_type %in% pair_list[[i]]),]

    df_map_pair <- cbind(df_map_pair, group = paste(pair_list[[i]][1], "vs.", pair_list[[i]][2]))
    df_map_pairs <- rbind(df_map_pairs, df_map_pair)
  }

  # Order by Location and change characters in group
  df_map_pairs <- df_map_pairs[order(df_map_pairs$Location),]
  df_map_pairs$group <- as.character(df_map_pairs$group)
  df_map_pairs$group_mod <- gsub("vs.", "\nvs.", df_map_pairs$group)

  return(df_map_pairs)
}

#' Plots a box and whisker diagram between paired samples types
#'
#' @param df_map_pairs_group A dataframe representing one Location with paired samples
#'
#' @return A ggplot object
#'
#' @export
#'
#' @import ggplot2
plotRPKM <- function(df_map_pairs_group){

  g <- ggplot(df_map_pairs_group, aes(sample_type, rpkm_log)) +
    geom_boxplot() +
    ggtitle(unique(df_map_pairs_group$Location)) +
    theme_classic() +
    theme(strip.text.x = element_blank(),
          legend.position = "none") +
    ylab("Abundance (log10[RPKM])") +
    xlab("")
  return(g)
}

#' Plots multiple box and whisker diagrams of total RPKM between paired sample types for every Location
#'
#' @param df_map_pairs A dataframe of paired sample types with a column named group i.e. saliva vs. stool
#'
#' @return A list of box and whisker ggplot objects.
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"), c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"), c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
#' df_map_pairs <- createPairedData(df_map, pair_list)
#'
#' g <- plotMultipleRPKM(df_map_pairs)
#'
#' @export
#'
#' @import ggplot2
plotMultipleRPKM <- function(df_map_pairs){

  graphs <- list()
  count <- 0
  unique_Location <- unique(df_map_pairs$Location)
  for(i in 1:length(unique_Location)){
    tmp <- df_map_pairs[df_map_pairs$Location == unique_Location[i],]
    unique_groups <- unique(tmp$group)
    for(j in 1:length(unique_groups)){
      count <- count + 1
      df_map_pairs_groups <- tmp[tmp$group == unique_groups[j],]
      graphs[[count]] <- plotRPKM(df_map_pairs_groups)
    }
  }
  return(graphs)
}

