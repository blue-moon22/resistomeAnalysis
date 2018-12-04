#' Plot scatter graph of proportion against DDD per 1000
#'
#' @param data_merge A dataframe of proportions containing ARG classes and confidence intervals, and resistanceMap DDD per 1000
#' @param country A string of a country e.g. "US"
#' @param class_cols A character vector of colours, named by ARG class
#'
#' @return None
#'
#' @examples
#' library(RColorBrewer)
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_pb_saliva_class <- joinProportionAndBootstrap(df_map, "Drug.Class", "saliva")
#' df_map_pb_dental_class <- joinProportionAndBootstrap(df_map, "Drug.Class", "dental")
#' rm_use_card <- readResistanceMapUseData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/ResistanceMap/resistanceMap_use.csv")
#' data_merge <- mergeData(list(saliva = df_map_pb_saliva_class, dental = df_map_pb_dental_class), rm_use_card)
#'
#' cols <- brewer.pal(length(unique(data_merge$class)), "Paired")
#' names(cols) <- unique(data_merge$class)
#' plotGraph(data_merge, "US", cols)
#'
#' @export
#'
#' @import ggplot2
plotGraph <- function(data_merge, country, class_cols){
  shape_values <- c(1,4)
  names(shape_values) <- unique(data_merge$sample_type)

  # Plot graph
  data_merge_country <- data_merge[data_merge$Country %in% country,]
  ggplot(data_merge_country, aes(sum.DDD.Per.1000.Pop, proportion,
                           ymin = CI_lb95, ymax = CI_ub95,
                           color = factor(class), shape = factor(sample_type))) +
    geom_point(size = 2) +
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    geom_errorbar(aes(ymin=CI_lb95, ymax=CI_ub95), width = max(data_merge_country$sum.DDD.Per.1000.Pop)/100) +
    xlab("Defined Daily Doses Per 1000 individuals") + ylab("% samples") +
    scale_color_manual("ARG Class", values = cols) +
    scale_shape_manual("Oral type", values = shape_values) +
    ggtitle(country)
}

#' Read resistanceMap use data
#'
#' @param filename A string of the csv filename of resistanceMap use data
#'
#' @return None
#'
#' @examples
#' rm_use_card <- readResistanceMapUseData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/ResistanceMap/resistanceMap_use.csv")
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @import dplyr
readResistanceMapUseData <- function(filename){
  # Read data
  rm_use <- read.csv(filename, stringsAsFactors = FALSE)

  # Sum DDD by CARD class
  rm_use_card <- rm_use %>% group_by(Country, CARD.Class) %>%
    summarise(sum.DDD.Per.1000.Pop = sum(DDD.Per.1000.Pop)) %>%
    dplyr::rename(class = CARD.Class) %>%
    as.data.frame()

  return(rm_use_card)
}

#' Merge dataframes proportion and confidence intervals of samples containing ARG class and resistanceMap use data
#'
#' @param percentages_data_list A list of dataframes of proportions containing ARG classes and confidence intervals for each sample type
#' @param resistance_map_data A dataframe containing resistanceMap use data
#'
#' @return A dataframe of proportions containing ARG classes and confidence intervals, and resistanceMap DDD per 1000
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_pb_saliva_class <- joinProportionAndBootstrap(df_map, "Drug.Class", "saliva")
#' df_map_pb_dental_class <- joinProportionAndBootstrap(df_map, "Drug.Class", "dental")
#' rm_use_card <- readResistanceMapUseData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/ResistanceMap/resistanceMap_use.csv")
#'
#' data_merge <- mergeData(list(saliva = df_map_pb_saliva_class, dental = df_map_pb_dental_class), rm_use_card)
#'
#' @export
#'
#' @importFrom dplyr full_join
mergeData <- function(percentages_data_list, resistance_map_data){
  # Merge all samples types percentages data
  for(i in 1:length(percentages_data_list)){
    percentages_data_list[[i]] <- percentages_data_list[[i]][percentages_data_list[[i]]$proportion > 0,]
    percentages_data_list[[i]]$sample_type <- names(percentages_data_list)[i]
  }
  perc_data_all <- do.call("rbind", percentages_data_list)
  names(perc_data_all)[names(perc_data_all) == "level"] <- "class"

  # Join metagenomic percentages data with resistance map data
  data_merge <- full_join(resistance_map_data, perc_data_all, by = c("Country", "class"))
  data_merge <- data_merge[!is.na(data_merge$proportion),]
  data_merge <- data_merge[!is.na(data_merge$sum.DDD.Per.1000.Pop),]
  return(data_merge)
}

#' Extract legend from ggplot object
#'
#' @param a.gplot A ggplot object
#'
#' @return A ggplot object legend
#'
#' @export
#'
#' @import ggplot2
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
