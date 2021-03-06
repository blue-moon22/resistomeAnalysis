#' Creates a dataframe that shows the proportion of samples that contain a unique property for every unique property for each Location
#' i.e. the proportion of samples that contain an ARG class for all ARG classes found in column Drug.Class for each Location
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param level A string of a column name e.g. Drug.Class
#'
#' @return A dataframe of proportions of samples that contain a unique property for each Location
#'
#' @export
#'
#' @import dplyr
#' @importFrom reshape2 melt
getProportion <- function(df_map, level){

  # Number of samples per Location
  sample_count <- df_map %>%
    group_by(Location, ID, Health) %>%
    summarise(sample_n = n())

  # Number of samples in a Location
  Location_count <- sample_count %>%
    group_by(Location, Health) %>%
    summarise(Location_n = n())
  Location_count$Location_health <- paste(Location_count$Location, Location_count$Health, sep="-")

  # Get unique AMR levels e.g. class or mechanism
  if(level %in% c("Drug.Class", "Resistance.Mechanism")){
    level_names <- unique(unlist(sapply(df_map[[level]], function(x) strsplit(x, ";")[[1]])))
  } else {
    level_names <- unique(df_map[[level]])
  }
  presence_matrix <- matrix(0, nrow = length(level_names), ncol = nrow(Location_count))

  # For each class name, find proportion for each Location
  for(i in 1:length(level_names)){
    if(level %in% c("Drug.Class", "Resistance.Mechanism")){
      tmp <- df_map[grep(level_names[i], df_map[[level]]), c("ID", "Location", "Health")]
    } else {
      tmp <- df_map[df_map[[level]] == level_names[i], c("ID", "Location", "Health")]
    }
    tmp <- tmp[!duplicated(tmp$ID),]
    tmp$Location_health <- paste(tmp$Location, tmp$Health, sep="-")
    for(j in 1:nrow(Location_count)){
      presence_matrix[i,j] <- sum(tmp$Location_health == Location_count$Location_health[j]) / (Location_count$Location_n[Location_count$Location_health == Location_count$Location_health[j]])
    }
  }

  # Clean presence matrix to df_map and melt
  presence_matrix <- presence_matrix*100
  presence_df_map <- as.data.frame(presence_matrix)
  names(presence_df_map) <- Location_count$Location_health
  rownames(presence_df_map) <- level_names
  presence_df_map <- presence_df_map[rowSums(presence_df_map) > 0,] # Remove no AMR classes
  presence_df_map$level_names <- rownames(presence_df_map)

  presence_df_map_m <- melt(presence_df_map)

  return(presence_df_map_m)
}

#' Bootstraps the samples and recalculates proportion of samples containing a property
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param level A string of a column name e.g. Drug.Class
#' @param B A numerical value of the number of types to sample in bootstrapping
#'
#' @return A dataframe of proportions of samples for each bootstrap iteration
#'
#' @export
#'
#' @importFrom purrr map_dbl
#' @importFrom reshape2 melt
bootstrap_percentage <- function(df_map, level="Drug.Class", B){

  # Unique levels (eg unique drug classes)
  names <- unique(unlist(sapply(df_map[[level]], function(x) strsplit(x, ";")[[1]])))

  # Function to search for name i.e. class name
  search_name <- function(bootstrap_data, name, level){
    if (level %in% c("Drug.Class", "Resistance.Mechanism")){
      match_name <- bootstrap_data$ID[grep(name, bootstrap_data[[level]])]
    } else {
      match_name <- bootstrap_data$ID[bootstrap_data[[level]] == name]
    }
    presence_name <- (length(unique(match_name)))/(length(unique(bootstrap_data$ID)))*100
    return(presence_name)
  }

  boot_perc <- matrix(NA, B*length(names), 2)
  individuals <- unique(df_map$ID)
  n <- length(individuals)
  set.seed(42)
  for(i in 1:B){
    sample_n <- sample(1:n, n, replace=T)
    boot_sample <- df_map[df_map$ID %in% individuals[sample_n],]
    presence_names <- map_dbl(names, function(x) search_name(boot_sample, x, level))
    presence_names <- cbind(presence_names, names)
    boot_perc[(i*length(names) - length(names) + 1):(i*length(names)), 1:2] <- as.matrix(presence_names)
  }
  boot_perc <- data.frame(boot_perc, stringsAsFactors = F)
  names(boot_perc) <- c("value", "level")
  boot_perc$value <- as.numeric(boot_perc$value)
  return(boot_perc)
}

#' Calculates the 95\% confidence intervals for each bootstrap iteration
#'
#' @param bootstrap_percentage_output The dataframe output from bootstrap_percentage function
#'
#' @return A dataframe of 95\% confidence intervals
#'
#' @export
#'
#' @import dplyr
bootstrap_CI <- function(bootstrap_percentage_output){
  boot_ci <- bootstrap_percentage_output %>%
    group_by(level) %>%
    summarise(proportion = mean(value),
              CI_lb=quantile(value, probs = 0.025),
              CI_ub=quantile(value, probs = 0.975))
  boot_ci$level <- as.character(boot_ci$level)
  return(boot_ci)
}

#' Merges dataframes of proportion of samples containing a certain property and confidence intervals of proportions from bootstrapping samples
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param level A string of a column name e.g. "Drug.Class"
#' @param B A numerical value of the number of types to sample in bootstrapping
#'
#' @return A dataframe of proportions and 95% confidence intervals of samples that contain a unique property for each Location
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_pb <- joinProportionAndBootstrap(df_map, "Drug.Class", "saliva")
#'
#' @export
#'
#' @importFrom dplyr left_join
joinProportionAndBootstrap <- function(df_map, level = "Drug.Class", B = 100){

  # # Get proportion of samples
  # prop <- getProportion(df_map, level = level)
  # names(prop) <- c("level", "Location_health", "proportion")

  # Bootstrap through countries and health states
  df_map$Location_health <- paste(df_map$Location, df_map$Health, sep="-")
  uniq_Location_health <- unique(df_map$Location_health)
  boot_ci <- data.frame()
  boot_perc <- data.frame()
  for(i in 1:length(uniq_Location_health)){
    df_map_filt <- df_map[df_map$Location_health %in% uniq_Location_health[i],]
    boot_perc_tmp <- bootstrap_percentage(df_map_filt, level, B=B)
    boot_ci_tmp <- bootstrap_CI(boot_perc_tmp)
    boot_perc_tmp$Location_health = uniq_Location_health[i]
    boot_ci_tmp$Location_health = uniq_Location_health[i]
    boot_perc <- rbind(boot_perc, boot_perc_tmp)
    boot_ci <- rbind(boot_ci, boot_ci_tmp)
  }

  # Join proportion and bootstrap CIs
  boot_ci$Location <- gsub("-.*", "", boot_ci$Location_health)
  boot_perc$Location <- gsub("-.*", "", boot_perc$Location_health)

  return(list(boot_perc = boot_perc, boot_ci = boot_ci))
}

#' Plots a bar chart of the porportions of samples that contain a property for each Location
#'
#' @param df_map_pb A dataframe of proportions and 95\% confidence intervals of samples that contain a unique property for each Location
#' @param colours A character vector of colours, named by Location
#'
#' @return None
#'
#' @examples
#' library(RColorBrewer)
#' library(ggplot2)
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' cols <- c("grey20", "grey50", "grey80", "white")
#' names(cols) <- c("China", "Philippines", "UK", "US")
#' df_map_pb_saliva_class <- joinProportionAndBootstrap(df_map, "Drug.Class", "saliva")
#'
#' plotPercentages(df_map_pb_saliva_class, cols) + xlab("ARG class") + ylab("% saliva samples containing ARG class")
#'
#' @export
#'
#' @import ggplot2
plotPercentages <- function(df_map_pb, cols){
  ggplot() +
    geom_bar(data = df_map_pb$boot_ci, mapping = aes(level, proportion, fill = Location),
             position = "dodge", colour="black", stat = "identity") +
    geom_point(data = df_map_pb$boot_perc, mapping = aes(level, value, group = Location), shape = 4, position=position_dodge(0.9)) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 11),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    guides(fill=guide_legend(title="Location")) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(limits = c(0,100), expand = c(0,0), breaks = seq(0, 100, 10)) +
    geom_errorbar(data = df_map_pb$boot_ci, mapping = aes(x = level, ymin = CI_lb, ymax = CI_ub, group = Location),
                  position = position_dodge()) +
    theme(axis.title=element_text(size=18), legend.text = element_text(size=14), legend.title = element_text(size=18))
}

#' Combines proportions and confidence intervals from non-subsampled and equivalent subsampled samples
#'
#' @param df_map A dataframe of combined non-subsampled mapping data and metadata
#' @param df_map_sub A dataframe of combined subsampled mapping data and metadata
#' @param level A string of a column name e.g. "Drug.Class"
#'
#' @return A dataframe with proportions and confidence intervals from non-subsampled and equivalent subsampled samples
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_sub_saliva <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/subsampled_saliva_merged.csv", without_US_duplicates = TRUE)
#'
#' df_map_pb_comb <- combineSubsampled(df_map, df_map_sub_saliva, "Drug.Class", "saliva")
#'
#' @export
combineSubsampled <- function(df_map, df_map_sub, level){

  # Joined proportion and bootstrap
  df_map_pb <- joinProportionAndBootstrap(df_map, level)
  df_map_pb_sub <- joinProportionAndBootstrap(df_map_sub, level)

  # Combine non-subsampled and subsampled
  df_map_pb_comb <- rbind(data.frame(df_map_pb, subsampled = "non-subsampled"), data.frame(df_map_pb_sub, subsampled = "subsampled"))
  return(df_map_pb_comb)
}

#' Plots a scatter plots of the porportions of samples that contain a property for each Location for subsampled and non-subsampled samples
#'
#' @param df_map_pb_comb A dataframe with proportions and confidence intervals from non-subsampled and equivalent subsampled samples
#'
#' @return None
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_sub_saliva <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/subsampled_saliva_merged.csv", without_US_duplicates = TRUE)
#' df_map_pb_comb <- combineSubsampled(df_map, df_map_sub_saliva, "Drug.Class", "saliva")
#'
#' plotScatter(df_map_pb_comb)
#'
#' @export
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
plotScatter <- function(df_map_pb_comb){
  # Remove 0 values
  df_map_pb_comb <- df_map_pb_comb[df_map_pb_comb$proportion != 0,]

  # Add labels
  df_map_pb_comb$class_Location_health <- paste(df_map_pb_comb$class, df_map_pb_comb$Location_health, sep = ": ")
  df_map_pb_comb$Location <- sapply(as.character(df_map_pb_comb$Location_health), function(x) strsplit(x, "-")[[1]][1])
  df_map_pb_comb$Location_health_subsampled <- paste(df_map_pb_comb$Location_health, df_map_pb_comb$subsampled, sep = ": ")
  df_map_pb_comb$Location_subsampled <- paste(df_map_pb_comb$Location, df_map_pb_comb$subsampled, sep = ": ")

  # Choose colours
  cols <- brewer.pal(length(unique(df_map_pb_comb$Location_health_subsampled)), "Paired")
  for(i in 1:(length(cols)/2)){
    cols[(2*i-1):(2*i)] <- cols[(2*i):(2*i-1)]
  }

  ggplot(df_map_pb_comb, aes(class_Location_health, proportion, col = Location_subsampled)) +
    geom_errorbar(aes(ymin=CI_lb95, ymax=CI_ub95), width=.7) +
    geom_point() +
    ylab("% samples") +
    facet_grid( ~ level, scales = "free", space = "free", switch = "x") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.text.x = element_text(angle = 90),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank()) +
    scale_color_manual("Location: non-subsampled/subsampled", values = cols) +
    scale_y_continuous(breaks = seq(0, 100, 10))
}

