#' Creates a dataframe that shows the proportion of samples that contain a unique property for every unique property for each country
#' i.e. the proportion of samples that contain an ARG class for all ARG classes found in column Drug.Class for each country
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param level A string of a column name e.g. Drug.Class
#'
#' @return A dataframe of proportions of samples that contain a unique property for each country
#'
#' @export
#'
#' @import dplyr
#' @importFrom reshape2 melt
getProportion <- function(df_map, level){

  # Number of samples per country
  sample_count <- df_map %>%
    group_by(Country, ID, Health) %>%
    summarise(sample_n = n())

  # Number of samples in a country
  country_count <- sample_count %>%
    group_by(Country, Health) %>%
    summarise(country_n = n())
  country_count$country_health <- paste(country_count$Country, country_count$Health, sep="-")

  # Get unique AMR levels e.g. class or mechanism
  level_names <- unique(unlist(sapply(df_map[[level]], function(x) strsplit(x, ";")[[1]])))
  level_names <- level_names[level_names != "N/A"] # Remove N/A level name e.g. class name
  presence_matrix <- matrix(0, nrow = length(level_names), ncol = nrow(country_count))

  # For each class name, find proportion for each country
  for(i in 1:length(level_names)){
    if(level %in% c("Drug.Class", "Resistance.Mechanism")){
      tmp <- df_map[grep(level_names[i], df_map[[level]]), c("ID", "Country", "Health")]
    } else {
      tmp <- df_map[df_map[[level]] == level_names[i], c("ID", "Country", "Health")]
    }
    tmp <- tmp[!duplicated(tmp$ID),]
    tmp$country_health <- paste(tmp$Country, tmp$Health, sep="-")
    for(j in 1:nrow(country_count)){
      presence_matrix[i,j] <- sum(tmp$country_health == country_count$country_health[j]) / (country_count$country_n[country_count$country_health == country_count$country_health[j]])
    }
  }

  # Clean presence matrix to df_map and melt
  presence_matrix <- presence_matrix*100
  presence_df_map <- as.data.frame(presence_matrix)
  names(presence_df_map) <- country_count$country_health
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
  for(i in 1:B){
    sample_n <- sample(1:n, n, replace=T)
    boot_sample <- df_map[df_map$ID %in% individuals[sample_n],]
    presence_names <- map_dbl(names, function(x) search_name(boot_sample, x, level))
    presence_names <- cbind(presence_names, names)
    boot_perc[(i*length(names) - length(names) + 1):(i*length(names)), 1:2] <- as.matrix(presence_names)
  }
  boot_perc <- data.frame(boot_perc, stringsAsFactors = F)
  names(boot_perc) <- c("value", "level_names")
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
    group_by(level_names) %>%
    summarise(CI_lb=quantile(value, probs = 0.025),
              CI_ub=quantile(value, probs = 0.975))
  boot_ci$level_names <- as.character(boot_ci$level_names)
  return(boot_ci)
}

#' Merges dataframes of proportion of samples containing a certain property and confidence intervals of proportions from bootstrapping samples
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param level A string of a column name e.g. "Drug.Class"
#' @param sample_type A string of a sample type e.g. "saliva"
#' @param B A numerical value of the number of types to sample in bootstrapping
#'
#' @return A dataframe of proportions and 95% confidence intervals of samples that contain a unique property for each country
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' df_map_pb <- joinProportionAndBootstrap(df_map, "Drug.Class", "saliva")
#'
#' @export
#'
#' @importFrom dplyr left_join
joinProportionAndBootstrap <- function(df_map, level = "Drug.Class", sample_type = "saliva", B = 100){

  # Filter sample type
  df_map <- df_map[df_map$sample_type %in% sample_type,]

  # Get proportion of samples
  prop <- getProportion(df_map, level = level)
  names(prop) <- c("level", "country_health", "proportion")

  # Bootstrap through countries and health states
  df_map$country_health <- paste(df_map$Country, df_map$Health, sep="-")
  uniq_country_health <- unique(df_map$country_health)
  boot_ci <- list()
  for(i in 1:length(uniq_country_health)){
    df_map_filt <- df_map[df_map$country_health %in% uniq_country_health[i],]
    boot_perc <- bootstrap_percentage(df_map_filt, level, B=B)
    boot_ci[[i]] <- data.frame(bootstrap_CI(boot_perc), country_health=uniq_country_health[i])
  }

  # Combine CIs
  boot_ci_comb <- do.call(rbind, boot_ci)
  boot_names <- names(boot_ci_comb)
  boot_names[1] <- "level"
  names(boot_ci_comb) <- boot_names

  # Join proportion and bootstrap CIs
  boot_ci_comb$country_health <- as.character(boot_ci_comb$country_health)
  prop$country_health <- as.character(prop$country_health)
  df_map_join <- left_join(prop, boot_ci_comb, by = c("country_health", "level"))

  # Get the 95% confidence interval
  df_map_join$CI_lb95 <- df_map_join$CI_lb + (df_map_join$proportion - df_map_join$CI_lb) * 0.05
  df_map_join$CI_ub95 <- df_map_join$CI_ub - (df_map_join$CI_ub - df_map_join$proportion) * 0.05

  df_map_join$Country <- sapply(df_map_join$country_health, function(x) strsplit(x, "-")[[1]][1])

  return(df_map_join)
}

#' Plots a bar chart of the porportions of samples that contain a property for each country
#'
#' @param df_map_pb A dataframe of proportions and 95\% confidence intervals of samples that contain a unique property for each country
#' @param colours A character vector of colours, named by country
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
plotPercentages <- function(df_map_pb, colours){
  ggplot(df_map_pb, aes(level, proportion, ymin = CI_lb95, ymax = CI_ub95, fill = Country)) +
    geom_bar(stat = "identity", position = "dodge", colour="black") +
    geom_errorbar(position = position_dodge()) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 11),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    guides(fill=guide_legend(title="Country")) +
    scale_fill_manual(values = colours) +
    scale_y_continuous(limits = c(0,100), expand = c(0,0), breaks = seq(0, 100, 10))
}

#' Combines proportions and confidence intervals from non-subsampled and equivalent subsampled samples
#'
#' @param df_map A dataframe of combined non-subsampled mapping data and metadata
#' @param df_map_sub A dataframe of combined subsampled mapping data and metadata
#' @param level A string of a column name e.g. "Drug.Class"
#' @param sample_type A string of a sample type e.g. "saliva"
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
combineSubsampled <- function(df_map, df_map_sub, level, sample_type){

  # Joined proportion and bootstrap
  df_map_pb <- joinProportionAndBootstrap(df_map, level, sample_type)
  df_map_pb_sub <- joinProportionAndBootstrap(df_map_sub, level, sample_type)

  # Combine non-subsampled and subsampled
  df_map_pb_comb <- rbind(data.frame(df_map_pb, subsampled = "non-subsampled"), data.frame(df_map_pb_sub, subsampled = "subsampled"))
  return(df_map_pb_comb)
}

#' Plots a scatter plots of the porportions of samples that contain a property for each country for subsampled and non-subsampled samples
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
  df_map_pb_comb$class_country_health <- paste(df_map_pb_comb$class, df_map_pb_comb$country_health, sep = ": ")
  df_map_pb_comb$country <- sapply(as.character(df_map_pb_comb$country_health), function(x) strsplit(x, "-")[[1]][1])
  df_map_pb_comb$country_health_subsampled <- paste(df_map_pb_comb$country_health, df_map_pb_comb$subsampled, sep = ": ")
  df_map_pb_comb$country_subsampled <- paste(df_map_pb_comb$country, df_map_pb_comb$subsampled, sep = ": ")

  # Choose colours
  cols <- brewer.pal(length(unique(df_map_pb_comb$country_health_subsampled)), "Paired")
  for(i in 1:(length(cols)/2)){
    cols[(2*i-1):(2*i)] <- cols[(2*i):(2*i-1)]
  }

  ggplot(df_map_pb_comb, aes(class_country_health, proportion, col = country_subsampled)) +
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
    scale_color_manual("Country: non-subsampled/subsampled", values = cols) +
    scale_y_continuous(breaks = seq(0, 100, 10))
}

