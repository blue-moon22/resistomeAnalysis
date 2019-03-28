#' Combines BLAST output data from CARD and PATRIC
#'
#' @param phages_data A dataframe containing combined phage blast data
#' @param args_data A dataframe containing combined ARG blast data
#'
#' @examples
#'
#' @export
#'
#' @import dplyr

combineARGsAndPhages <- function(phages_data, args_data, min_percentage) {

  # Combine phage and arg results
  args_data$qseqid_mod <- gsub("_length.*", "", args_data$qseqid)
  phages_args <- inner_join(phages_data, args_data, by = c("ID", "qseqid_mod"), suffix = c(".phage", ".arg"))

  # Proportion of samples that contain phage and ARG
  phages_args_pairs <- phages_args %>%
    group_by(Location.phage, sample_type.phage) %>%
    mutate(num_samples = n_distinct(ID)) %>%
    group_by(Location.phage, sample_type.phage, num_samples, phage, ARO.Name) %>%
    summarise(n_ID = n_distinct(ID)) %>%
    mutate(percentage = n_ID/num_samples * 100) %>%
    filter(!is.na(ARO.Name) & !is.na(phage)) %>%
    filter(percentage >= min_percentage)

  # Percentage samples that contain phage (for phages that have an ARG)
  perc_phages <- phages_data %>% group_by(Location, sample_type) %>%
    mutate(num_samples = n_distinct(ID)) %>%
    group_by(Location, sample_type, phage, num_samples) %>%
    summarise(num_samples_phages = n_distinct(ID)) %>%
    mutate(perc_samples_phages = num_samples_phages/num_samples*100)

  # Combine phage perc and phage-ARG perc
  phages_args_pairs <- left_join(phages_args_pairs, perc_phages, by = c("Location.phage" = "Location", "sample_type.phage" = "sample_type", "phage",
                                                                        "num_samples"))
  names(phages_args_pairs) <- c("Location", "sample_type", "num_samples", "phage", "arg_phage", "n_ID", "percentage", "num_samples_phage", "perc_samples_phages")
  perc_phage_split <- phages_args_pairs[,c("Location", "sample_type", "num_samples", "phage", "arg_phage", "num_samples_phage", "perc_samples_phages")]
  perc_phage_arg_split <- phages_args_pairs[,c("Location", "sample_type", "num_samples", "phage", "arg_phage", "n_ID", "percentage")]

  names(perc_phage_split) <- c("Location", "sample_type", "num_samples", "phage", "arg_phage", "n_ID", "percentage")
  perc_phage_split$arg_phage <- perc_phage_split$phage

  phages_args_pairs <- rbind(perc_phage_split, perc_phage_arg_split)

  # Add classes
  args_class <- args_data %>% select(ARO.Name, Drug.Class)
  args_class <- args_class[!(duplicated(args_class$ARO.Name)),]
  phages_args_pairs <- left_join(phages_args_pairs, args_class, by = c("arg_phage" = "ARO.Name"))
  phages_args_pairs$Drug.Class[is.na(phages_args_pairs$Drug.Class)] <- "Phage (No ARG Class)"
  phages_args_pairs$arg_phage <- factor(phages_args_pairs$arg_phage, levels = unique(phages_args_pairs$arg_phage))

  return(phages_args_pairs)
}

#' Plot bar charts of Phage-ARG combinations
#'
#' @param phage_args_comb A dataframe containing the percentage of samples containing each phage and ARG
#' @param phages_args_pairs A dataframe containing the percentage of samples containing phage-ARG pairs
#' @param sample_type A string of the sample type i.e. saliva
#' @param Location A string of the location i.e. "China"
#' @param colours A character vector of colours named by drug class
#'
#' @examples
#'
#' @export
#'
#' @import ggplot2
plotARGPhage <- function(phages_args_pairs, colours) {

  # Phages by sample type
  #phage_args_site <- unique(phages_args_pairs$phage[phages_args_pairs$sample_type %in% sample_type])
  #phage_args_comb_site <- phage_args_comb[phage_args_comb$phage %in% phage_args_site & phage_args_comb$sample_type %in% sample_type,]

  # Plot
  ggplot(phages_args_pairs, aes(arg_phage, percentage, fill = Drug.Class,
                                     label = arg_phage)) +
    geom_bar(stat="identity", position="dodge") +
    geom_text(stat = "summary", fun.y = max, hjust = -0.1, size = 2) +
    coord_flip() + xlab("Phage") +
    theme_classic() +
    scale_fill_manual("ARG Class", values = colours) +
    scale_y_continuous(limits = c(0, 125), breaks = seq(0, 125, 10)) +
    theme(axis.text.x = element_text(angle=70, hjust = 1),
          strip.background = element_blank(), strip.text.y = element_text(angle = 180),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

