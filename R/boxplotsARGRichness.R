#' Runs a t-test between groups from a dataframe of paired samples
#'
#' @param df_map_subsampled_pairs A dataframe of paired samples containing a column named "group" i.e. saliva vs. stool
#'
#' @return A dataframe of p-values, asterisks, countries and groups.
#'
#' @examples
#' df_map_subsampled <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/subsampled_argrich_merged.csv", without_US_duplicates = TRUE)
#' pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"), c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"), c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
#' df_map_subsampled_pairs <- createPairedData(df_map_subsampled, pair_list)
#'
#' ttest_groups <- runTtest(df_map_subsampled_pairs)
#'
#' @export
runTtest <- function(df_map_subsampled_pairs){

  df_map_argrich <- df_map_subsampled_pairs[!duplicated(paste0(df_map_subsampled_pairs$ID, df_map_subsampled_pairs$group)),]

  # T-test
  p_values <- c()
  ttest_groups <- df_map_argrich[!duplicated(paste0(df_map_argrich$Location, df_map_argrich$group)),]
  for(i in 1:nrow(ttest_groups)){
    y <- df_map_argrich[df_map_argrich$Location == ttest_groups$Location[i] & df_map_argrich$group == ttest_groups$group[i],]
    y <- y[order(y$Sample.name),]
    type1 <- strsplit(ttest_groups$group[i], " vs. ")[[1]][1]
    type2 <- strsplit(ttest_groups$group[i], " vs. ")[[1]][2]
    y1 <- y$arg_richness[y$sample_type == type1]
    y2 <- y$arg_richness[y$sample_type == type2]
    p_values <- c(p_values, wilcox.test(y1, y2, paired = TRUE, alternative = "two.sided")$p.value)
  }
  ttest_groups$pvalue <- p_values
  asterisk <- rep(NA, length(p_values))
  asterisk[p_values < 0.05] <- "*"
  asterisk[p_values < 0.01] <- "**"
  asterisk[p_values < 0.001] <- "***"
  ttest_groups$asterisk <- asterisk

  # Order ttest_groups
  ttest_groups <- ttest_groups[order(ttest_groups$Location, ttest_groups$group),]
  return(ttest_groups)
}

#' Plots multiple box and whisker diagrams of the ARG richness comparison between paired sample types for every Location
#'
#' @param ttest_groups A dataframe of p-values, asterisks, countries and groups.
#' @param df_map_subsampled_pairs A dataframe of paired samples containing group column.
#'
#' @return A list of box and whisker ggplot objects.
#'
#' @examples
#' df_map_subsampled <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/subsampled_argrich_merged.csv", without_US_duplicates = TRUE)
#' pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"), c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"), c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
#' df_map_subsampled_pairs <- createPairedData(df_map_subsampled, pair_list)
#' ttest_groups <- runTtest(df_map_subsampled_pairs)
#'
#' g <- plotMultipleARGRichnessGraphs(ttest_groups, df_map_subsampled_pairs)
#'
#' @export
#'
#' @import ggplot2
plotMultipleARGRichnessGraphs <- function(ttest_groups, df_map_subsampled_pairs){

  df_map_argrich <- df_map_subsampled_pairs[!duplicated(paste0(df_map_subsampled_pairs$ID, df_map_subsampled_pairs$group)),]
  g <- list()
  set.seed(42)
  for(i in 1:nrow(ttest_groups)){
    df_map_argrich_group <- df_map_argrich[df_map_argrich$Location == ttest_groups$Location[i] & df_map_argrich$group == ttest_groups$group[i],]
    g[[i]] <- ggplot(df_map_argrich_group, aes(sample_type, arg_richness)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 0.8) +
      theme_classic() +
      ylab("ARG Richness") +
      xlab("") +
      ylim(c(0, max(df_map_argrich$arg_richness) + 20)) +
      ggtitle(paste(ttest_groups[i,]$Location, sep = "\n")) +
      geom_text(data = ttest_groups[i,], aes(label=asterisk),
                x = 1.5, y = max(df_map_argrich_group$arg_richness)+10, size = 7,
                inherit.aes = FALSE)
  }
  return(g)
}
