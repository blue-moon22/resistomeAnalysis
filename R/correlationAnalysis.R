#' Get a dataframe of correlated ARGs and taxa using Spearman's Correlation
#'
#' @param data A dataframe of combined non-subsampled mapping data and metadata
#' @param ids A character vector of sample IDs
#' @param taxon_level A string representing the taxon level to apply Spearman's Correlation, where k is Kingdom, p is Phylum,
#' c is Class, o is Order, f is Family, g is Genus, s is Species and t is Strain
#' @param taxon_ignore A string representing the taxon level to not include in Spearman's Correlation, where k is Kingdom, p is Phylum,
#' c is Class, o is Order, f is Family, g is Genus, s is Species, t is Strain and "" is don't ignore any taxa.
#' It has to be a taxon after taxon_level.
#'
#' @return A dataframe containing correlated ARGs and taxa, rho and adjusted p-values
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' metaphlan <- read.csv("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/METAPHLAN/metaphlan.csv", stringsAsFactors = FALSE, row.names = 1)
#' metaphlan_rpkm <- combineMetaphlanandARG(df_map, metaphlan)
#'
#' high_cor_uk_oral <- getSpearmanCorrelation(metaphlan_rpkm, ids = unique(df_map$ID[df_map$Cohort == "Twin" & df_map$sample_type == "saliva"]), taxon_level = "t", taxon_ignore = "@@@")
#'
#' @export
#'
#' @importFrom magrittr %>%
getSpearmanCorrelation <- function(data, ids = NA, taxon_level, taxon_ignore){

  if(!any(is.na(ids))){
    data <- data[,names(data) %in% ids]
  }

  arg_data <- data[!grepl("__", row.names(data)),]
  arg_data <- arg_data[rowSums(arg_data != 0) > ncol(arg_data)/2,]
  species_data <- data[grepl(paste0(taxon_level, "__"), row.names(data)) & !grepl(paste0(taxon_ignore, "__"), row.names(data)),]
  species_data <- species_data[rowSums(species_data != 0) > ncol(species_data)/2,]
  cor_data <- c()
  cor_test <- c()
  x <- c()
  y <- c()
  # Spearman Correlation p values
  for(i in 1:nrow(species_data)){
    for(j in 1:nrow(arg_data)){
      cor <- cor.test(as.numeric(species_data[i,]), as.numeric(arg_data[j,]), method = "spearman")
      cor_data <- c(cor_data, cor$estimate)
      cor_test <- c(cor_test, cor$p.value)
      x <- c(x, row.names(species_data)[i])
      y <- c(y, row.names(arg_data)[j])
    }
  }

  # Remove repeat tests and direct comparisons
  cor_df <- data.frame(x = x, y = y, rho = cor_data, p_val = cor_test)
  cor_df <- cor_df[paste0(x, y) != paste0(y, x),]

  # Benjamini-Hochberg correction
  cor_df$FDR <- p.adjust(cor_df$p_val, method = "BH")
  cor_df <- cor_df[cor_df$FDR < 0.05,]

  # Re-label
  cor_df$phylum <- sapply(as.character(cor_df$x), function(x) strsplit(x, "p__")[[1]][2]) %>% sapply(function(x) strsplit(x, "\\|")[[1]][1])
  cor_df$x <- gsub("Streptococcus_mitis_oralis_pneumoniae", "Streptococcus_mitis/oralis/pneumoniae", cor_df$x)
  cor_df_tmp_x <- cor_df$x
  cor_df_tmp_x[grep("unclassified", cor_df$x)] <- sapply(cor_df$x[grep("unclassified", cor_df$x)], function(x) strsplit(x, "\\|")[[1]][length(strsplit(x, "\\|")[[1]])-1]) %>%
    sapply(function(x) strsplit(x, "__")[[1]][2])
  cor_df$x[grep("__", cor_df$x)] <- gsub("\\|", "_", sapply(gsub(paste0(taxon_level, "__"), "", cor_df_tmp_x[grep("__", cor_df$x)]), function(x) strsplit(x, "__")[[1]][length(strsplit(x, "__")[[1]])]))

  return(cor_df)
}

#' Draws a heatmap from a dataframe of correlated ARGs and taxa
#'
#' @param cor_df A dataframe containing correlated ARGs and taxa, rho and adjusted p-values
#' @param bottom_margin A numeric value in millimeters for height of column names
#' @param padding A numeric value in millimeters for width for row names
#' @param phyla A character vector of phyla names for labelling columns of heatmap
#'
#' @return None
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' metaphlan <- read.csv("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/METAPHLAN/metaphlan.csv", stringsAsFactors = FALSE, row.names = 1)
#' metaphlan_rpkm <- combineMetaphlanandARG(df_map, metaphlan)
#' high_cor_uk_oral <- getSpearmanCorrelation(metaphlan_rpkm, ids = unique(df_map$ID[df_map$Cohort == "Twin" & df_map$sample_type == "saliva"]), taxon_level = "t", taxon_ignore = "@@@")
#' phyla <- c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Candidatus_Saccharibacteria", "Fusobacteria", "Spirochaetes", "Verrucomicrobia", "Ascomycota", "Synergistetes")
#' h <- drawCorrelationHeatmap(high_cor_uk_oral, 100, 5, phyla)
#'
#' @export
#'
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap draw
#' @importFrom reshape2 dcast
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
drawCorrelationHeatmap <- function(cor_df, bottom_margin = 0, left_margin = 0, phyla){
  cor_df <- cor_df[cor_df$rho > 0,]
  sparse_matrix <- dcast(data = cor_df, formula = y ~ x + phylum, fun.aggregate = sum, value.var = "rho")
  row.names(sparse_matrix) <- sparse_matrix$y
  sparse_matrix <- sparse_matrix[,-1]

  row_cols <- c(brewer.pal(length(phyla), "Spectral"))
  names(row_cols) <- phyla

  phyla_samples <- rep(NA, ncol(sparse_matrix))
  phyla_cols <- rep(NA, ncol(sparse_matrix))
  for(i in 1:length(phyla)){
    phyla_samples[grep(phyla[i], names(sparse_matrix))] <- phyla[i]
    phyla_cols[grep(phyla[i], names(sparse_matrix))] <- row_cols[names(row_cols) == phyla[i]]
    names(sparse_matrix) <- gsub(paste0("_", phyla[i]), "", names(sparse_matrix))
  }
  names(phyla_cols) <- phyla_samples
  row.names(sparse_matrix) <- gsub(" \\[.*", "", row.names(sparse_matrix))
  colnames(sparse_matrix) <- gsub(" G.*", "", gsub("_", " ", colnames(sparse_matrix)))

  ha = HeatmapAnnotation(type = phyla_samples, col = list(type = phyla_cols),
                         annotation_legend_param = list(type = list(title = "Phylum",
                                                                    title_gp = gpar(fontsize = 20),
                                                                    labels_gp = gpar(fontsize = 16, fontface = "italic"),
                                                                    grid_height = unit(8, "mm"))),
                         show_annotation_name = FALSE)

  set.seed(42)
  ht <- Heatmap(sparse_matrix, top_annotation = ha, name = "rho", show_row_names = TRUE, cluster_rows = TRUE,
          col = colorRamp2(c(0,1),  c("white", "red")), column_names_max_height = unit(bottom_margin, "mm"),
          row_title_rot = 0, row_title_gp = gpar(fontsize = 5), show_column_names = TRUE, row_names_max_width = unit(left_margin, "mm"),
          heatmap_legend_param = list(color_bar = "continuous", title_gp = gpar(fontsize = 20),
                                      labels_gp = gpar(fontsize=16), legend_height = unit(8, "cm")),
          row_names_gp = gpar(fontface = "italic"), column_names_gp = gpar(fontface = "italic"))
  draw(ht)
}


