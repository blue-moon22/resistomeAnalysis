#' Draws a heatmap from a dataframe of ARG abundance for every sample
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param Location A string of the study Location e.g. "China"
#' @param col_vector A character vector of colours, named by sample type
#' @param show_column_names Logical: If TRUE, shows column names on heatmap
#'
#' @return None
#'
#' @examples
#' library(RColorBrewer)
#'
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' top_col <- brewer.pal(4, "Paired")
#' col_vector <- c("stool" = top_col[1], "saliva" = top_col[3])
#' createRPKMHeatmap(df_map, "China", col_vector, show_column_names = TRUE)
#'
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom reshape2 dcast
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
createRPKMHeatmap <- function(df_map, Location, col_vector, show_column_names = FALSE) {

  # Remove antibiotic word from Drug Class
  df_map$Drug.Class <- gsub(" antibiotic", "", df_map$Drug.Class)
  df_map$Drug.Class <- gsub("antibiotic", "", df_map$Drug.Class)

  # Select Location
  df_map_Location <- df_map[df_map$Location == Location,]

  # Matrix for heatmap
  Location_rpkm <- dcast(data = df_map_Location, formula = V1 ~ Sample.name + sample_type, fun.aggregate = sum, value.var = "rpkm")
  row.names(Location_rpkm) <- Location_rpkm$V1
  Location_rpkm <- Location_rpkm[,-1]
  Location_rpkm_log <- log10(Location_rpkm+1)

  # Create labels for top group annotations
  labels <- colnames(Location_rpkm)
  top_labels <- sapply(labels, function(x) strsplit(x, "_")[[1]][length(strsplit(x, "_")[[1]])])
  ha = HeatmapAnnotation(type = top_labels,
                         col = list(type = col_vector),
                         annotation_legend_param = list(type = list(title = "Sample Type")))

  # Change Location_rpkm colnames
  colnames_sample <- sapply(labels, function(x) paste0(strsplit(x, "_")[[1]][-length(strsplit(x, "_")[[1]])], collapse = "_"))
  colnames(Location_rpkm_log) <- as.character(colnames_sample)

  # Create ARG Class labels for row annotations
  args <- row.names(Location_rpkm)
  drug_class <- sapply(args, function(x) df_map_Location$Drug.Class[df_map_Location$V1 == x[1]][1])

  # Create heatmap
  set.seed(1)
  Heatmap(Location_rpkm_log, na_col = "#f2f2f2", top_annotation = ha, name = "log10(RPKM+1)", show_row_names = FALSE, cluster_rows = FALSE,
          col = colorRamp2(c(0, min(Location_rpkm_log[Location_rpkm_log != 0]), max(Location_rpkm_log, na.rm = TRUE)),  c("#f2f2f2", "#f7db04", "red4")),
          split = drug_class, row_title_rot = 0, row_title_gp = gpar(fontsize = 5), show_column_names = show_column_names,
          heatmap_legend_param = list(color_bar = "continuous"))
}
