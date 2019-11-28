#' Draws a dendrogram of hierachical clustering of ARG rpkm for longitudinal US samples
#'
#' @param df_map_dup A dataframe of combined non-subsampled or subsampled mapping data and metadata, including longitudinal samples
#' @param coloured_labels A list of colours named by sample type
#'
#' @return None
#'
#' @examples
#' library(RColorBrewer)
#' df_map_dup <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = FALSE)
#' top_col <- rev(brewer.pal(8, "Spectral"))
#' coloured_labels <- c("stool" = top_col[1], "dorsum of tongue" = top_col[2], "buccal mucosa" = top_col[3], "dental" = top_col[4])
#' createUSDendrogram(df_map_dup, coloured_labels)
#'
#' @export
#'
#' @importFrom ape as.phylo
#' @importFrom reshape2 dcast
createUSDendrogram <-  function(df_map_dup, coloured_labels){
  # Select hmp
  df_map_hmp <- df_map_dup[df_map_dup$Location == "US",]
  body_sites <- unique(df_map_hmp$sample_type)
  for(i in 1:length(body_sites)){
    df_map_hmp_body_site <- df_map_hmp[df_map_hmp$sample_type == body_sites[i],]
    df_samples_one <- unique(df_map_hmp_body_site$Sample.name[df_map_hmp_body_site$Visit_Number == 1])
    df_samples_two <- unique(df_map_hmp_body_site$Sample.name[df_map_hmp_body_site$Visit_Number == 2])
    df_samples_rm <- df_samples_one[!(df_samples_one %in% df_samples_two)]
    df_ids_rm <- df_map_hmp_body_site$ID[df_map_hmp_body_site$Sample.name %in% df_samples_rm &
                                           df_map_hmp_body_site$sample_type == body_sites[i]]
    df_map_hmp <- df_map_hmp[!(df_map_hmp$ID %in% df_ids_rm),]
  }

  # Matrix for heatmap
  hmp_rpkm <- dcast(data = df_map_hmp, formula = V1 ~ ID + Sample.name + sample_type, fun.aggregate = sum, value.var = "rpkm")
  row.names(hmp_rpkm) <- hmp_rpkm$V1
  hmp_rpkm <- hmp_rpkm[,-1]
  hmp_rpkm_log <- log10(hmp_rpkm+1)

  # Remove non-duplicates
  colnames_mod <- as.character(sapply(names(hmp_rpkm_log), function(x) paste0(strsplit(x, "_")[[1]][c(2,3)], collapse = "_")))
  hmp_rpkm_log_filt <- hmp_rpkm_log[,colnames_mod %in% unique(colnames_mod[duplicated(colnames_mod)])]
  sample_type_labels <- sapply(names(hmp_rpkm_log_filt), function(x) strsplit(x, "_")[[1]][3])
  names(hmp_rpkm_log_filt) <- sapply(names(hmp_rpkm_log_filt), function(x) strsplit(x, "_")[[1]][2])

  # Hierarchical cluster
  set.seed(1)
  hc <- hclust(dist(t(hmp_rpkm_log_filt)))
  sample_type_colours <- rep(NA, length(sample_type_labels))
  for(i in 1:length(coloured_labels)){
    sample_type_colours[grepl(names(coloured_labels[i]), sample_type_labels)] <- coloured_labels[i]
  }

  # Draw dendrogram
  par(mar = c(0, 1, 0, 0))
  plot(as.phylo(hc), type = "fan", tip.color = sample_type_colours, cex = 0.5)
  legend("bottomright", legend = names(coloured_labels), fill = coloured_labels)

}
