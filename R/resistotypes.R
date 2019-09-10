#' Run principal coordinate analysis on mapping data and metadata
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#'
#' @return A dataframe of multidimensional scaled coordinates
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' mds <- runPrincipleCoordinateAnalysis(df_map)
#'
#' @export
#'
#' @import dplyr
#' @importFrom reshape2 dcast
#' @importFrom vegan wcmdscale
runPrincipleCoordinateAnalysis <- function(df_map){

  # Create rpkm matrix for PCoA
  rpkm <- dcast(data = df_map, formula = V1 ~ ID, fun.aggregate = sum, value.var = "rpkm")
  row.names(rpkm) <- rpkm$V1
  rpkm <- rpkm[,-1]
  rpkm <- t(rpkm)

  # PCoA
  dst <- dist(rpkm, method="binary")
  mds <- wcmdscale(dst, w=rep(1,nrow(rpkm)))
  mds <- as.data.frame(mds)
  mds$sample_type <- df_map[!(duplicated(df_map$ID)),]$sample_type
  mds$Location <- df_map[!(duplicated(df_map$ID)),]$Location
  mds$group <- paste(mds$sample_type, mds$Location, sep = " - ")

  set.seed(1)
  # Hierarchical clustering
  dst=dist(mds[,!(names(mds) %in% c("sample_type", "Location", "group"))])
  dst_hclust <- hclust(dst)
  # Silhoette analysis
  avg_sil <- numeric(4)
  for(k in 2:(length(avg_sil)+1)) {
    tmp <- silhouette(cutree(dst_hclust, k = k), dst)
    avg_sil[k-1] <- mean(tmp[,3])
  }

  # Show clusters on MDS
  mds$clusters <- paste0("R", as.character(cutree(dst_hclust, which.max(avg_sil)+1)))

  return(mds)
}

#' Plot principal coordinates with resistotype labels
#'
#' @param mds A dataframe of multidimensional scaled coordinates
#'
#' @return None
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' mds <- runPrincipleCoordinateAnalysis(df_map)
#' plotResistotypes(mds)
#'
#' @export
#'
#' @import ggplot2
#' @importFrom cluster silhouette
#' @importFrom RColorBrewer brewer.pal
plotResistotypes <- function(mds){

  cols <- brewer.pal(length(unique(mds$clusters)), "Dark2")
  ggplot(mds, aes(V1, V2, color = clusters)) +
    geom_point(alpha = 0.6) +
    xlab("PCo 1") + ylab("PCo 2") +
    scale_color_manual(name = "Resistotypes", values = cols) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
}
