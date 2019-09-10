#' Combine a dataframe containing mapping data and metadata data with a dataframe of Metaphlan relative abundances
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param metaphlan A dataframe of Metaphlan relative abundances of all samples
#'
#' @return A dataframe combining rpkm from mapping data and relative abundance from metaphlan with row names as IDs
#'
#' @export
#'
#' @importFrom reshape2 dcast
combineMetaphlanandARG <- function(df_map, metaphlan){

  # Create rpkm matrix for PCoA
  rpkm <- dcast(data = df_map, formula = ARO.Name ~ ID, fun.aggregate = sum, value.var = "rpkm")
  row.names(rpkm) <- rpkm$ARO.Name
  rpkm <- rpkm[,-1]

  # Samples of metaphlan and args
  metaphlan_samples <- colnames(metaphlan)
  arg_samples <- colnames(rpkm)
  intersect_samples <- Reduce(intersect, list(metaphlan_samples, arg_samples))
  rpkm <- rpkm[,names(rpkm) %in% intersect_samples]
  metaphlan <- metaphlan[,names(metaphlan) %in% intersect_samples]
  df_map_comb <- rbind(rpkm, metaphlan)
  df_map_comb <- df_map_comb[rowSums(df_map_comb) != 0,]

  return(df_map_comb)
}

