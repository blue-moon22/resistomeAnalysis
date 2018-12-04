#' Combine the Metaphlan relative abundances for all samples
#'
#' @param filenames A character vector of filenames of relative abundance output data from Metaphlan
#' @param output_file A string of the file path of the output csv file
#'
#' @return None
#'
#' @export
#'
#' @importFrom dplyr full_join
combineMetaphlanSamples <- function(filenames, output_file){
  metaphlan <- NULL
  for(i in 1:length(filenames)){
    df_map <- read.delim(filenames[i], stringsAsFactors = FALSE)
    sample <- gsub("_profile.txt", "", strsplit(filenames[i], "/")[[1]][4])
    names(df_map) <- c("Taxa", sample)
    if(is.null(metaphlan)){
      metaphlan <- df_map
    } else {
      metaphlan <- full_join(metaphlan, df_map, by = "Taxa")
    }
  }
  metaphlan[is.na(metaphlan)] <- 0

  write.csv(metaphlan, output_file, row.names = FALSE)

}
