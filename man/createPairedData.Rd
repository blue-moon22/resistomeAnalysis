% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boxplotsRPKM.R
\name{createPairedData}
\alias{createPairedData}
\title{Creates a dataframe of paired sample types with a column named group}
\usage{
createPairedData(df_map, pair_list)
}
\arguments{
\item{df_map}{A dataframe of combined non-subsampled or subsampled mapping data and metadata}

\item{pair_list}{A list of vectors of length two containing paired sample types}
}
\value{
A dataframe of paired sample types with a column named group i.e. saliva vs. stool
}
\description{
Creates a dataframe of paired sample types with a column named group
}
\examples{
df_map_subsampled <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/subsampled_argrich_merged.csv", without_US_duplicates = TRUE)
pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"), c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"), c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
df_map_subsampled_pairs <- createPairedData(df_map_subsampled, pair_list)

}
