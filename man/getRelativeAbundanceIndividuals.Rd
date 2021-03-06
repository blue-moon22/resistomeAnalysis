% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stackedBarChartsRPKM.R
\name{getRelativeAbundanceIndividuals}
\alias{getRelativeAbundanceIndividuals}
\title{Get relative abundance of total rpkm per individual, sample type and drug class}
\usage{
getRelativeAbundanceIndividuals(df_map)
}
\arguments{
\item{df_map}{A dataframe of combined non-subsampled or subsampled mapping data and metadata}
}
\value{
A dataframe of relative abundances of drug classes for every individual and sample type
}
\description{
Get relative abundance of total rpkm per individual, sample type and drug class
}
\examples{
df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
df_map_rel_ind <- getRelativeAbundanceIndividuals(df_map)

}
