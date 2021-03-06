% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/resistanceMapUse.R
\name{plotGraph}
\alias{plotGraph}
\title{Plot scatter graph of proportion against DDD per 1000}
\usage{
plotGraph(data_merge, Location, class_cols, shape_values)
}
\arguments{
\item{data_merge}{A dataframe of proportions containing ARG classes and confidence intervals, and resistanceMap DDD per 1000}

\item{Location}{A string of a Location e.g. "US"}

\item{class_cols}{A character vector of colours, named by ARG class}
}
\value{
None
}
\description{
Plot scatter graph of proportion against DDD per 1000
}
\examples{
library(RColorBrewer)
df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
df_map_pb_saliva_class <- joinProportionAndBootstrap(df_map, "Drug.Class", "saliva")
df_map_pb_dental_class <- joinProportionAndBootstrap(df_map, "Drug.Class", "dental")
rm_use_card <- readResistanceMapUseData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/ResistanceMap/resistanceMap_use.csv")
data_merge <- mergeData(list(saliva = df_map_pb_saliva_class, dental = df_map_pb_dental_class), rm_use_card)

cols <- brewer.pal(length(unique(data_merge$class)), "Paired")
names(cols) <- unique(data_merge$class)
plotGraph(data_merge, "US", cols)

}
