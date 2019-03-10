#' Runs differential analysis using DESeq2 for a Location between paired sample types
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param Location A string of a Location that is found in the Location column of df_map
#' @param compare_samples A character vector of length two of sample types to compare
#'
#' @return A list of length two of a dataframe of DESeq2 results and a character vector of sample types compared in DESeq2
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#'
#' output_list <- runDESeq2(df_map, "US", c("stool", "dorsum of tongue"))
#'
#' @export
#'
#' @importFrom reshape2 dcast
#' @import DESeq2
runDESeq2 <- function(df_map, Location, compare_samples = c("stool", "dorsum of tongue")){

  # Extract Location
  df_map_Location <- df_map[df_map$Location == Location,]

  # Get sample names and samples types that contain paired samples
  id_type_list <- sapply(compare_samples, function(x) df_map_Location$Sample.name[df_map_Location$sample_type == x])
  sample_names <- Reduce(intersect, id_type_list)
  df_map_Location <- df_map_Location[df_map_Location$Sample.name %in% sample_names,]

  df_map_Location <- df_map_Location[df_map_Location$sample_type %in% compare_samples,]

  # Create count matrix
  Location_count <- dcast(data = df_map_Location, formula = V1 ~ ID, fun.aggregate = sum, value.var = "V4")
  arg_names <- Location_count$V1
  row.names(Location_count) <- arg_names
  Location_count <- Location_count[,names(Location_count) != "V1"]
  Location_count <- as.matrix(Location_count)

  # Create Location sample metadata
  Location_meta <- df_map_Location[,c("ID", "Sample.name", "sample_type")]
  Location_meta <- Location_meta[!duplicated(Location_meta$ID),]
  Location_meta$sample_type <- factor(Location_meta$sample_type)
  Location_meta$Sample.name <- factor(Location_meta$Sample.name)

  # Create model matrix
  Location_model <- model.matrix( ~ Sample.name + sample_type, Location_meta)

  # Run DESeq2

  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  dds <- DESeqDataSetFromMatrix(countData = Location_count, colData = Location_meta, design = ~ sample_type)
  dds <- dds[rowSums(counts(dds)) > 5,]
  geoMeans = apply(counts(dds), 1, gm_mean)
  dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  dds = DESeq(dds, fitType="local")

  contrast_list <- resultsNames(dds)
  sample_comparisons <- c(contrast_list[length(contrast_list)-1], contrast_list[length(contrast_list)])
  sample_comp <- strsplit(sample_comparisons[2], "sample_type_")[[1]][2]
  Location_results <- results(dds, contrast=list(sample_comparisons))

  return(list(results = Location_results, sample_comparison = sample_comp))
}

#' Plot Volcano plots of differential analysis
#'
#' @param Location_results A dataframe of DESeq2 results for one Location
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param title A string of the Location name for title of volcano plot
#'
#' @return A ggplot2 object of a volcano plot
#'
#' @examples
#' df_map <- readMappingData("/home/vicky/Documents/CHMI/Resistome-paper/resistomeAnalysis/db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
#' output_list <- runDESeq2(df_map, "US", c("stool", "dorsum of tongue"))
#' g <- plotVolcano(output_list$results, df_map, "stool vs. saliva", "China")
#'
#' @export
#'
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @import RColorBrewer
plotVolcano <- function(Location_results, df_map, title){

  # Extract as data frame
  topT <- as.data.frame(Location_results)

  # Extract ARG labels
  topT$ARG <- sapply(row.names(Location_results), function(x) strsplit(x, "\\|")[[1]][length(strsplit(x, "\\|")[[1]])])

  # Rename some ARG labels
  topT$ARG[topT$ARG == "Bifidobacterium_ileS_conferring_resistance_to_mupirocin"] <- "ileS"

  # Remove p adjust with NA
  topT <- topT[!is.na(topT$padj),]

  # Get class colours
  topT$drug_class <- sapply(row.names(topT), function(x) unique(df_map$Drug.Class_mod[df_map$V1 == x]))

  # Volcano plot
  volcano_plot <- ggplot() +
  geom_point(data = topT[topT$padj >= 0.05 & abs(topT$log2FoldChange) <= 2,],
             aes(log2FoldChange, -log10(padj)), colour = "black", shape = 1, fill = NA) +
  geom_point(data = topT[topT$padj < 0.05 & abs(topT$log2FoldChange) > 2,],
             aes(log2FoldChange, -log10(padj), color = "red"),
             size = 2, shape = 1, fill = NA) +
  geom_label_repel(data = topT[topT$padj < 0.05 & abs(topT$log2FoldChange) > 2,],
            aes(x=log2FoldChange, y=-log10(padj), label=ARG, fill = drug_class)) +
  xlim(c(-25,25)) +
  ylab("-Log10(adjusted p-value)\n") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=14),
        title = element_text(size=20)) +
  xlab("Log2 fold change") +
  ggtitle(title) +
  scale_fill_manual("ARG Class", breaks = names(class_colours), values = class_colours) +
  guides(colour=FALSE) +
  geom_hline(aes(yintercept = -log10(max(topT$padj[topT$padj<0.05], na.rm=TRUE))), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = -2), linetype = "dashed") +
  geom_vline(aes(xintercept = 2), linetype = "dashed")
  return(volcano_plot)
}
