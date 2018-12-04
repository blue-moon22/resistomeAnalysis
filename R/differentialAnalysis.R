#' Runs differential analysis using DESeq2 for a country between paired sample types
#'
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param Country A string of a country that is found in the Country column of df_map
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
runDESeq2 <- function(df_map, Country, compare_samples = c("stool", "dorsum of tongue")){

  # Extract Country
  df_map_Country <- df_map[df_map$Country == Country,]

  # Get sample names and samples types that contain paired samples
  id_type_list <- sapply(compare_samples, function(x) df_map_Country$Sample.name[df_map_Country$sample_type == x])
  sample_names <- Reduce(intersect, id_type_list)
  df_map_Country <- df_map_Country[df_map_Country$Sample.name %in% sample_names,]

  df_map_Country <- df_map_Country[df_map_Country$sample_type %in% compare_samples,]

  # Create count matrix
  Country_count <- dcast(data = df_map_Country, formula = V1 ~ ID, fun.aggregate = sum, value.var = "V4")
  arg_names <- Country_count$V1
  row.names(Country_count) <- arg_names
  Country_count <- Country_count[,names(Country_count) != "V1"]
  Country_count <- as.matrix(Country_count)

  # Create Country sample metadata
  Country_meta <- df_map_Country[,c("ID", "Sample.name", "sample_type")]
  Country_meta <- Country_meta[!duplicated(Country_meta$ID),]
  Country_meta$sample_type <- factor(Country_meta$sample_type)
  Country_meta$Sample.name <- factor(Country_meta$Sample.name)

  # Create model matrix
  Country_model <- model.matrix( ~ Sample.name + sample_type, Country_meta)

  # Run DESeq2

  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  dds <- DESeqDataSetFromMatrix(countData = Country_count, colData = Country_meta, design = ~ sample_type)
  dds <- dds[rowSums(counts(dds)) > 5,]
  geoMeans = apply(counts(dds), 1, gm_mean)
  dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  dds = DESeq(dds, fitType="local")

  contrast_list <- resultsNames(dds)
  sample_comparisons <- c(contrast_list[length(contrast_list)-1], contrast_list[length(contrast_list)])
  sample_comp <- strsplit(sample_comparisons[2], "sample_type_")[[1]][2]
  Country_results <- results(dds, contrast=list(sample_comparisons))

  return(list(results = Country_results, sample_comparison = sample_comp))
}

#' Plot Volcano plots of differential analysis
#'
#' @param Country_results A dataframe of DESeq2 results for one country
#' @param df_map A dataframe of combined non-subsampled or subsampled mapping data and metadata
#' @param sample_comp A string of sample type comparison e.g. "saliva vs. stool"
#' @param title A string of the country name for title of volcano plot
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
plotVolcano <- function(Country_results, df_map, sample_comp, title){

  # Extract as data frame
  topT <- as.data.frame(Country_results)

  # Extract ARG labels
  topT$ARG <- sapply(row.names(Country_results), function(x) strsplit(x, "\\|")[[1]][6])

  # Rename some ARG labels
  topT$ARG[topT$ARG == "Bifidobacterium_ileS_conferring_resistance_to_mupirocin"] <- "ileS"

  # Remove p adjust with NA
  topT <- topT[!is.na(topT$padj),]

  # Get class colours
  topT$drug_class <- sapply(row.names(topT), function(x) unique(df_map$Drug.Class_mod[df_map$V1 == x]))
  class_names <- unique(df_map$Drug.Class_mod)
  n <- length(class_names)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  class_colours = rev(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  names(class_colours) <- class_names

  # Volcano plot
  volcano_plot <- ggplot() +
  geom_point(data = topT[topT$padj >= 0.05 & abs(topT$log2FoldChange) <= 2,],
             aes(log2FoldChange, -log10(padj)), colour = "black", shape = 1, fill = NA) +
  geom_point(data = topT[topT$padj < 0.05 & abs(topT$log2FoldChange) > 2,],
             aes(log2FoldChange, -log10(padj), color = "red"),
             size = 2, shape = 1, fill = NA) +
  geom_label_repel(data = topT[topT$padj < 0.05 & abs(topT$log2FoldChange) > 2,],
            aes(x=log2FoldChange, y=-log10(padj), label=ARG, fill = drug_class)) +
  annotate(geom = "text", x = c(-20, 20), y = -20, label = sample_comp, size = 6) +
  xlim(c(-25,25)) +
  ylab("-Log10(adjusted p-value)\n") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=14),
        title = element_text(size=20)) +
  xlab("") +
  ggtitle(title) +
  scale_fill_manual("ARG Class", breaks = names(class_colours), values = class_colours) +
  guides(colour=FALSE) +
  geom_hline(aes(yintercept = -log10(max(topT$padj[topT$padj<0.05], na.rm=TRUE))), linetype = "dashed") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_vline(aes(xintercept = -2), linetype = "dashed") +
  geom_vline(aes(xintercept = 2), linetype = "dashed")
  return(volcano_plot)
}
