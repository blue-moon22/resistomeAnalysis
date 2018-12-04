# Load packages
library(resistomeAnalysis)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(purrr)
library(grid)

#### Preprocessing - Merge bedtools and metaphlan2 outputs ####
metadata_file <- "db/SAMPLES/metadata/metadata_healthy.csv"
card_dir <- "db/CARD_DB/card-data/"
# Merge non-subsampled mapping data
mapping_data_filenames <- list.files("db/MAPPING_DATA/NONSUBSAMPLED", full.names = TRUE)
output_file <- "db/MAPPING_DATA/nonsubsampled_merged.csv"
combineBedtools(mapping_data_filenames, output_file, metadata_file, card_dir, subsampled = FALSE)

# Merge subsampled saliva mapping data for percentage of ARG classes analysis
mapping_data_filenames <- list.files("db/MAPPING_DATA/SUBSAMPLE_SALIVA", full.names = TRUE)
output_file <- "db/MAPPING_DATA/subsampled_saliva_merged.csv"
combineBedtools(mapping_data_filenames, output_file, metadata_file, card_dir, subsampled = TRUE)

# Merge subsampled dental mapping data for percentage of ARG classes analysis
mapping_data_filenames <- list.files("db/MAPPING_DATA/SUBSAMPLE_DENTAL", full.names = TRUE)
output_file <- "db/MAPPING_DATA/subsampled_dental_merged.csv"
combineBedtools(mapping_data_filenames, output_file, metadata_file, card_dir, subsampled = TRUE)

# Merge subsampled mapping data for ARG richness analysis
mapping_data_filenames <- list.files("db/MAPPING_DATA/SUBSAMPLE_ARGRICH", full.names = TRUE)
output_file <- "db/MAPPING_DATA/subsampled_argrich_merged.csv"
combineBedtools(mapping_data_filenames, output_file, metadata_file, card_dir, subsampled = TRUE)

# Merge metaphlan2 output
taxa_filenames <- list.files("db/METAPHLAN/output_files", full.names = TRUE)
combineMetaphlanSamples(taxa_filenames, output_file = "db/METAPHLAN/metaphlan.csv")

#### Read mapping data ####
# Read non-subsampled mapping data
df_map <- readMappingData("db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = TRUE)
df_map_dup <- readMappingData("db/MAPPING_DATA/nonsubsampled_merged.csv", without_US_duplicates = FALSE)
# Read subsampled mapping data for ARG richness
df_map_subsampled <- readMappingData("db/MAPPING_DATA/subsampled_argrich_merged.csv", without_US_duplicates = TRUE)
# Read subsampled saliva mapping data
df_map_sub_saliva <- readMappingData("db/MAPPING_DATA/subsampled_saliva_merged.csv", without_US_duplicates = TRUE)
# Read subsampled dental mapping data
df_map_sub_dental <- readMappingData("db/MAPPING_DATA/subsampled_dental_merged.csv", without_US_duplicates = TRUE)

#### Figure 1a ####
# Select colours for graph
cols <- c("grey20", "grey50", "grey80", "white")
names(cols) <- c("China", "Philippines", "UK", "US")

# Percentage saliva samples containing ARG class
df_map_pb_saliva_class <- joinProportionAndBootstrap(df_map, "Drug.Class", "saliva")

# Generate figure
tiff(filename = "figures/Figure1a.tiff", width = 2100, height = 1000, res = 180)
plotPercentages(df_map_pb_saliva_class, cols) + xlab("ARG class") + ylab("% saliva samples") +
  theme(axis.title=element_text(size=18), legend.text = element_text(size=14), legend.title = element_text(size=18))
dev.off()

#### Figure 1b ####
# Percentage saliva samples containing ARG mechanism
df_map_pb <- joinProportionAndBootstrap(df_map, "Resistance.Mechanism", "saliva")

# Generate figure
tiff(filename = "figures/Figure1b.tiff", width = 800, height = 1000, res = 180)
plotPercentages(df_map_pb, cols) + xlab("ARG mechanism") + ylab("% saliva samples") +
  theme(axis.title=element_text(size=18), legend.text = element_text(size=14), legend.title = element_text(size=18))
dev.off()

#### Figure 1c ####
# Percentage dental samples contraining ARG class
df_map_pb_dental_class <- joinProportionAndBootstrap(df_map, "Drug.Class", "dental")

# Generate figure
tiff(filename = "figures/Figure1c.tiff", width = 1250, height = 1000, res = 170)
plotPercentages(df_map_pb_dental_class, cols) + xlab("ARG class") + ylab("% dental samples") +
  theme(axis.title=element_text(size=18), legend.text = element_text(size=14), legend.title = element_text(size=18))
dev.off()

#### Figure 1d ####
# Percentage dental samples containing ARG mechanism
df_map_pb <- joinProportionAndBootstrap(df_map, "Resistance.Mechanism", "dental")

# Generate figure
tiff(filename = "figures/Figure1d.tiff", width = 700, height = 1000, res = 180)
plotPercentages(df_map_pb, cols) + xlab("ARG mechanism") + ylab("% dental samples") +
  theme(axis.title=element_text(size=18), legend.text = element_text(size=14), legend.title = element_text(size=18))
dev.off()

#### Figure 1e ####
# Read resistance map use data
rm_use_card <- readResistanceMapUseData("db/ResistanceMap/resistanceMap_use.csv")

# Merge percentages data for saliva and dental classes
data_merge <- mergeData(list(saliva = df_map_pb_saliva_class, dental = df_map_pb_dental_class), rm_use_card)

# Get class colours
cols <- brewer.pal(length(unique(data_merge$class)), "Paired")
names(cols) <- unique(data_merge$class)

# Create graphs
countries <- unique(data_merge$Country)
g <- list()
for(i in 1:length(countries)){
  g[[i]] <- plotGraph(data_merge, countries[i], cols) +
    theme(legend.position = "none", axis.title.y = element_text(size=14), plot.title = element_text(size=18))
}

# Add legend to grob list
g[[length(g)+1]] <- g_legend(plotGraph(data_merge, unique(data_merge$Country), cols)
                             + theme(legend.title = element_text(size=16), legend.text = element_text(size=12)))

# Generate figure
tiff("figures/Figure1e.tiff", width = 1500, height = 750, res = 150)
grid.arrange(grobs = g, layout_matrix = rbind(c(1,2,5),c(3,4,5)))
dev.off()

#### Figure 2a ####
# Run principle coordinate analysis
mds <- runPrincipleCoordinateAnalysis(df_map)

# Plot resisostypes
tiff("figures/Figure2a.tiff", width = 1600, height = 1000, res = 220)
plotResistotypes(mds) + theme(axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=10))
dev.off()

#### Figure 2b ####
# Label by sample types
tiff("figures/Figure2b.tiff", width = 1600, height = 1000, res = 220)
plotSampleTypes(mds) + theme(axis.title = element_text(size=16), legend.title = element_text(size=16), legend.text = element_text(size=10))
dev.off()

#### Figure 3a ####
# Create dataset of paired samples
pair_list <- list(c("stool", "dental"), c("stool", "saliva"), c("dental", "saliva"),
                  c("stool", "dorsum of tongue"), c("stool", "buccal mucosa"),
                  c("dorsum of tongue", "buccal mucosa"), c("dorsum of tongue", "dental"), c("buccal mucosa", "dental"))
df_map_pairs <- createPairedData(df_map, pair_list)

# Plot boxplots which includes comparison with stool
plotMultipleRPKM <- function(df_map_pairs){
  graphs <- list()
  count <- 0
  unique_Country <- unique(df_map_pairs$Country)
  for(i in 1:length(unique_Country)){
    tmp <- df_map_pairs[df_map_pairs$Country == unique_Country[i],]
    unique_groups <- unique(tmp$group)
    for(j in 1:length(unique_groups)){
      count <- count + 1
      df_map_pairs_groups <- tmp[tmp$group == unique_groups[j],]
      graphs[[count]] <- plotRPKM(df_map_pairs_groups)
    }
  }
  return(graphs)
}

stool_comp_graphs <- plotMultipleRPKM(df_map_pairs[df_map_pairs$group_mod %in% c("stool \nvs. dental", "stool \nvs. saliva", "stool \nvs. buccal mucosa", "stool \nvs. dorsum of tongue"),])
lay <- rbind(c(1, 2, 3), c(4, 5, 6))
tiff("figures/Figure3a.tiff", width = 1800, height = 1000, res = 230)
grid.arrange(grobs = stool_comp_graphs, layout_matrix = lay)
dev.off()

#### Figure 3b ####

# Get relative abundance for each country and sample type with paired oral and stool samples
df_map_rel <- getRelativeAbundance(df_map[df_map$Cohort %in% c("Twin", "Chinese", "HMP1"),])

# Create muliple graphs of relative abundance
uniq_country <- unique(df_map_rel$Country)
cols <- brewer.pal(length(unique(df_map_rel$Drug.Class.Alt)), "Set3")
names(cols) <- unique(df_map_rel$Drug.Class.Alt)
g <- list()
for(i in 1:length(uniq_country)){
  g[[i]] <- plotARGClassAbundance(df_map_rel[df_map_rel$Country == uniq_country[i],], cols) + ggtitle(gsub(" - ", "\n", uniq_country[i])) +
    theme(legend.position = "none", plot.title = element_text(size=18), axis.title = element_text(size=16), axis.text.x = element_text(size=14))
}
g[[length(g)+1]] <- g_legend(plotARGClassAbundance(df_map_rel, cols) + theme(legend.title = element_text(size=18), legend.text = element_text(size=12)))

# Plot relative abundance
lay <- rbind(c(1, 2), c(3, 4))
tiff("figures/Figure3b.tiff", width = 2800, height = 2000, res = 250)
grid.arrange(grobs = g, layout_matrix = lay)
dev.off()

#### Figure 4 ####
# Read subsampled data for ARG richness analysis
# Create dataset of paired samples
df_map_subsampled_pairs <- createPairedData(df_map_subsampled, pair_list)

# Run t-test
ttest_groups <- runTtest(df_map_subsampled_pairs)

# Plot boxplots which includes comparison with stool
graphs_stool_comp <- plotMultipleARGRichnessGraphs(ttest_groups[ttest_groups$group %in% c("stool vs. dental", "stool vs. saliva", "stool vs. buccal mucosa", "stool vs. dorsum of tongue"),],
                                        df_map_subsampled_pairs[df_map_subsampled_pairs$group %in% c("stool vs. dental", "stool vs. saliva", "stool vs. buccal mucosa", "stool vs. dorsum of tongue"),])
tiff("figures/Figure4.tiff", width = 1600, height = 1000, res = 200)
grid.arrange(grobs = graphs_stool_comp, ncol = 3)
dev.off()

#### Figure 5a ####
# Combine metaphlan and rpkm for each sample
metaphlan <- read.csv("db/METAPHLAN/metaphlan.csv", stringsAsFactors = FALSE, row.names = 1)
metaphlan_rpkm <- combineMetaphlanandARG(df_map, metaphlan)

# Chinese healthy saliva with paired stool
sample_ids <- Reduce(intersect, list(unique(df_map$Sample.name[df_map$Country == "China" & df_map$sample_type == "saliva"]),
                                     unique(df_map$Sample.name[df_map$Country == "China" & df_map$sample_type == "stool"])))


# Get spearman's correlation
high_cor_china_saliva <- getSpearmanCorrelation(metaphlan_rpkm, ids = unique(df_map$ID[df_map$Sample.name %in% sample_ids & df_map$sample_type == "saliva"]), taxon_level = "t", taxon_ignore = "@@@")

# Plot heatmap
phyla <- c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Candidatus_Saccharibacteria", "Fusobacteria", "Spirochaetes",
            "Verrucomicrobia", "Ascomycota", "Synergistetes")
tiff("figures/Figure5a.tiff", width = 4500, height = 4000, res = 400)
drawCorrelationHeatmap(high_cor_china_saliva, 150, phyla = phyla)
dev.off()

#### Figure 5b ####
# Spearman's correlation of UK saliva
high_cor_uk_oral <- getSpearmanCorrelation(metaphlan_rpkm, ids = unique(df_map$ID[df_map$Cohort == "Twin" & df_map$sample_type == "saliva"]), taxon_level = "t", taxon_ignore = "@@@")

# Plot heatmap
tiff("figures/Figure5b.tiff", width = 2500, height = 2500, res = 300)
drawCorrelationHeatmap(high_cor_uk_oral, 150, 5, phyla = phyla)
dev.off()

#### Figure 6c ####
# Spearman's correlation of Philippines saliva
high_cor_philippines_saliva <- getSpearmanCorrelation(metaphlan_rpkm, ids = unique(df_map$ID[df_map$Cohort == "Philippines" & df_map$sample_type == "saliva"]), taxon_level = "t", taxon_ignore = "@@@")

# Plot heatmap
tiff("figures/Figure5c.tiff", width = 3500, height = 2800, res = 300)
drawCorrelationHeatmap(high_cor_philippines_saliva, 150, phyla = phyla)
dev.off()

#### Supplementary Figure 1a ####
# Percentage saliva samples containing ARG class
# Get proportion and bootstrap from subsampled and non-subsampled
df_map_pb_comb <- combineSubsampled(df_map, df_map_sub_saliva, level = "Drug.Class", sample_type="saliva")

# Plot
tiff("figures/Supplementary_Figure1a.tiff", width = 2000, height = 1000, res = 180)
plotScatter(df_map_pb_comb)
dev.off()

#### Supplementary Figure 1b ####
# Percentage saliva samples containing ARG mechanism
# Get proportion and bootstrap from subsampled and non-subsampled
df_map_pb_comb_m <- combineSubsampled(df_map, df_map_sub_saliva, level = "Resistance.Mechanism", sample_type="saliva")

# Plot
tiff("figures/Supplementary_Figure1b.tiff", width = 1200, height = 1000, res = 180)
plotScatter(df_map_pb_comb_m)
dev.off()

#### Supplementary Figure 1c ####
# Percentage dental samples containing ARG class
# Get proportion and bootstrap from subsampled and non-subsampled
df_map_pb_comb_dental <- combineSubsampled(df_map, df_map_sub_dental, level = "Drug.Class", sample_type="dental")

# Plot
tiff("figures/Supplementary_Figure1c.tiff", width = 2000, height = 1000, res = 180)
plotScatter(df_map_pb_comb_dental)
dev.off()

#### Supplementary Figure 1d ####
# Percentage dental samples containing ARG mechanism
# Get proportion and bootstrap from subsampled and non-subsampled
df_map_pb_comb_dental_m <- combineSubsampled(df_map, df_map_sub_dental, level = "Resistance.Mechanism", sample_type="dental")

# Plot
tiff("figures/Supplementary_Figure1d.tiff", width = 1200, height = 1000, res = 180)
plotScatter(df_map_pb_comb_dental_m)
dev.off()

#### Supplementary Figure 2a ####
# Create RPKM heatmap for China
top_col <- brewer.pal(6, "Paired")
col_vector <- c("stool" = top_col[1], "saliva" = top_col[3], "dental" = top_col[5])

tiff(filename = "figures/Supplementary_Figure2a.tiff", width = 3000, height = 1500, res = 300)
createRPKMHeatmap(df_map, "Chinese", col_vector)
dev.off()

#### Supplementary Figure 2b ####
# Create RPKM heatmap for US (including duplicates)
top_col <- rev(brewer.pal(8, "Spectral"))
col_vector <- c("stool" = top_col[1], "dorsum of tongue" = top_col[2], "buccal mucosa" = top_col[3], "dental" = top_col[4])

tiff(filename = "figures/Supplementary_Figure2b.tiff", width = 3000, height = 1500, res = 300)
createRPKMHeatmap(df_map_dup, "HMP1", col_vector)
dev.off()

#### Supplementary Figure 2c ####
# Create RPKM heatmap for UK
top_col <- brewer.pal(4, "Paired")
col_vector <- c("stool" = top_col[1], "saliva" = top_col[3])

tiff(filename = "figures/Supplementary_Figure2c.tiff", width = 4000, height = 1500, res = 300)
createRPKMHeatmap(df_map, "Twin", col_vector, show_column_names = TRUE)
dev.off()

#### Supplementary Figure 3 ####
# Create dendrogram of US duplicates
top_col <- rev(brewer.pal(8, "Spectral"))
coloured_labels <- c("stool" = top_col[1], "dorsum of tongue" = top_col[2], "buccal mucosa" = top_col[3], "dental" = top_col[4])

tiff(filename = "figures/Supplementary_Figure3.tiff", width = 3000, height = 3000, res = 200)
createUSDendrogram(df_map_dup, coloured_labels)
dev.off()

#### Supplementary Figure 4 ####
# Get individual relative abundance of ARGs by ARG class
df_map_rel_ind <- getRelativeAbundanceIndividuals(df_map[df_map$Cohort %in% c("Twin", "Chinese", "HMP1"),])

# Create multiple plots of relative abundance
uniq_Country <- unique(df_map_rel_ind$Country)
cols <- brewer.pal(length(unique(df_map_rel_ind$Drug.Class.Alt)), "Set3")
names(cols) <- unique(df_map_rel_ind$Drug.Class.Alt)

# Plot muliple graphs
g <- list()
for(i in 1:length(uniq_Country)){
  g[[i]] <- plotIndividualAbundance(df_map_rel_ind[df_map_rel_ind$Country == uniq_Country[i],], cols) + ggtitle(gsub(" - ", "\n", uniq_Country[i])) +
    theme(legend.position = "none", plot.title = element_text(size=18), axis.title = element_text(size=16))
}
g[[length(g)+1]] <- g_legend(plotIndividualAbundance(df_map_rel_ind, cols) + theme(legend.title = element_text(size=18), legend.text = element_text(size=12)))

# Plot relative abundance
lay <- rbind(c(1, 2, 3, 4))
tiff("figures/Supplementary_Figure4.tiff", width = 3500, height = 1000, res = 160)
grid.arrange(grobs = g, layout_matrix = lay)
dev.off()

#### Supplementary Figure 5a->j ####

list_deseq <- list(China = list(sample_type = c("saliva", "stool", "dental")),
                   US = list(sample_type = c("stool", "buccal mucosa", "dorsum of tongue", "dental")),
                   UK = list(sample_type = c("saliva", "stool")))

# DeSeq2 analysis through each comparison
results_list <- list()
name <- list()
sample_comps <- list()
count <- 0
for(i in 1:length(list_deseq)){
  for(k in 1:(length(list_deseq[[i]]$sample_type)-1)){
    for(l in (k+1):length(list_deseq[[i]]$sample_type)){
      count = count + 1
      deseq_output <- runDESeq2(df_map[df_map$Cohort %in% c("Chinese", "HMP1", "Twin"),], Country = names(list_deseq[i]),
                                compare_samples = c(list_deseq[[i]]$sample_type[k], list_deseq[[i]]$sample_type[l]))
      results_list[[count]] <- deseq_output$results
      sample_comps[[count]] <- deseq_output$sample_comparison
      name[[count]] <- paste(names(list_deseq[i]), list_deseq[[i]]$sample_type[k], list_deseq[[i]]$sample_type[l], sep= "_")
    }
  }
}

# Create volcano plots
countries <- sapply(unlist(name), function(x) strsplit(x, "_")[[1]][1])
sample_comps <- gsub("\\.", " ", sample_comps)
sample_comp_labels <- lapply(sample_comps, function(x) strsplit(x, "_vs_")[[1]][c(2,1)])
sample_comps <- gsub("_", " ", sample_comps)
titles <- map2_chr(.x=countries, .y=sample_comps, .f = function(.x, .y) paste0(.x, "\n", .y))
p <- list() # Print plots out individually and save in grob list
for(i in 1:length(results_list)){
  p[[i]] <- plotVolcano(results_list[[i]], df_map[df_map$Cohort %in% c("Chinese", "HMP1", "Twin"),], sample_comp_labels[[i]], titles[i])
}

# Print them to files
for(i in 1:length(p)){
  filename <- paste0("figures/Supplementary_Figure5", letters[i], ".tiff")
  tiff(filename, width = 5000, height = 3000, res = 250)
  print(grid.arrange(arrangeGrob = p[[i]], layout_matrix = matrix(1, 1, 1),
                     bottom=textGrob("Log2 fold change", just = c("centre", "bottom"), vjust=1, gp = gpar(fontsize=20)),
                     vp=viewport(width=0.95, height=0.95)))
  dev.off()
}

#### Supplementary Figure 6 ####
# Plot ARG richness comparisons that do not contain a stool comparison
graphs_oral_comp <- plotMultipleARGRichnessGraphs(ttest_groups[ttest_groups$group %in% c("dental vs. saliva", "buccal mucosa vs. dental",
                                                                              "dorsum of tongue vs. buccal mucosa", "dorsum of tongue vs. dental"),],
                                       df_map_subsampled_pairs[df_map_subsampled_pairs$group %in% c("dental vs. saliva", "buccal mucosa vs. dental",
                                                                          "dorsum of tongue vs. buccal mucosa", "dorsum of tongue vs. dental"),])
tiff("figures/Supplementary_Figure6.tiff", width = 2200, height = 500, res = 200)
grid.arrange(grobs = graphs_oral_comp, nrow = 1)
dev.off()

#### Supplementary Figure 7 ####
# Get spearman's correlation
high_cor_china_stool <- getSpearmanCorrelation(metaphlan_rpkm, ids = unique(df_map$ID[df_map$Sample.name %in% sample_ids & df_map$sample_type == "stool"]), taxon_level = "t", taxon_ignore = "@@@")

# Plot heatmap
tiff("figures/Supplementary_Figure7.tiff", width = 3000, height = 5000, res = 400)
drawCorrelationHeatmap(high_cor_china_stool, 150, phyla = phyla)
dev.off()
