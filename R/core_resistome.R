#' Create a circular graph showing the core ARGs for each ARG Class
#'
#' @param data A dataframe with the percentage of samples that contain each ARG
#' @param by_cohort A boolean stating whether the bars for each cohort should be side-by-side (TRUE) stacked (FALSE)
#' @param hjust_val A numerical vector with the position of the ARG class text for each ARG class
#' @param distance_labels A numeric stating distance from the ARG class text and the x-axis
#'
#' @examples
#'
#' @export
#'
#' @import dplyr viridis
createCircularGraph <- function(data_arg, hjust_val = rep(0, length(unique(data_arg$boot_ci$Drug.Class))), distance_labels = -25, bar_label_size, group_label_size){

  data_main <- data_arg$boot_ci
  data_points <- data_arg$boot_perc

  # Set a number of 'empty bar' to add at the end of each group
  empty_bar=length(unique(data_main$Location)) - 1
  nObsType=nlevels(as.factor(data_main$Location))
  data_main$Drug.Class <- as.factor(data_main$Drug.Class)
  data_points$Drug.Class <- as.factor(data_points$Drug.Class)

  to_add = data.frame( matrix(NA, empty_bar*nlevels(data_main$Drug.Class), ncol(data_main)) )
  colnames(to_add) = colnames(data_main)
  to_add$Drug.Class=rep(levels(data_main$Drug.Class), each=1 )

  data_main=rbind(data_main, to_add)
  data_main=data_main %>% arrange(Drug.Class, level)
  data_main$id <- seq(1:nrow(data_main))

  to_add = data.frame( matrix(NA, empty_bar*nlevels(data_points$Drug.Class), ncol(data_points)) )
  colnames(to_add) <- colnames(data_points)
  to_add$Drug.Class=rep(levels(data_points$Drug.Class), each=1 )

  data_points <- rbind(data_points, to_add)
  data_points <- data_points %>% arrange(Drug.Class, level)
  data_points <- left_join(data_points, data_main, by = c("Drug.Class", "level", "Location_health", "Location")) %>%
    select(c("value", "Drug.Class", "level", "Location_health", "Location", "id"))

  df_labels= data_main %>% group_by(id, level, Location, labels) %>% summarize(tot=sum(CI_ub))
  number_of_bar=nrow(df_labels)
  angle= 90 - 360 * (df_labels$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  df_labels$hjust<-ifelse( angle < -90, 1, 0)
  df_labels$angle<-ifelse(angle < -90, angle+180, angle)

  # prepare a data frame for base lines
  base_data=data_main %>%
    group_by(Drug.Class) %>%
    summarize(start=min(id), end=max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))
  log <- base_data$start == base_data$end
  base_data$start[log] <- base_data$start[log] - 0.5
  base_data$end[log] <- base_data$end[log] + 0.5

  # remove word antibiotic and semi-colon from basedata
  base_data$Drug.Class <- gsub("antibiotic", "", base_data$Drug.Class)
  base_data$Drug.Class <- gsub(";", "\n", base_data$Drug.Class)

  # prepare a data frame for grid (scales)
  grid_data = base_data
  grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start = grid_data$start - 1
  grid_data=grid_data[-1,]
  log <- grid_data$start == grid_data$end
  grid_data$start[log] <- grid_data$start[log] - 0.5
  grid_data$end[log] <- grid_data$end[log] + 0.5

  # Core resistome across healthy individuals across countries
  p <- ggplot() +
    geom_bar(data = data_main, mapping = aes(x=as.factor(id), y=proportion, fill = Location),
             stat = "identity", position = "dodge", alpha=0.5) +
    geom_errorbar(data = data_main, mapping = aes(x = as.factor(id), ymin = CI_lb, ymax = CI_ub, group = Location),
                  position = position_dodge(), colour = "blue") +
    geom_point(data = data_points, mapping = aes(x=as.factor(id), y = value, group = Location), shape = 4, position=position_dodge(0.9)) +
    scale_fill_viridis(discrete=TRUE, breaks = sort(unique(data_main$Location[!is.na(data_main$Location)])))  +

    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +

    # Add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(max(df_labels$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +

    ylim(-150,max(df_labels$tot, na.rm=T)+20) +
    theme_minimal() +
    theme(
      #legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(2,4), "cm")
    ) +
    coord_polar() +

    # Add labels on top of each bar
    geom_text(data=df_labels, aes(x=id, y=tot+10, label=labels, hjust=hjust), color="black", fontface="bold", alpha=0.6, size=bar_label_size, angle=df_labels$angle, inherit.aes=FALSE, parse = TRUE) +

    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6, inherit.aes=FALSE) +
    geom_text(data=base_data, aes(x = title, y = distance_labels, label=Drug.Class), hjust=hjust_val, colour = "black", alpha=0.8, size=group_label_size, fontface="bold", inherit.aes = FALSE)

  return(p)
}

#' Draw a circular graph showing the core ARGs for each ARG Class
#'
#' @param data A dataframe of the mapping data
#' @param hjust_val A numerical vector with the position of the ARG class text for each ARG class
#' @bar_label_size Text size of bar label
#' @group_label_size Text size of group label
#' @B A numerical value of the number of types to sample in bootstrapping
#'
#' @examples
#'
#' @export
#'
#' @import dplyr tidyr
drawCoreARGs <- function(data, hjust_val, bar_label_size, group_label_size, B = 100) {
  data_arg <- joinProportionAndBootstrap(data, "ARO.Name", B)
  df_class <- data %>% select(ARO.Name, Drug.Class)
  df_class <- df_class[!(duplicated(df_class$ARO.Name)),]
  data_arg$boot_ci <- left_join(data_arg$boot_ci, df_class, by = c("level" = "ARO.Name"))
  data_arg$boot_perc <- left_join(data_arg$boot_perc, df_class, by = c("level" = "ARO.Name"))
  data_arg$boot_ci$labels <- data_arg$boot_ci$level
  #data_arg <- data_arg %>% transform(Drug.Class = strsplit(Drug.Class, ";")) %>% unnest(Drug.Class)

  # Get core
  data_arg$boot_ci <- data_arg$boot_ci[data_arg$boot_ci$proportion >= 70,]
  data_arg$boot_ci$CI_lb[is.na(data_arg$boot_ci$CI_lb)] <- data_arg$boot_ci$proportion[is.na(data_arg$boot_ci$CI_lb)]
  data_arg$boot_ci$CI_ub[is.na(data_arg$boot_ci$CI_lb)] <- data_arg$boot_ci$proportion[is.na(data_arg$boot_ci$CI_lb)]
  data_arg$boot_ci$labels <- paste0('italic("', gsub(" \\[.*", "", data_arg$boot_ci$labels), '")')

  data_arg$boot_perc <- data_arg$boot_perc[paste(data_arg$boot_perc$level, data_arg$boot_perc$Location) %in% paste(data_arg$boot_ci$level, data_arg$boot_ci$Location),]

  p <- createCircularGraph(data_arg, hjust_val = hjust_val, bar_label_size = bar_label_size, group_label_size = group_label_size)
  return(p)
}
