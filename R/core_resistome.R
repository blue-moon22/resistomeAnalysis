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
createCircularGraph <- function(data, by_cohort = FALSE, hjust_val = rep(0, length(unique(data$Drug.Class))), distance_labels = -25){
  # Scale proportion
  if(!by_cohort){
    data$proportion <- data$proportion / length(unique(data$Location))
  }

  # Set a number of 'empty bar' to add at the end of each group
  empty_bar=length(unique(data$Location)) - 1
  nObsType=nlevels(as.factor(data$Location))
  data$Drug.Class <- as.factor(data$Drug.Class)
  if(by_cohort){
    to_add = data.frame( matrix(NA, empty_bar*nlevels(data$Drug.Class), ncol(data)) )
    colnames(to_add) = colnames(data)
    to_add$Drug.Class=rep(levels(data$Drug.Class), each=1 )
  } else {
    to_add = data.frame( matrix(NA, empty_bar*nlevels(data$Drug.Class)*nObsType, ncol(data)) )
    colnames(to_add) = colnames(data)
    to_add$Drug.Class=rep(levels(data$Drug.Class), each=empty_bar*nObsType )
  }
  data=rbind(data, to_add)
  data=data %>% arrange(Drug.Class, level)
  if(by_cohort){
    data$id <- seq(1:nrow(data))
  } else {
    data$id=rep( seq(1, nrow(data)/nObsType) , each=nObsType)
  }

  # Get the name and the y position of each label
  if(by_cohort){
    df_labels= data %>% group_by(id, level, Location, labels) %>% summarize(tot=sum(CI_ub95))
  } else {
    df_labels= data %>% group_by(id, level, labels) %>% summarize(tot=sum(proportion))
  }
  number_of_bar=nrow(df_labels)
  angle= 90 - 360 * (df_labels$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  df_labels$hjust<-ifelse( angle < -90, 1, 0)
  df_labels$angle<-ifelse(angle < -90, angle+180, angle)

  # prepare a data frame for base lines
  base_data=data %>%
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
  if(by_cohort){
    p <- ggplot(data, aes(x=as.factor(id), y=proportion, fill = Location, ymin = CI_lb95, ymax = CI_ub95)) +
      geom_bar(stat = "identity", position = "dodge", alpha=0.5)
  } else {
    p <- ggplot(data, aes(x=as.factor(id), y=proportion, fill = Location)) +
      geom_bar(stat = "identity", alpha=0.5)
  }
  p <- p +
    geom_errorbar(position = position_dodge()) +
    scale_fill_viridis(discrete=TRUE, breaks = sort(unique(data$Location[!is.na(data$Location)])))  +

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
      plot.margin = unit(rep(-1,4), "cm")
    ) +
    coord_polar() +

    # Add labels on top of each bar
    geom_text(data=df_labels, aes(x=id, y=tot+10, label=labels, hjust=hjust), color="black", fontface="bold", alpha=0.6, size=3, angle=df_labels$angle, inherit.aes=FALSE) +

    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6, inherit.aes=FALSE) +
    geom_text(data=base_data, aes(x = title, y = distance_labels, label=Drug.Class), hjust=hjust_val, colour = "black", alpha=0.8, size=3, fontface="bold", inherit.aes = FALSE)

  return(p)
}

#' Draw a circular graph showing the core ARGs for each ARG Class
#'
#' @param data A dataframe of the mapping data
#' @param sample_type A string stating the sample type e.g. "saliva"
#' @param hjust_val A numerical vector with the position of the ARG class text for each ARG class
#'
#' @examples
#'
#' @export
#'
#' @import dplyr tidyr
drawCoreARGs <- function(data, sample_type, hjust_val) {
  data_arg <- joinProportionAndBootstrap(data, "V1", sample_type)
  df_class <- data %>% select(V1, Drug.Class)
  df_class <- df_class[!(duplicated(df_class$V1)),]
  data_arg <- left_join(data_arg, df_class, by = c("level" = "V1"))
  data_arg <- data_arg %>% transform(Drug.Class = strsplit(Drug.Class, ";")) %>% unnest(Drug.Class)

  # Change label names
  data_arg$labels <- sapply(data_arg$level, function(x) strsplit(x, "\\|")[[1]][length(strsplit(x, "\\|")[[1]])])

  # Get core
  data_arg <- data_arg[data_arg$proportion >= 70,]

  p <- createCircularGraph(data_arg, by_cohort = TRUE, hjust_val = hjust_val)
  return(p)
}
