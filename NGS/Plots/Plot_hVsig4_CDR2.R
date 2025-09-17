# ===============================
# List of experiments
exp_list <- c("R1_500nM", "R2_500nM", "R2_50nM", "R3_50nM", "R3_5nM", "R4_5nM")
# ===============================

# Set working directory
setwd("~/Desktop/R code/Vsig4NGS/CDR23/For github")

# Load libraries
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)

# Define color mapping function
get_color <- function(id) {
  if (grepl("^Rob_", id)) {
    return("Red")
  } else if (id == "Vsig4_WT") {
    return("Blue")
  } else {
    return("Black")
  }
}

# Fixed color and shape levels
color_levels <- c("Black", "Blue", "Red")
color_values <- c(Black = "black", Blue = "blue", Red = "red")

shape_levels <- c("others", "Hits")
shape_values <- c(others = 1, Hits = 16)  # hollow circle for others, filled for hits

# Open combined PDF
pdf("All_hVsig4_CDR2_plots.pdf", width = 8, height = 6)

for (exp_name in exp_list) {
  
  file_name <- paste0("h", exp_name, "_CDR2.csv")
  if (!file.exists(file_name)) {
    warning("File not found: ", file_name)
    next
  }
  
  data <- read.csv(file_name)
  
  # Apply color and shape mappings
  data$color_mapping <- sapply(data$ID, get_color)
  data$color_mapping <- factor(data$color_mapping, levels = color_levels)
  
  data$shape_mapping <- ifelse(data$ID %in% c("Rob_1633", "Rob_166", "Rob_1005",
                                              "Rob_4071", "Rob_2839", "Vsig4_WT"),
                               "Hits", "others")
  data$shape_mapping <- factor(data$shape_mapping, levels = shape_levels)
  
  # Label mapping
  data$label_mapping <- ifelse(data$ID == "Rob_1633", "VHH_h1",
                               ifelse(data$ID == "Rob_166", "VHH_h2",
                                      ifelse(data$ID == "Rob_1005", "VHH_h3",
                                             ifelse(data$ID == "Rob_4071", "VHH_h4",
                                                    ifelse(data$ID == "Rob_2839", "VHH_h5",
                                                           ifelse(data$ID == "Vsig4_WT", "VHH_WT", ""))))))
  
  # Plot
  plot <- ggplot(data, aes(x = Enrichment_factor, y = Count_binding)) +
    geom_point(aes(color = color_mapping, shape = shape_mapping), size = 1, alpha = 1) +
    geom_text_repel(aes(label = label_mapping), max.overlaps = Inf) +
    scale_shape_manual(values = shape_values) +
    scale_color_manual(values = color_values) +
    labs(color = "Point Color",
         title = paste("Human Vsig4", exp_name, "enriched CDR2")) +
    theme(legend.position = "none",
          axis.text = element_text(size = 16, color = "black", face = "bold"),
          axis.title = element_text(size = 18, face = "bold"),
          axis.line = element_line(size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
          plot.title.position = "plot") +
    scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_x_log10(breaks = c(10, 100, 1000, 10000, 100000),
                  labels = trans_format("log10", math_format(10^.x)))
  
  # Print into the combined PDF
  print(plot)
  
  message("Added plot to combined PDF: ", exp_name)
}

# Close combined PDF
dev.off()
