# ===============================
# List of experiments for mVsig4 CDR3
exp_list <- c("R1_25nM", "R2_25nM", "R2_5nM", "R3_5nM", "R3_0.5nM", "R4_0.5nM")
# ===============================

setwd("~/Desktop/R code/Vsig4NGS/CDR23/For github")

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
pdf("All_mVsig4_CDR3_plots.pdf", width = 8, height = 6)

for (exp_name in exp_list) {
  
  # Read CSV
  file_name <- paste0("m", exp_name, "_CDR3.csv")
  if (!file.exists(file_name)) {
    warning("File not found: ", file_name)
    next
  }
  
  data <- read.csv(file_name)
  
  # Apply color and shape mappings (mVsig4 CDR3-specific Hits)
  data$color_mapping <- sapply(data$ID, get_color)
  data$color_mapping <- factor(data$color_mapping, levels = color_levels)
  
  data$shape_mapping <- ifelse(data$ID %in% c("Rob_6731", "Rob_6732", "Vsig4_WT"),
                               "Hits", "others")
  data$shape_mapping <- factor(data$shape_mapping, levels = shape_levels)
  
  # Label mapping (mVsig4 CDR3-specific)
  data$label_mapping <- ifelse(data$ID == "Rob_6731", "VHH_m4",
                               ifelse(data$ID == "Rob_6732", "VHH_m2",
                                      ifelse(data$ID == "Vsig4_WT", "VHH_WT", "")))
  
  # Plot
  plot <- ggplot(data, aes(x = Enrichment_factor, y = Count_binding)) +
    geom_point(aes(color = color_mapping, shape = shape_mapping), size = 1, alpha = 1) +
    geom_text_repel(aes(label = label_mapping), max.overlaps = Inf) +
    scale_shape_manual(values = shape_values) +
    scale_color_manual(values = color_values) +
    labs(color = "Point Color",
         title = paste("Mouse Vsig4", exp_name, "enriched CDR3")) +
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
