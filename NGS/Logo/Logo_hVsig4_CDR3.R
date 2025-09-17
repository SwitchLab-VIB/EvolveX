# Set up a work directory
setwd("~/Desktop/R code/Vsig4NGS/CDR23/For github/Logo")

# Install required packages
install.packages("ggseqlogo")
install.packages("ggplot2")
install.packages("cowplot")

# Load required packages
library(ggseqlogo)
library(ggplot2)
library(cowplot)  # Load the cowplot package

# Replace 'input_file.csv' with the actual file path
hR1NC_CDR3 <- "hR1NC_CDR3.csv"
WT <- "WT.csv"
hR1_CDR3 <- "hR1_CDR3.csv"
input <-"input.csv"


# Load the CSV file
data_wt <- read.csv(WT)
data_hR1_CDR3 <- read.csv(hR1_CDR3)
data_hR1NC_CDR3 <- read.csv(hR1NC_CDR3)
data_input <- read.csv(input)

csRasmol <- make_col_scheme(chars=c("D","E","C","M","K","R","S","T","F","Y","N","Q","G","L","V","I","A","W","H","P"),
                            cols=c("red","red","yellow","yellow","blue","blue","orange","orange","midnightblue","midnightblue","cyan","cyan","gray50","green3","green3","green3","gray30","pink1","cornflowerblue","magenta"))

# Set the limits for the x-axis to ensure equal width for each amino acid position
x_limits <- c(0.5, max(10, 17.5))  # Adjust based on your data


# Create a sequence logo for CDR3 without amino acid legend
p1 <- ggseqlogo(data_hR1_CDR3$Sequence, method = 'prob',col_scheme = csRasmol) +
  scale_x_continuous(breaks = 1:16, labels = 99:114, limits = x_limits) +
  theme(legend.position = "none",
        axis.title.x = element_text(vjust = -0.5, hjust = 0.45)) +  # Adjust the vjust value as needed
  labs(x = "Position", title = "Human Vsig4 R1 115 enriched designs, CDR3") +
  theme(plot.title = element_text(size = 14), axis.text.x = element_text(size = 8))

p2 <- ggseqlogo(data_hR1NC_CDR3$Sequence, method = 'prob',col_scheme = csRasmol) +
  scale_x_continuous(breaks = 1:16, labels = 99:114, limits = x_limits) +
  theme(legend.position = "none",
        axis.title.x = element_text(vjust = -0.5, hjust = 0.45)) +  # Adjust the vjust value as needed
  labs(x = "Position", title = "Human Vsig4 R1 115 negative designs, CDR3") +
  theme(plot.title = element_text(size = 14), axis.text.x = element_text(size = 8))

p3 <- ggseqlogo(data_wt$CDR3, method = 'prob',col_scheme = csRasmol) +
  scale_x_continuous(breaks = 1:16, labels = 99:114, limits = x_limits) +
  theme(legend.position = "none",
        axis.title.x = element_text(vjust = -0.5, hjust = 0.45)) +  # Adjust the vjust value as needed
  labs(x = "Position", title = "VHH_WT, CDR3") +
  theme(plot.title = element_text(size = 14), axis.text.x = element_text(size = 8))

p4 <- ggseqlogo(data_input$CDR3, method = 'prob',col_scheme = csRasmol) +
  scale_x_continuous(breaks = 1:16, labels = 99:114, limits = x_limits) +
  theme(legend.position = "none",
        axis.title.x = element_text(vjust = -0.5, hjust = 0.45)) +  # Adjust the vjust value as needed
  labs(x = "Position", title = "Input designed sequences, CDR3") +
  theme(plot.title = element_text(size = 14), axis.text.x = element_text(size = 8))

# Arrange the sequence logos vertically on the left and right using cowplot

Human_sequence_logos <- plot_grid(p3, p4, p1, p2, ncol = 1)

# Display the Human sequence logos
print(Human_sequence_logos)

ggsave("Human_sequence_logos.pdf", Human_sequence_logos, dpi = 600, width = 8, height = 6)
