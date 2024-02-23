setwd("C:/Users/alex9/Documents/PhD/Publications/ddPCR_microfossils/Supplement")
rm(list = ls())

# Load packags
library(readxl)
library(writexl)
library(tidyverse)
library(gridExtra)

# load data

ddpcr_raw <- read_xlsx("Sample_info_MUC_GC_ddPCR.xlsx", sheet = "ddPCR")

head(ddpcr_raw, n=2)

# Remove rows that include "NTC" and "Extraction Blank"
filtered<- ddpcr_raw %>%
  filter(!grepl("NTC|Extraction Blank", tag))


# Select only the tag and the concentration columns
selected_data <- filtered %>%
  select(tag, replicate,`conc. FAM\r\n[copies/µL]`, `conc. HEX\r\n[copies/µL]`)

# account for samples that were ran twice with replicates
# Add a suffix to the tag for the first and second set of replicates
selected_data <- selected_data %>%
  group_by(tag) %>%
  mutate(replicate_set = ifelse(row_number() <= 3, "_a", "_b"),
         tag = paste(tag, replicate_set, sep = ""))

# account samples that have less thann 3 replicates
selected_data <- selected_data %>%
  group_by(tag) %>%
  filter(n() == 3)

# Reshape the data
reshaped_data <- selected_data %>%
  pivot_wider(names_from = replicate, values_from = c(`conc. FAM\r\n[copies/µL]`, `conc. HEX\r\n[copies/µL]`))

# Separate FAM and HEX data
fam_data <- reshaped_data %>% select(starts_with("conc. FAM"))
hex_data <- reshaped_data %>% select(starts_with("conc. HEX"))

# Calculate correlations and p-values for FAM
fam_cor12 <- cor.test(fam_data$`conc. FAM\r\n[copies/µL]_1`, fam_data$`conc. FAM\r\n[copies/µL]_2`)
fam_cor13 <- cor.test(fam_data$`conc. FAM\r\n[copies/µL]_1`, fam_data$`conc. FAM\r\n[copies/µL]_3`)
fam_cor23 <- cor.test(fam_data$`conc. FAM\r\n[copies/µL]_2`, fam_data$`conc. FAM\r\n[copies/µL]_3`)

# Calculate correlations and p-values for HEX
hex_cor12 <- cor.test(hex_data$`conc. HEX\r\n[copies/µL]_1`, hex_data$`conc. HEX\r\n[copies/µL]_2`)
hex_cor13 <- cor.test(hex_data$`conc. HEX\r\n[copies/µL]_1`, hex_data$`conc. HEX\r\n[copies/µL]_3`)
hex_cor23 <- cor.test(hex_data$`conc. HEX\r\n[copies/µL]_2`, hex_data$`conc. HEX\r\n[copies/µL]_3`)

# Create data frames for FAM and HEX correlations and p-values
fam_table <- data.frame(
  pairs = c("1-2", "1-3", "2-3"),
  correlation = c(fam_cor12$estimate, fam_cor13$estimate, fam_cor23$estimate),
  p_value = c(fam_cor12$p.value, fam_cor13$p.value, fam_cor23$p.value)
)

hex_table <- data.frame(
  pairs = c("1-2", "1-3", "2-3"),
  correlation = c(hex_cor12$estimate, hex_cor13$estimate, hex_cor23$estimate),
  p_value = c(hex_cor12$p.value, hex_cor13$p.value, hex_cor23$p.value)
)

# Print the tables
print(fam_table)
print(hex_table)


# Create a function to plot actual values with double square root transformation and a trend line
plot_values <- function(data, x, y, title) {
  ggplot(data, aes_string(x = paste0("sqrt(sqrt(", x, "))"), y = paste0("sqrt(sqrt(", y, "))"))) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "darkred") +
    labs(title = title, x = x, y = y) +
    theme_minimal()
}


# Create plots for FAM
fam_plot12 <- plot_values(fam_data, "`conc. FAM\r\n[copies/µL]_1`", "`conc. FAM\r\n[copies/µL]_2`", "FAM Replicates 1 vs 2")
fam_plot13 <- plot_values(fam_data, "`conc. FAM\r\n[copies/µL]_1`", "`conc. FAM\r\n[copies/µL]_3`", "FAM Replicates 1 vs 3")
fam_plot23 <- plot_values(fam_data, "`conc. FAM\r\n[copies/µL]_2`", "`conc. FAM\r\n[copies/µL]_3`", "FAM Replicates 2 vs 3")

# Create plots for HEX
hex_plot12 <- plot_values(hex_data, "`conc. HEX\r\n[copies/µL]_1`", "`conc. HEX\r\n[copies/µL]_2`", "HEX Replicates 1 vs 2")
hex_plot13 <- plot_values(hex_data, "`conc. HEX\r\n[copies/µL]_1`", "`conc. HEX\r\n[copies/µL]_3`", "HEX Replicates 1 vs 3")
hex_plot23 <- plot_values(hex_data, "`conc. HEX\r\n[copies/µL]_2`", "`conc. HEX\r\n[copies/µL]_3`", "HEX Replicates 2 vs 3")

# Combine all plots into one figure
tiff(file="Fig_replicate_consistancy.tiff", units="in", width = 10, height = 7, res =300, compression = "lzw")
grid.arrange(fam_plot12, fam_plot13, fam_plot23, hex_plot12, hex_plot13, hex_plot23, ncol = 2)
dev.off()



