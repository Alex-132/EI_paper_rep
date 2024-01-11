setwd("G:/.shortcut-targets-by-id/1mLh3UzlnjPm_h7mFbD5TlPOG-84Bw4vm/phytoark_google_drive/cores_data_jerome/published environmental data")

# clean environment
rm(list = ls())

library(tidyverse)
library(readxl)
library(writexl)
library(gridExtra)


# load data
sal <- read_xlsx("salinity datasets Holocene Baltic Sea.xlsx", skip=1)
temp <- read_xlsx("temperature datasets Holocene Baltic Sea.xlsx", sheet="Lakes", skip=1)
oxy <- read_xlsx("../../env_factors_history/History_of_the_baltic_sea.xlsx", sheet="Oxygen")


# add Year CE time scale
temp$`Year CE` <- 1950 - temp$`Age (yr cal BP)`
sal$`Year CE` <- 1950 - sal$`year (cal BP)`

# filter oxy values EGB as salinity is also from EGB
oxy <- oxy[oxy$location == "EGB",]


# Rename the columns
names(oxy)[names(oxy) == "Mn/Ti"] <- "Mn_Ti"
names(temp)[names(temp) == "air Temp (°C)"] <- "temp"
names(sal)[names(sal) == "freshwater"] <- "freshwater_sp_proportion"

# Select only the columns you need
oxy <- oxy[, c("Year CE", "Mn_Ti")]
temp <- temp[, c("Year CE", "temp")]
sal <- sal[, c("Year CE", "freshwater_sp_proportion")]

# Combine all data by "Year CE"
env_data <- merge(oxy, temp, by = "Year CE", all = TRUE)
env_data <- merge(env_data, sal, by = "Year CE", all = TRUE)



# Create the plots
# temperature
plot_temp <- ggplot(temp, aes(x = temp, y = `Year CE`)) +
  #geom_point(color = "#4A708B") +
  geom_path(color = "#4A708B") +
  labs(x = "Air temperature \n[°C]", y = "Year CE") +
  scale_y_continuous(limits = c(-10000, 2000), breaks = seq(-10000, 2000, 1000), expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size = 14),  
        axis.title.x = element_text(size = 16))

# salinity
plot_freshwater <- ggplot(sal, aes(x = freshwater_sp_proportion, y = `Year CE`)) +
  #geom_point(color = "#CDAD00") +
  geom_path(color = "#CDAD00") +
  labs(x = "Freshwater species \n[%]") +
  scale_y_continuous(limits = c(-10000, 2000), breaks = seq(-10000, 2000, 1000), expand = c(0, 0)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size = 14),  
        axis.title.x = element_text(size = 16),  
        axis.text.y = element_text(size = 14),  
        axis.title.y = element_text(size = 16))  
# oxygen
plot_Mn_Ti <- ggplot(oxy, aes(x = Mn_Ti, y = `Year CE`)) +
  #geom_point(color = "#556B2F") +
  geom_path(color = "#556B2F") +
  labs(x = "Oxygen \n[Mn/Ti]") +
  scale_y_continuous(limits = c(-10000, 2000), breaks = seq(-10000, 2000, 1000), expand = c(0, 0)) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size = 14),  
        axis.title.x = element_text(size = 16))

# Combine the plots
Fig_1 <- grid.arrange(plot_freshwater, plot_temp, plot_Mn_Ti, ncol = 3)

# Save the plots to a TIFF file
ggsave("C:/Users/alex9/Documents/PhD/R_stuff/environmental_data/literature/Fig_1.tiff", grid.arrange(plot_freshwater, plot_temp, plot_Mn_Ti, ncol = 3), width = 10, height = 7, dpi = 300)

# Save the plots to a PNG file
ggsave("C:/Users/alex9/Documents/PhD/R_stuff/environmental_data/literature/Fig_1.png", grid.arrange(plot_freshwater, plot_temp, plot_Mn_Ti, ncol = 3), width = 10, height = 7, dpi = 300)


# different orientation
# Create the plots
# temperature
plot_temp <- ggplot(temp, aes(y = temp, x = `Year CE`)) +
  geom_line(color = "#4A708B") +
  labs(y = "Air temperature \n[°C]", x = NULL) +  # Remove x-axis title
  scale_x_continuous(limits = c(-10000, 2000), breaks = seq(-10000, 2000, 1000), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
        axis.text.y = element_text(size = 14),  
        axis.title.y = element_text(size = 18),
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points"))

# salinity
plot_freshwater <- ggplot(sal, aes(y = freshwater_sp_proportion, x = `Year CE`)) +
  geom_line(color = "#CDAD00") +
  labs(y = "Freshwater species \n[%]", x = NULL) +  # Remove x-axis title
  scale_x_continuous(limits = c(-10000, 2000), breaks = seq(-10000, 2000, 1000), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
        axis.text.y = element_text(size = 14),  
        axis.title.y = element_text(size = 18),
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points"))

# oxygen
plot_Mn_Ti <- ggplot(oxy, aes(y = Mn_Ti, x = `Year CE`)) +
  geom_line(color = "#556B2F") +
  labs(y = "Oxygen \n[Mn/Ti]") +
  scale_x_continuous(limits = c(-10000, 2000), breaks = seq(-10000, 2000, 1000), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.y = element_line(colour = "black"), axis.line.x = element_line(colour = "black"),
        axis.text.x = element_text(size = 14, angle =90),  
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),  
        axis.title.y = element_text(size = 18),
        plot.margin = margin(5.5, 20, 5.5, 5.5, "points"))

Fig_2 <- grid.arrange(plot_freshwater, plot_temp, plot_Mn_Ti, nrow = 3)

# Save the plots to a TIFF file
ggsave("C:/Users/alex9/Documents/PhD/R_stuff/environmental_data/literature/Fig_2.tiff", grid.arrange(plot_freshwater, plot_temp, plot_Mn_Ti, nrow = 3), width = 12, height = 8, dpi = 300)

# Save the plots to a PNG file
ggsave("C:/Users/alex9/Documents/PhD/R_stuff/environmental_data/literature/Fig_2.png", grid.arrange(plot_freshwater, plot_temp, plot_Mn_Ti, nrow = 3), width = 12, height = 8, dpi = 300)



