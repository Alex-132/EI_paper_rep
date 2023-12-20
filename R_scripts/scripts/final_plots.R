setwd("R_scripts/scripts")
rm(list = ls())

library(corrplot)
library(RColorBrewer)
library(devtools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(writexl)
library(ggimage)
library(ggpubr)
library(scales)
library(gridExtra)
library(cowplot)
library(segmented)


#load data
all_data <- read_xlsx("ddPCR_biomon_mb_data.xlsx")

all_data_na <- read_xlsx("ddPCR_biomon_mb_data_na.xlsx")
gc_data <- read_xlsx("ddPCR_meta_GC.xlsx")
temp <- read_xlsx("../data/Temp_time_MUC.xlsx")


# set datasets
# set data
all_data$relative_abund_sp2 <- ifelse(is.na(all_data_na$relative_abund_sp), NA, all_data$relative_abund_sp)
all_data$relative_abund2 <- ifelse(is.na(all_data_na$relative_abund), NA, all_data$relative_abund)


# Replace values in the MUC column
all_data <- all_data %>%
  mutate(
    MUC = case_when(
      MUC == "EMB262_12_2_MUC" ~ "GOF",
      MUC == "EMB262_6_28_MUC" ~ "EGB",
      MUC == "EMB262_13_8_MUC" ~ "LD",
      TRUE ~ MUC
    )
  )

all_data_na <- all_data_na %>%
  mutate(
    MUC = case_when(
      MUC == "EMB262_12_2_MUC" ~ "GOF",
      MUC == "EMB262_6_28_MUC" ~ "EGB",
      MUC == "EMB262_13_8_MUC" ~ "LD",
      TRUE ~ MUC
    )
  )

gc_data <- gc_data %>%
  mutate(
    GC = case_when(
      GC == "EMB262_12_2_GC" ~ "GOF",
      GC == "EMB262_6_30_GC" ~ "EGB",
      TRUE ~ GC
    )
  )


# for Skeletonema marinoi
Sm <- all_data %>% filter(Interest=="Skeletonema")


#for Apocalathium malmogiense
Am <- all_data %>% filter(Interest=="Apocalathium")


### plot biomonitoring, ddPCR and metabarcoding data ###
## create data frame with realitve abundances
all_data_rel <- all_data %>%
  mutate(
    conc = (conc / sum(conc, na.rm = TRUE)) * 100,
    relative_abund = (relative_abund / sum(relative_abund, na.rm = TRUE)) * 100,
    relative_abund_sp = (relative_abund_sp / sum(relative_abund_sp, na.rm = TRUE)) * 100,
    norm.reads = (norm.reads / sum(norm.reads, na.rm = TRUE)) * 100,
    relative_abund_sp2 = (relative_abund_sp2 / sum(relative_abund_sp2, na.rm = TRUE)) * 100,
    relative_abund2 = (relative_abund2 / sum(relative_abund2, na.rm = TRUE)) * 100
  )


# Filter out data points for relative_abund when relative_abund2 shows NAs
all_data_rel$relative_abund[is.na(all_data$relative_abund2)] <- NA

# Filter out data points for relative_abund_sp when relative_abund_sp2 shows NAs
all_data_rel$relative_abund_sp[is.na(all_data_rel$relative_abund_sp2)] <- NA

# Reshape the data into long format
all_data_rel_long <- all_data_rel %>%
  pivot_longer(
    cols = c(conc, relative_abund_sp, relative_abund, norm.reads),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = case_when(
      variable == "conc" ~ "ddPCR",
      variable == "relative_abund_sp" ~ "biomontoring species",
      variable == "relative_abund" ~ "biomontoring phylum",
      variable == "norm.reads" ~ "metabarcoding"
    )
  )

# for Skeletonema marinoi
Sm_rel_long <- all_data_rel_long %>% filter(Interest=="Skeletonema")


#for Apocalathium malmogiense
Am_rel_long <- all_data_rel_long %>% filter(Interest=="Apocalathium")


#### do the same for 5-year intervals ####
b_all_data_rel <- all_data_rel %>%
  mutate(interval_years = cut(Year, breaks = seq(min(Year), max(Year) + 5, by = 5), right = FALSE)) %>%
  group_by(MUC, interval_years, Interest) %>%
  summarise(Year = first(Year),
            conc = ifelse(all(is.na(conc)), NA_real_, sum(conc, na.rm = TRUE)),
            relative_abund = ifelse(all(is.na(relative_abund)), NA_real_, sum(relative_abund, na.rm = TRUE)),
            relative_abund_sp = ifelse(all(is.na(relative_abund_sp)), NA_real_, sum(relative_abund_sp, na.rm = TRUE)),
            norm.reads = ifelse(all(is.na(norm.reads)), NA_real_, sum(norm.reads, na.rm = TRUE))) %>%
  mutate(interval_years = gsub("\\[|\\)", "", interval_years),
         interval_years = gsub(",", "-", interval_years))

# Reshape the data into long format
b_all_data_rel_long <- b_all_data_rel %>%
  pivot_longer(
    cols = c(conc, relative_abund_sp, relative_abund, norm.reads),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = case_when(
      variable == "conc" ~ "ddPCR",
      variable == "relative_abund_sp" ~ "biomontoring species",
      variable == "relative_abund" ~ "biomontoring phylum",
      variable == "norm.reads" ~ "metabarcoding"
    )
  )

# for Skeletonema marinoi
b_Sm_rel_long <- b_all_data_rel_long %>% filter(Interest=="Skeletonema")


#for Apocalathium malmogiense
b_Am_rel_long <- all_data_rel_long %>% filter(Interest=="Apocalathium")

# Skeletonema
Fig_2a_1 <- ggplot(b_Sm_rel_long, aes(x = Year, y = sqrt(value), color = variable)) +
  geom_point() +
  geom_line() +
  scale_color_manual(
    values = c(
      "ddPCR" = "#694705",
      "biomontoring species" = "#00CDCD",
      "biomontoring phylum" = "#0000FF",
      "metabarcoding" = "salmon"
    )
  ) +
  facet_wrap(~MUC, ncol=1) +
  ggtitle(expression(atop(paste("Relative abundance of ", italic("S. marinoi")), ""))) +
  ylab("square.root(relative abundance)") +
  xlab("Year CE") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size=16),
    axis.title.y = element_text(size=18),
    axis.title.x = element_text(size=18),
    axis.text =element_text(size=14)
  ) +
  labs(color = "")


## Apocalathium
Fig_2b_1 <- ggplot(b_Am_rel_long, aes(x = Year, y = sqrt(value), color = variable)) +
  geom_point() +
  geom_line() +
  scale_color_manual(
    values = c(
      "ddPCR" = "#694705",
      "biomontoring species" = "#00CDCD",
      "biomontoring phylum" = "#0000FF",
      "metabarcoding" = "salmon"
    )
  ) +
  facet_wrap(~MUC, ncol=1) +
  ggtitle(expression(atop(paste("Relative abundance of ", italic("A. malmogiense")), ""))) +
  ylab("square.root(relative abundance)") +
  xlab("Year CE") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size=16),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=18),
    axis.text =element_text(size=14)
  ) +
  labs(color = "")


# combine both plots
Fig_2 <- ggarrange(Fig_2a_1, Fig_2b_1, labels=c("A","B"), ncol = 2, common.legend=TRUE)

tiff(file = "Results_paper/suppl/Fig_S2.tiff", units="in", width = 10, height = 7, res =300, compression = "lzw")
Fig_2
dev.off()

png(file = "Results_paper/suppl/Fig_S2.png", units="in", width = 10, height = 7, res =300)
Fig_2
dev.off()

###### correlation analyses ######

# create 5-year intervals
c_all_data <- all_data_na %>%
  mutate(interval_years = cut(Year, breaks = seq(min(Year), max(Year) + 5, by = 5), right = FALSE)) %>%
  group_by(MUC, interval_years, Interest) %>%
  summarise(Year = first(Year),
            conc = ifelse(all(is.na(conc)), NA_real_, sum(conc, na.rm = TRUE)),
            relative_abund = ifelse(all(is.na(relative_abund)), NA_real_, sum(relative_abund, na.rm = TRUE)),
            relative_abund_sp = ifelse(all(is.na(relative_abund_sp)), NA_real_, sum(relative_abund_sp, na.rm = TRUE)),
            norm.reads = ifelse(all(is.na(norm.reads)), NA_real_, sum(norm.reads, na.rm = TRUE))) %>%
  mutate(interval_years = gsub("\\[|\\)", "", interval_years),
         interval_years = gsub(",", "-", interval_years))

# only keep years that when monitoring was performed
filtered <- c_all_data %>% filter(!is.na(relative_abund_sp))
filtered <- filtered %>% filter(!is.na(conc))

# only NAs in norm.reads stay. replace.
filtered[is.na(filtered)] <- 0


head(filtered)


# Define a function to remove outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Remove outliers from the data
filtered_no_outliers <- filtered %>%
  mutate(
    conc = remove_outliers(conc),
    relative_abund_sp = remove_outliers(relative_abund_sp),
    relative_abund = remove_outliers(relative_abund),
    norm.reads = remove_outliers(norm.reads)
  )

# Calculate the correlation coefficients and p-values for each group of Interest
correlation <- filtered_no_outliers %>%
  group_by(Interest) %>%
  summarize(
    cor_conc_sp = cor(conc, relative_abund_sp, use = "complete.obs"),
    cor_conc_phylum = cor(conc, relative_abund, use = "complete.obs"),
    cor_norm_sp = cor(norm.reads, relative_abund_sp, use = "complete.obs"),
    cor_norm_phylum = cor(norm.reads, relative_abund, use = "complete.obs"),
    cor_conc_norm = cor(conc, norm.reads, use = "complete.obs"),
    cor_phylum_sp = cor(relative_abund, relative_abund_sp, use= "complete.obs"),
    p_conc_sp = cor.test(conc, relative_abund_sp, method = "pearson")$p.value,
    p_conc_phylum = cor.test(conc, relative_abund, method = "pearson")$p.value,
    p_norm_sp = cor.test(norm.reads, relative_abund_sp, method = "pearson")$p.value,
    p_norm_phylum = cor.test(norm.reads, relative_abund, method = "pearson")$p.value,
    p_conc_norm = cor.test(conc, norm.reads, method = "pearson")$p.value,
    p_phylum_sp = cor.test(relative_abund, relative_abund_sp, method = "pearson")$p.value
  )

#### corrplot ####
# Split the data by Interest
data_split <- split(filtered_no_outliers[,c("conc","relative_abund","relative_abund_sp","norm.reads")], filtered_no_outliers$Interest)


# Create a corrplot for each Interest and save it as a .tiff file in the specified output folder
output_folder <- "Results_paper/final/"
for (i in seq_along(data_split)) {
  # Open a new graphics device with increased size
  tiff(file = paste0(output_folder, "corrplot_", names(data_split)[i], ".tiff"), width = 6, height = 6, units = 'in', res = 300, compression = "lzw")
  
  # Increase the top margin of the plot
  par(mar=c(5,4,4,2)+0.1)
  
  # Create the corrplot with a title and reduced font size
  cat("\nCorrplot for Interest:", names(data_split)[i], "\n")
  corrplot(cor(data_split[[i]]), method="circle", col=brewer.pal(n=8, name="RdYlBu"), tl.pos="d", tl.cex=0.6, addCoef.col="black", number.cex=0.8, type="lower")
  
  # Add a title manually
  title(main=names(data_split)[i], line=+0.4)
  
  # Close the graphics device
  dev.off()
}


### plot adjustments ###
# Define the colors and labels
colors <- c("#005f69", "#9C6B30")
becco <- c("#00688B","#FF3030", "#CDAD00")
# Create labels
labels <- c(expression(italic("A. malmogiense")), expression(italic("S. marinoi")))
# Create a named vector that maps the group names to the color values
color_map <- setNames(colors, correlation$Interest)


# Define a custom cube root transformation
cuberoot_trans <- function() {
  trans_new(
    name = "cuberoot",
    transform = function(x) sign(x) * abs(x)^(1/3),
    inverse = function(x) x^3
  )
}



#### correlation of abundance with years ####
correlation_years <- c_all_data %>%
  group_by(Interest) %>%
  summarize(
    cor_conc_year = cor(conc, Year, use = "complete.obs"),
    p_conc_year = cor.test(conc, Year, method = "pearson")$p.value
  )


# plot
# calculate the maximum x and y values of the plot area to position the text labels in the upper right corner
max_x <- max(c_all_data$Year, na.rm = TRUE)
max_y <- max(c_all_data$conc, na.rm = TRUE)

# calculate the distance between the columns
column_distance1 <- 0.2 * (max_x - min(c_all_data$Year, na.rm = TRUE))
# decrease max_x to move text labels towards the left
max_x <- max_x - 0.4 * (max_x - min(c_all_data$Year, na.rm = TRUE))
max_y <- max_y - 0.2 * (max_y - min(c_all_data$conc, na.rm = TRUE))

# create a data frame with the text labels and their positions
text_labels <- data.frame(
  Interest = correlation_years$Interest,
  x = ifelse(correlation_years$Interest == "Apocalathium", max_x, max_x + column_distance1),
  y = max_y - seq(0, by = 0.0 * (max_y - min(c_all_data$conc, na.rm = TRUE)), length.out = nrow(correlation_years)),
  label = paste0("R=", round(correlation_years$cor_conc_year,2), "\np=", round(correlation_years$p_conc_year,4))
)

# Split the data by 'Interest'
split_data <- split(c_all_data, c_all_data$Interest)

# Calculate correlation for each split part
correlation_results <- lapply(split_data, function(df) {
  lm_model <- lm(conc ~ Year, data = df)
  seg_model <- segmented(lm_model)
  bp <- seg_model$psi[,1]
  
  list(
    before = cor.test(df$Year[df$Year < bp], df$conc[df$Year < bp]),
    after = cor.test(df$Year[df$Year >= bp], df$conc[df$Year >= bp])
  )
})

# Extract correlation and p-value for each part
correlation_summary <- lapply(correlation_results, function(ct) {
  c(before_correlation = ct$before$estimate, before_p.value = ct$before$p.value,
    after_correlation = ct$after$estimate, after_p.value = ct$after$p.value)
})

# Convert list to data frame
correlation_split <- do.call(rbind, correlation_summary)

# Ensure correlation_split is always a data frame
if(is.vector(correlation_split)) {
  correlation_split <- data.frame(t(correlation_split))
  names(correlation_split) <- c("before_correlation", "before_p.value", "after_correlation", "after_p.value")
}

# Print the correlation for each 'Interest' group with p-value
print(correlation_split)



# Fit a segmented regression for each group and create a data frame for the segmented lines
seg_data <- do.call(rbind, lapply(split_data, function(df) {
  lm_model <- lm(conc ~ Year, data = df)
  seg_model <- segmented(lm_model)
  bp <- seg_model$psi[,1]
  data.frame(
    Year = c(min(df$Year), bp, max(df$Year)),
    conc = predict(seg_model, newdata = data.frame(Year = c(min(df$Year), bp, max(df$Year)))),
    Interest = unique(df$Interest)
  )
}))

# Create your ggplot object without the 'add = "reg.line"' argument
Fig_3a <- ggscatter(c_all_data, x = "Year", y = "conc", color = "Interest", palette = "colors",
                    conf.int = FALSE) +
  scale_y_continuous(trans = cuberoot_trans()) +
  labs(y = "cube.root(ddPCR concentration)",
       x = "Year CE",
       title = "ddPCR conc. against Years") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.title=element_text(size=18),
        plot.title=element_text(size=20),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 0.5), axis.text.y = element_text(size = 16),
        legend.text = element_text(size = 18)) +
  scale_color_manual(values=colors, labels=labels, name=NULL)


# add the text labels using the annotate function
for (i in seq_len(nrow(text_labels))) {
  Fig_3a <- Fig_3a + annotate("text", x=text_labels$x[i], y=text_labels$y[i], label=text_labels$label[i], color=color_map[text_labels$Interest[i]], size=5, hjust=0)
}


# Add the segmented lines to your plot
Fig_3a <- Fig_3a + geom_line(data = seg_data, aes(x = Year, y = conc, color = Interest), linewidth=1)
Fig_3a

## 5 year intervals
# build index with 5 year intervals
c_all_data <- c_all_data %>%
  group_by(MUC, Year) %>%
  mutate(index = sum(conc[Interest == "Skeletonema"]) / sum(conc[Interest %in% c("Skeletonema", "Apocalathium")])) %>%
  ungroup()

Fig_3c <-  ggplot(c_all_data, mapping = aes(x=Year, y=index, colour=MUC)) +
  geom_point(size=2) +
  geom_line() +
  scale_colour_manual(values = becco) +
  ylim(0,1) +
  xlim(1790, 2020) +
  #geom_smooth(method = lm, se=FALSE) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "#00a8bb", size=1.5) +
  #annotate("text", x=2017, y=0.6, label="GES", color="#00a8bb") +
  facet_wrap(~MUC, nrow=3) +
  labs(x = "Year CE",
       y = "ske/apo-index",
       title="ddPCR ske/apo-index over time",
       subtitle = "5-year intervals") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title=element_text(size=18),
    axis.text = element_text(size=16),
    plot.title=element_text(size=20),
    plot.subtitle = element_text(size=16),
    legend.text = element_text(size = 16),
    legend.position = "non")


# Dia & Dino biomonitoring 
# 5-year intervals
c_all_data <- c_all_data %>%
  group_by(MUC, Year) %>%
  mutate(index_bio_sp = sum(relative_abund_sp[Interest == "Skeletonema"]) / sum(relative_abund_sp[Interest %in% c("Skeletonema", "Apocalathium")])) %>%
  ungroup()

c_all_data <- c_all_data %>%
  group_by(MUC, Year) %>%
  mutate(index_bio_phy = sum(relative_abund[Interest == "Skeletonema"]) / sum(relative_abund[Interest %in% c("Skeletonema", "Apocalathium")])) %>%
  ungroup()

Fig_5b <-  ggplot(c_all_data, mapping = aes(x=Year, y=index_bio_phy, colour=MUC)) +
  geom_point(size=2) +
  geom_line() +
  scale_colour_manual(values = becco) +
  ylim(0,1) +
  xlim(1970, 2020) +
  #geom_smooth(method = lm, se=FALSE) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "#00a8bb", size=1.5) +
  #annotate("text", x=2017, y=0.6, label="GES", color="#00a8bb") +
  facet_wrap(~MUC, nrow=3) +
  labs(x = "Year CE",
       y = "diatom/dinoflagellate-index",
       title="Biomonitoring diatom/dinoflagellate-index over time",
       subtitle = "5-year intervals") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.title=element_text(size=18),
    axis.text = element_text(size=16),
    legend.text = element_text(size = 16),
    plot.title=element_text(size=20),
    plot.subtitle = element_text(size=16),
    legend.position = "non")


Fig_3_3 <- ggarrange(Fig_3a, Fig_3c, Fig_5b, labels=c("A","B", "C"), ncol = 3)

tiff(file = "Results_paper/final/Fig_4.tiff", units="in", width = 18, height = 7, res =300, compression = "lzw")
Fig_3_3
dev.off()

png(file = "Results_paper/final/Fig_4.png", units="in", width = 18, height = 7, res =300)
Fig_3_3
dev.off()




############################
### GC analyses ###
# Create a sequence of breaks that are evenly spaced across the range of your data
breaks <- seq(min(gc_data$numeric_Year), max(gc_data$numeric_Year), length.out = 10)

# Match these breaks to the closest values in numeric_Year
closest_breaks <- sapply(breaks, function(x) gc_data$numeric_Year[which.min(abs(gc_data$numeric_Year - x))])

# Get the corresponding labels for these breaks
label <- gc_data$Year[match(closest_breaks, gc_data$numeric_Year)]

Fig_dis_GC <- ggplot(gc_data, aes(x=numeric_Year, y=sqrt(sqrt(conc)), colour=Interest)) +
  geom_point(size=1) +
  geom_line() +
  scale_colour_manual(values = colors, labels=labels, name=NULL) +
  scale_x_continuous(breaks = closest_breaks, labels = label) + 
  facet_wrap(~GC, ncol=1) +
  labs(x = "Year CE",
       y = "double.square.root\n(average ddPCR concentration)",
       title = "ddPCR concentratiion through the core") +
  theme_bw() + theme(legend.position="bottom",
                     strip.text = element_text(size = 12, face = "bold"),
                     axis.title=element_text(size=18),
                     legend.text = element_text(size = 16),
                     axis.text.x = element_text(size= 14, angle = 90),
                     axis.text.y = element_text(size= 14)) 



### GC index ###

Fig_3_GC <- ggplot(gc_data, mapping = aes(x=numeric_Year, y=index, colour=GC)) +
  geom_point(size=2) +
  geom_line() +
  scale_colour_manual(values = becco, name=NULL) +
  scale_x_continuous(breaks = closest_breaks, labels = label) + 
  ylim(0,1) +
  #geom_smooth(method = lm, se=FALSE) +
  #geom_hline(yintercept=0.5, linetype="dashed", color = "#00a8bb", size=1.5) +
  #annotate("text", x=2017, y=0.6, label="GES", color="#00a8bb") +
  facet_wrap(~GC, nrow=3) +
  labs(x = "Year CE",
       y = "ske/ske-index",
       title="ddPCR ske/apo-index over time") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    axis.title=element_text(size=18),
    legend.text = element_text(size = 16),
    axis.text.x = element_text(size= 14, angle = 90),
    axis.text.y = element_text(size= 14))


## save for paper
Fig_GC <- ggarrange(Fig_dis_GC, Fig_3_GC, labels=c("A","B"), ncol = 2)
Fig_GC <-  annotate_figure(Fig_GC, top = text_grob("GC data ddPCR concentration and index", size = 14, face = "bold"))

tiff(file = "Results_paper/final/Fig_5.tiff", units="in", width = 12, height = 6, res =300, compression = "lzw")
Fig_GC
dev.off()

png(file = "Results_paper/final/Fig_5.png", units="in", width = 12, height = 6, res =300)
Fig_GC
dev.off()



#########################
#### Supplement ####

#correlation between ddPCR concenctrations
cor_species <- cor.test(Am$conc, Sm$conc, method = "pearson")

p_value5 <- format(cor_species$p.value, digits = 2)
r_value5 <- format(cor_species$estimate, digits = 2)

# Create subsets for 'Apocalathium' and 'Skeletonema'
apocalathium_data <- subset(all_data, Interest == "Apocalathium")
skeletonema_data <- subset(all_data, Interest == "Skeletonema")


# Create a new data frame with 'conc' from both subsets
plot_data <- data.frame(Apocalathium = apocalathium_data$conc, Skeletonema = skeletonema_data$conc)

# Create the scatter plot
Fig_8 <- ggplot(plot_data, aes(x = Apocalathium, y = Skeletonema)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "#8B0000", size=1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Apocalathium ddPCR conc.", y = "Skeletonema ddPCR conc.", 
       title = "Concentration of Apocalathium vs Skeletonema") +
  theme_bw() +
  geom_text(x = 3, y = 2, hjust = 1.1, vjust = 2,
            label = paste("R =", r_value5, "\np =", p_value5))


tiff(file = "Results_paper/suppl/Fig_S4.tiff", units="in", width = 6, height = 6, res =300, compression = "lzw")
Fig_8
dev.off()

png(file = "Results_paper/suppl/Fig_S4.png", units="in", width = 6, height = 6, res =300)
Fig_8
dev.off()

# 5-year intervals
correlation_indices3 <- cor.test(c_all_data$index, c_all_data$index_bio_sp, method = "pearson")
p_value3 <- format(correlation_indices3$p.value, digits = 2)
r_value3 <- format(correlation_indices3$estimate, digits = 2)

correlation_indices4 <- cor.test(c_all_data$index, c_all_data$index_bio_phy, method = "pearson")
p_value4 <- format(correlation_indices4$p.value, digits = 2)
r_value4 <- format(correlation_indices4$estimate, digits = 2)

# Skeletonema & Apocalathium

# Fit an exponential model
exp_model <- nls(index ~ a * exp(b * index_bio_sp), data = c_all_data, start = list(a = 0.5, b = 0.5), control = nls.control(maxiter = 300))


# log transformation
ln_model <- lm(log(index) ~ index_bio_sp, data = c_all_data)


# plot
Fig_7a <- ggscatter(c_all_data, y = "index", x = "index_bio_sp") +
  geom_line(aes(x = index_bio_sp, y = exp(predict(ln_model, newdata = data.frame(index_bio_sp = index_bio_sp)))), color = "#8B0000", size=1) +
  scale_x_continuous(trans = cuberoot_trans()) +
  scale_y_continuous(trans = cuberoot_trans()) +
  labs(x = "Index biomonitoring species",
       y = "Index ddPCR",
       title = "ddPCR ske/apo-index. against biomonitoring ske/apo-index",
       subtitle = "5-year intervals") +
  theme_bw() +
  theme(legend.position="bottom",
        axis.title=element_text(size=14),
        plot.title=element_text(size=14))+
  scale_color_manual(values=colors, labels=labels, name=NULL)+
  geom_text(x = 0.95, y = 1.1, hjust = 1.1, vjust = 2,
            label = paste("R =", r_value3, "\np =", p_value3))


tiff(file = "Results_paper/suppl/Fig_S5.tiff", units="in", width = 7, height = 6, res =300, compression = "lzw")
Fig_7a
dev.off()

png(file = "Results_paper/suppl/Fig_S5.png", units="in", width = 7, height = 6, res =300)
Fig_7a
dev.off()



