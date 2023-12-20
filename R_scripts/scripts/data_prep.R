setwd("C:/Users/alex9/Documents/PhD/R_stuff/ddPCR_screening")
rm(list = ls())


library(devtools)
library(tidyverse)
library(readxl)
library(writexl)

### load data ###
#load MUC data
muc_ages <- read_xlsx("core_sampling_data_MUC_ages.xlsx")

#load GC data
gc_ages <- read_xlsx("core_sampling_data_GC_ages.xlsx")
# load ddPCR results
ddpcr <- read_xlsx("ddPCR_abundances.xlsx")

# load biomonitoring data for Skeletonema and Apocalathium
biomon_march_to_june <- read.csv("../biomonitoring/Biomonitoring_MUCs_march-to-june__JR.csv")
biomon_march_to_may <- read.csv("../biomonitoring/Biomonitoring_MUCs_march-to-may__JR.csv")
biomon_april_to_may <- read.csv("../biomonitoring/Biomonitoring_MUCs_april-may__JR.csv")

# load biomonitoring data for Diatoms and Dinoflagellates
biomon_march_to_may_dia_dino <- read.csv("../biomonitoring/Biomonitoring_MUCs_march-to-may__dino-dia__JR.csv")

# load metabarcoding
mb_sp <- read.csv("Data_Skeletonema_Apocalathium__all_cores__07june2023.csv")
mb_dia_dino <- read.csv("../biomonitoring/sedaDNA_MUCs_metabarcoding_Euka02__dino-dia__JR.csv")


### adjust data ###

# muc age table combine with ddPCR data
muc_ages %>%
  full_join(ddpcr, by=c("tag")) %>%
  dplyr::select("cruise", "tag", "sample_id", "relative\r\nconc. FAM", "X̅ conc. Am\r\n[copies/µL]", "relative\r\nconc. HEX", "X̅ conc. Sm\r\n [copies/µL]", "depth","age", "latitude [°N]", "longditude [°E]", "water depth [m bsl]","delta_Am", "delta_Sm", "dia/dino-index") -> muc_ages

names(muc_ages) <- c("cruise", "tag", "sample_id", "Am_conc", "average_Am_conc", "Sm_conc", "average_Sm_conc", "depth", "age", "latitude [°N]", "longditude [°E]", "water depth [m bsl]", "delta_Am", "delta_Sm", "index")
muc_ages <- muc_ages[- grep("Mix", muc_ages$cruise),]
muc_ages <- muc_ages[- grep("NTC", muc_ages$cruise),]
muc_ages <- muc_ages[- grep("EB", muc_ages$cruise),]
muc_ages %>%
  drop_na(cruise) %>%
  drop_na(depth) -> muc_ages
muc_ages <- unique(muc_ages)

muc_age_new <- muc_ages[- grep("EMB262_3_10_MUC", muc_ages$cruise),]
muc_age_new <- muc_age_new %>%
  rename(Apocalathium_malmogiense = average_Am_conc,
         Skeletonema_marinoi = average_Sm_conc)

muc_age_new <- muc_age_new %>%
  pivot_longer(cols = c(Apocalathium_malmogiense, Skeletonema_marinoi),
               names_to = "species",
               values_to = "conc")
muc_age_new %>%
  dplyr::select("cruise", "age", "conc", "species", "index") -> muc_age_new
# Create a new column in the first data frame that maps the values in the species column to the corresponding values in the Interest column
muc_age_new <- muc_age_new %>% dplyr::mutate(Interest = case_when(
  species == "Apocalathium_malmogiense" ~ "Apocalathium",
  species == "Skeletonema_marinoi" ~ "Skeletonema",
  TRUE ~ as.character(species)
))
muc_age_new <- muc_age_new %>% dplyr::select(-species)
# rename the columns of ddPCR dataframe to match the column names in the biomonitoring dataframes
muc_age_new <- muc_age_new %>% rename(MUC = cruise,  Year = age)
# Merge the first and second data frames
muc_age_new <- muc_age_new %>%
  mutate(phylum = case_when(
    Interest == "Apocalathium" ~ "Dinoflagellates",
    Interest == "Skeletonema" ~ "Diatoms",
    TRUE ~ NA_character_
  ))

muc_age_new <- unique(muc_age_new)

# GC age table combine with ddPCR data
# get the tag right, so it fits the ddPCR data
gc_ages$tag <- gsub("EMB262_12_3_GC_", "EMB262_12_2_GC_", gc_ages$tag)
gc_ages %>%
  full_join(ddpcr, by=c("tag")) %>%
  dplyr::select("cruise", "tag", "sample_id", "relative\r\nconc. FAM", "X̅ conc. Am\r\n[copies/µL]", "relative\r\nconc. HEX", "X̅ conc. Sm\r\n [copies/µL]", "depth","age(BP)", "age","numeric_age(BP)", "numeric_age", "latitude [°N]", "longditude [°E]", "water depth [m bsl]","delta_Am", "delta_Sm", "dia/dino-index") -> gc_ages

names(gc_ages) <- c("cruise", "tag", "sample_id", "Am_conc", "average_Am_conc", "Sm_conc", "average_Sm_conc", "depth", "age(BP)", "age" ,"numeric_age(BP)", "numeric_age", "latitude [°N]", "longditude [°E]", "water depth [m bsl]", "delta_Am", "delta_Sm", "index")
gc_ages <- gc_ages[- grep("Mix", gc_ages$cruise),]
gc_ages <- gc_ages[- grep("NTC", gc_ages$cruise),]
gc_ages <- gc_ages[- grep("EB", gc_ages$cruise),]
gc_ages %>%
  drop_na(cruise) %>%
  drop_na(depth) -> gc_ages
gc_ages <- unique(gc_ages)

gc_age_new <- gc_ages%>%
  rename(Apocalathium_malmogiense = average_Am_conc,
         Skeletonema_marinoi = average_Sm_conc)


gc_age_new <- gc_age_new %>%
  pivot_longer(cols = c(Apocalathium_malmogiense, Skeletonema_marinoi),
               names_to = "species",
               values_to = "conc")
gc_age_new %>%
  dplyr::select("cruise", "age", "age(BP)", "numeric_age(BP)", "numeric_age", "depth", "conc", "species", "index") -> gc_age_new
# Create a new column in the first data frame that maps the values in the species column to the corresponding values in the Interest column
gc_age_new <- gc_age_new %>% dplyr::mutate(Interest = case_when(
  species == "Apocalathium_malmogiense" ~ "Apocalathium",
  species == "Skeletonema_marinoi" ~ "Skeletonema",
  TRUE ~ as.character(species)
))
gc_age_new <- gc_age_new %>% dplyr::select(-species)
# rename the columns of ddPCR dataframe to match the column names in the biomonitoring dataframes
gc_age_new <- gc_age_new %>% rename(GC = cruise,  Year = age, numeric_Year = numeric_age)
# Merge the first and second data frames
gc_age_new <- gc_age_new %>%
  mutate(phylum = case_when(
    Interest == "Apocalathium" ~ "Dinoflagellates",
    Interest == "Skeletonema" ~ "Diatoms",
    TRUE ~ NA_character_
  ))

gc_age_new <- unique(gc_age_new)

# save
write_xlsx(gc_age_new, "ddPCR_meta_GC.xlsx")


# prepare the biomonitoring data

a_biomon_march_to_may <- biomon_march_to_may %>%
  filter(Interest != "Pericrocotus brevirostris") %>%
  group_by(MUC, Year, Interest) %>%
  dplyr::mutate(total_abund = sum(total.ABUNDNR),
                relative_abund = total_abund / tot.biomass * 100) %>%
  dplyr::select(MUC, Year, Interest, relative_abund) %>%
  dplyr::distinct()

a_biomon_march_to_may <- a_biomon_march_to_may %>%
  mutate(phylum = case_when(
    Interest == "Apocalathium"~"Dinoflagellates",
    Interest == "Skeletonema"~"Diatoms",
    TRUE ~ NA_character_
  ))


# prepare the biomonitoring data dia dino

a_biomon_march_to_may_dia_dino <- biomon_march_to_may_dia_dino %>%
  group_by(MUC, Year, phylum) %>%
  dplyr::mutate(total_abund = sum(total.ABUNDNR),
                relative_abund = total_abund / tot.biomass * 100) %>%
  dplyr::select(MUC, Year, phylum, relative_abund) %>%
  dplyr::distinct()

# adjust phylum columns
a_biomon_march_to_may_dia_dino <- a_biomon_march_to_may_dia_dino %>%
  dplyr::mutate(phylum = case_when(
    phylum == "Bacillariophyta" ~ "Diatoms",
    phylum == "Myzozoa" ~ "Dinoflagellates",
    TRUE ~ phylum
  ))


# prepare the metabarcoding data
a_mb_dia_dino <- mb_dia_dino %>%
  rename(MUC = station,
         Year = age,
         phylum = Interest) %>%
  dplyr::select(MUC, Year, phylum, norm.reads) %>%
  mutate(phylum = case_when(
    phylum == "Diatom" ~ "Diatoms",
    phylum == "Dinoflagellata" ~ "Dinoflagellates",
    TRUE ~ phylum
  ))

# remove interest column from muc_ddPCR and the biomonitoring datasets to join
muc_age_new_1 <- muc_age_new %>%
  dplyr::select(-Interest)

a_biomon_march_to_may <- a_biomon_march_to_may %>%
  ungroup() %>%
  dplyr::select(-Interest) %>%
  rename(relative_abund_sp = relative_abund)
  
# combine all the data in 1 dataset
all_data <- muc_age_new_1 %>%
  full_join(a_biomon_march_to_may, by= c("MUC", "Year", "phylum")) %>%
  full_join(a_biomon_march_to_may_dia_dino, by = c("MUC", "Year", "phylum")) %>%
  full_join(a_mb_dia_dino, by=c("MUC", "Year", "phylum"))

# add interest column
all_data <- all_data %>%
  mutate(Interest = case_when(
    phylum == "Dinoflagellates" ~ "Apocalathium",
    phylum == "Diatoms" ~ "Skeletonema",
    TRUE ~ NA_character_
  ))

write_xlsx(all_data, "ddPCR_biomon_mb_data_na.xlsx")

# replace NAs with 0
all_data[is.na(all_data)] <- 0

# store in excel sheet
write_xlsx(all_data, "ddPCR_biomon_mb_data.xlsx")

