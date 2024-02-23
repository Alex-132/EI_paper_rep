# Data preparation and Downstream analyses

This repository contains several R scripts and data used for the data preparation, analysis, and visualization. Here's a brief overview of each script and their dependencies:

## 1. Data Preparation (`data_prep.R`)

This script is used to prepare the data and combine everything into one matrix.

## 2. Data Analysis and Visualization (`replicate_consistency.R`,`final_plots.R`)

This scripts are used to perform all analyses run and used in the paper and supplement.

## 3. Paleoenvironmental Data Analysis (`env_data_publ.R`)

This script is used to analyze and plot the paleoenvironmental data.

## 4. Metabarcoding Data Cleaning (`Script_Euka02_cleaning.Rmd`)

This script is used to remove the metabarcoding data set from PCR error, possible contamination and similar.

## 5. Metabarcoding Data Preparation (`Scripts_MUC__metabarcoding.R`)

This script is used to combine the prepare the metabarcoding data for the final data analysis.

## 6. Biomonitoring Data Preparation on Taxon-Level (`Script_MUCs_dia_dino__ICES.Rmd`)

This script is used to prepare and reduce the ICES monitoring data for the final data analysis on taxon-level (Diatom & Dinoflagellates)

## 6. Biomonitoring Data Preparation on Genus-Level (`Script_MUCs_skele_apo__ICES.Rmd`)

This script is used to prepare and reduce the ICES monitoring data for the final data analysis on genus-level (Apocalathium & Skeletonema)

## Dependencies

The following R packages are required to run the scripts:

```R
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
library(tidyverse)
library(reshape)
library(viridis)
library(seqinr)
library(worrms)
library(hrbrthemes)
library(ggthemes)
```

Please ensure that you have these packages installed before running the scripts. You can install any missing packages using `install.packages("package-name")` in R.
