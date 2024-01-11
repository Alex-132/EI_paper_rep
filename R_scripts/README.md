# Data preparation and Downstream analyses

This repository contains several R scripts and data used for the data preparation, analysis, and visualization. Here's a brief overview of each script and their dependencies:

## 1. Data Preparation (`data_prep.R`)

This script is used to prepare the data and combine everything into one matrix.

## 2. Data Analysis and Visualization (`final_plots.R`)

This script is used to perform all analyses run and used in the paper and supplement.

## 3. Paleoenvironmental Data Analysis (`env_data_publ.R`)

This script is used to analyze and plot the paleoenvironmental data.

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
```

Please ensure that you have these packages installed before running the scripts. You can install any missing packages using `install.packages("package-name")` in R.
