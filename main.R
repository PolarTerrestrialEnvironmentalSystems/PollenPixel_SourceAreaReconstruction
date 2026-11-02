#main script to run complete source area construction
rm(list = ls())

setwd(getSrcDirectory(function(){})[1])

source("scripts/00_data_prep.R")

# Set workers after data preparation, which clears the global environment.
NODES <- 5

source("scripts/01_Preprocessing_LS.R")

source("scripts/02_PixelFlow_parallel_noOptim_LS.R")
