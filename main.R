#main script to run complete source area construction
rm(list = ls())

setwd(getSrcDirectory(function(){})[1])

# set the amount of cores you want to use for parallel computing
NODES <- 5

source("scripts/00_data_prep.R")

source("scripts/01_Preprocessing_LS.R")

source("scripts/02_PixelFlow_parallel_noOptim_LS.R")