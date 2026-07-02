#main script to run complete source area construction
rm(list = ls())

# set the amount of cores you want to use for parallel computing
NODES <- 20

source("00_data_prep.R")

source("01_Preprocessing_LS.R")

source("02_PixelFlow_parallel_noOptim_LS.R")