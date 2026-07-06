This repository contains code to conduct the first two steps of the PollenPixel workflow (pollen to land cover and land cover to source area).

The pipeline can be executed by just downloading the repository and running the script main.R. Input data and packages will be downloaded and installed automatically. 

Please note that you can adjust the number of nodes to use for parallelization in the file scripts/02_PixelFlow_parallel_noOptim_LS.R. The variable is called nodes and located in line 8.

Results will be stored in a newly created directory "landcoverchange".