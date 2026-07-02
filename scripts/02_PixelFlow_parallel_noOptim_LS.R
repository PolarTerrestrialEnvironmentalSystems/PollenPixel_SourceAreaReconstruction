#################################
### Template: LandCoverChange ###
### 22.10.2023                ###
###  Linux 5/1                ###
#################################
rm(list = ls())

suppressMessages({
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(readr)
  library(readxl)
  library(sf); sf_use_s2(FALSE)
  library(stars)
  library(ggplot2)
  library(patchwork)
  library(transport)
  library(lpSolve)
  library(vegan)
})

#############
### Setup ###
#############

### Project folder
project_dir <- "landcoverchange"
rda_dir     <- "landcoverchange"
data_dir    <- "input"


### LogFile
#read alaska sites
ala_sites <- readLines(glue::glue("{metadata_path}/ala_sites.txt"))
{
  #get directory names from rda files
  rda_files    <- tibble(Name = list.files(rda_dir)) %>% 
    mutate(PSA = grepl("PSA", Name, fixed = TRUE)) %>%
    filter(PSA) %>% 
    mutate(ID = sapply(strsplit(Name, "_"), 
                       function(x) as.numeric(x[[2]]))) %>%
    dplyr::select(-PSA) %>%
    mutate(lat = as.numeric(sapply(strsplit(Name, "_"), 
                                   function(x) gsub("N", "", x[[length(x)]])))) %>% 
    filter(ID %in% ala_sites)
}

#function do the pixel flow per site
site_pixel_flow <- function(run) {
  ID      <- rda_files$ID[run]
  name    <- rda_files$Name[run]
  dir_out <- glue::glue("{project_dir}{name}")
  
  ### Create folder and copy files
  {
    if (!file.exists(dir_out)) {
      dir.create(dir_out, recursive = TRUE)
      invisible(
        lapply(
          c("psaChange", "summaryResults"),
          function(x) dir.create(glue::glue("{dir_out}/{x}"))
        )
      )
    }
    
    ## Move files
    fls <- list.files(
      glue::glue("{rda_dir}/{rda_files$Name[run]}"),
      recursive = TRUE
    )
    
    invisible(
      lapply(
        c("map_init.png", "map_init_withHuman", "Map.png", "rda_summary.png"),
        function(x) {
          if (any(grepl(x, fls, fixed = TRUE))) {
            file.copy(
              glue::glue("{rda_dir}/{rda_files$Name[run]}/{fls[which(grepl(x, fls, fixed = TRUE))]}"),
              glue::glue("{project_dir}/{rda_files$Name[run]}/{fls[which(grepl(x, fls, fixed = TRUE))]}"),
              overwrite = TRUE
            )
          }
        }
      )
    )
  }
  
  runOut <- suppressWarnings(
    suppressMessages(
      tryCatch({
        
        #### class definitions
        {
          ## Class defs
          class_defs <- readxl::read_excel(
            glue::glue("{data_dir}/settings/landcover_class_definitions_mod.xlsx")
          )
          
          ## Native resolution = 309.2208 meter, possible options: 300, 3000, 5000, 10000
          resolution <- 300
          
          ## Pixel per class for RDA
          pxl_extract_PSA        <- 800
          pxl_extract_PSA_buffer <- 50
          
          ## Maximal spatial buffer around PSA
          max_buffer <- 1000
          
          ## Flow cost parameters
          envs  <- c(0, 100)
          dists <- c(0, 100)
          distF <- 0.05
          
          ## Human classes will be transfered with (TRUE) or without (FALSE)
          ## scaling in relation to landcover distribution in PSA
          scaleDensity <- TRUE
        }
        
        ###################################
        ### Pixel flow                  ###
        ###################################
        
        #### 4.2 Pixel Flow
        if (!file.exists(glue::glue("{dir_out}/psaChange/LcovOut_new.rda"))) {
          
          # loading data created during preprocessing
          load(glue::glue("{rda_dir}/{name}/summaryResults/initRast.rda"))
          load(glue::glue("{rda_dir}/{name}/summaryResults/rdaOut.rda"))
          load(glue::glue("{rda_dir}/{name}/summaryResults/psaLCover_new.rda"))
          load(glue::glue("{rda_dir}/{name}/psaOvlp/psaOvlp.rda"))
          load(glue::glue("{rda_dir}/{name}/psaEnv/env_stars.rda"))
          
          ### Landcover classes
          lcovs     <- st_dimensions(psaOvlp@densities)$attributes$values
          lcovs_num <- as.numeric(gsub("LC_", "", lcovs))
          
          ### Flow table
          psaFlow <- psaLCover@landcover_ts %>% 
            filter(Dataset_ID == ID, Age > 500) %>%
            group_split(Age) %>%
            Reduce("rbind", .)
          
          # point sf with pixel_ID and landcover
          init_crds <- initRast[1] %>% 
            st_as_sf(na.rm = FALSE) %>% 
            st_centroid() %>%
            rowid_to_column(var = "pxlID") %>% 
            relocate("pxlID", .before = "Landcover") %>%
            filter(Landcover %in% as.numeric(names(psaFlow)[-c(1, 2)])) %>% 
            suppressWarnings()
          
          #df with environmental vars for each pixel
          init_env <- st_extract(env_stars, init_crds) %>% 
            st_as_sf() %>% 
            st_drop_geometry() %>%
            apply(., 2, function(x) {
              x[is.na(x)] <- median(x, na.rm = TRUE)
              x
            }) %>% 
            suppressMessages()
          
          # rda coordinates for each pixel
          init_rda <- as.data.frame(
            predict(
              psaOvlp@rda_mod[[1]],
              init_env %>% as.data.frame(),
              type = "wa",
              scaling = 1
            )
          ) %>% 
            dplyr::select("RDA1", "RDA2") %>% 
            st_as_sf(coords = c("RDA1", "RDA2"))
          
          # landcover density at pixel
          init_ovlp <- init_rda %>%
            st_extract(
              lapply(1:dim(psaOvlp@densities)[3], function(x) {
                split(psaOvlp@densities)[x, ] %>% setNames("densities")
              }) %>% 
                do.call("c", .) %>% 
                setNames(lcovs) %>% 
                merge(),
              .
            ) %>% 
            st_as_sf() %>% 
            st_drop_geometry()
          
          # get centers for land cover classes in RDA space
          distCentre <- lapply(1:dim(psaOvlp@densities)[3], function(x) {
            tmp <- split(psaOvlp@densities)[x]
            min <- tmp %>% 
              st_as_sf() %>% 
              pull(1) %>% 
              which.max()
            tmp %>% 
              st_as_sf() %>% 
              slice(min) %>% 
              st_centroid() %>% 
              st_geometry() %>% 
              suppressWarnings()
          })
          
          # distance of each pixel to the centers of landcover classes
          init_dist <- distCentre %>% 
            lapply(function(x) {
              as.numeric(st_distance(init_rda, x))
            }) %>% 
            Reduce("cbind", .) %>% 
            as_tibble() %>% 
            setNames(names(init_ovlp)) %>% 
            suppressWarnings()
          
          ## costs (with scaling)
          if (!is.null(envs)) {
            envs_scale <- max(envs) - scales::rescale(init_ovlp %>% as.matrix(), envs)
          } else {
            envs_scale <- max(init_ovlp %>% as.matrix(), na.rm = TRUE) - init_ovlp %>% as.matrix()
          }
          
          if (!is.null(dists)) {
            dist_scale <- scales::rescale(init_dist %>% as.matrix(), dists)
          } else {
            dist_scale <- init_dist %>% as.matrix()
          }
          
          costs <- abind::abind(envs_scale, dist_scale * distF, along = 3) %>%
            apply(., 1:2, sum, na.rm = TRUE)
          
          classFlow <- matrix(
            nrow = nrow(init_env),
            ncol = length(psaFlow$Age %>% unique()) + 1
          )
          
          classFlow[, 1] <- init_crds$Landcover
          
          #### start FLOW ####
          NA_survey <- data.frame(
            col = numeric(),
            from = numeric(),
            na_count = numeric()
          )
          
          for (y in 2:ncol(classFlow)) {
            print(y)
            
            old_p <- ((tibble(lcov = classFlow[, y - 1]) %>% 
                         group_by(lcov) %>% 
                         summarize(p = n()) %>%
                         right_join(tibble(lcov = lcovs_num), by = "lcov") %>% 
                         mutate(p = ifelse(is.na(p), 0, p)) %>%
                         arrange(lcov) %>% 
                         pull(p)) / nrow(classFlow)) * 100
            
            new_p <- as.numeric(unlist(psaFlow %>% filter(Age == (psaFlow$Age %>% unique())[y-1]))[-c(1,2)])
            
            costM <- 1 - psaOvlp@overlapp[, , 2]
            
            flow <- transport(old_p, 
                              new_p, 
                              costM, 
                              fullreturn = FALSE) %>% 
              as_tibble() %>%
              mutate(from = lcovs_num[from], to = lcovs_num[to]) %>% 
              filter(from != to) %>%
              mutate(
                ident = apply(
                  .,
                  1,
                  function(x) glue::glue("{min(as.numeric(x[1:2]))}_{max(as.numeric(x[1:2]))}")
                )
              ) %>%
              group_split(ident) %>% 
              lapply(function(x) {
                if (nrow(x) > 1) {
                  diff <- x %>% 
                    pull(mass) %>% 
                    diff() %>% 
                    abs()
                  x[which.max(x %>% pull(mass)), ] %>% 
                    mutate(mass = diff)
                } else {
                  x
                }
              }) %>% 
              Reduce("rbind", .) %>% 
              dplyr::select(-ident) %>%
              mutate(nrPxl = floor(nrow(classFlow) * mass / 100)) %>% 
              arrange(from) %>% 
              group_split(from)
            
            flip_pixels <- function(x) {
              cost_sub <- tibble(index = 1:nrow(costs)) %>% 
                bind_cols(costs) %>%     
                filter(classFlow[, y - 1] == unique(x$from)) %>%
                dplyr::select(
                  index,
                  lcovs[lcovs_num %in% (x %>% pull(to))]
                ) %>%
                arrange(apply(.[, -1], 1, sum)) %>% 
                mutate(newClass = NA, .after = index)
              
              trans_remain <- x
              
              for (rep in 1:5) {
                preselect <- apply(
                  cost_sub[, glue::glue("LC_{trans_remain$to}")],
                  1,
                  function(x) which.min(x)
                )
                
                cost_sub   <- cost_sub %>% left_join((trans_remain %>% group_split(to) %>% lapply(function(y) {
                  cost_sub %>% 
                    mutate(newClass_sub = trans_remain$to[preselect]) %>%
                    filter(newClass_sub == y$to) %>% arrange(glue::glue("LC_{y$to}")) %>%
                    slice(1:(min(c(y$nrPxl, nrow(.)))))
                }) %>% Reduce("rbind", .) %>% dplyr::select(index, newClass_sub)), by = 'index') %>%
                  mutate(newClass = ifelse(!is.na(newClass_sub) & is.na(newClass), newClass_sub, newClass)) %>%
                  dplyr::select(-newClass_sub)
                
                trans_remain <- x %>% 
                  group_split(to) %>% 
                  lapply(function(y2) {
                    y2 %>% mutate(nrPxl = nrPxl - sum(cost_sub$newClass == y2$to, na.rm = TRUE))
                  }) %>% 
                  Reduce("rbind", .) %>% 
                  filter(nrPxl > 0)
                
                if (nrow(trans_remain) < 1) break
              }
              
              NAS <- sum(is.na(cost_sub$newClass))
              
              survey <- c(
                col = y,
                from = unique(x$from),
                na_count = NAS
              )
              
              NA_survey[nrow(NA_survey) + 1, ] <<- survey
              
              cost_sub %>% dplyr::select(index, newClass)
            }
            
            index_change <- lapply(flow, flip_pixels) %>% 
              Reduce("rbind", .) %>% 
              filter(!is.na(newClass))
            
            classFlow[, y] <- classFlow[, y - 1]
            classFlow[index_change %>% pull(index), y] <- index_change %>% pull(newClass)
          }
          
          
          #### end FLOW ####
          NA_survey %>% 
            group_by(col) %>% 
            summarize(na_count = sum(na_count)) %>% 
            print(n = 23)
          
          test_out <- lapply(1:ncol(classFlow), function(x) {
            (as_tibble(classFlow)[, x]) %>%
              table(.) %>% 
              as_tibble() %>% 
              mutate(P = (n / nrow(classFlow)) * 100) %>%
              dplyr::select(1, 3) %>% 
              setNames(c("LCov", "P")) %>%
              full_join(tibble(LCov = as.character(lcovs_num)), by = "LCov") %>%
              arrange(LCov) %>% 
              mutate(P = ifelse(is.na(P), 0, P)) %>%
              mutate(Age = as.numeric(c(-1000, psaFlow$Age)[x]), Type = "pred") %>%
              bind_rows(
                psaFlow %>% 
                  filter(Age == as.numeric(c(-1000, psaFlow$Age)[x])) %>%
                  pivot_longer(-c(Age, Dataset_ID)) %>%
                  dplyr::select(name, value, Age) %>%
                  mutate(type = "orig") %>% 
                  setNames(c("LCov", "P", "Age", "Type"))
              )
          }) %>% 
            Reduce("rbind", .)
          
          pl <- ggplot(
            test_out,
            aes(x = P, y = Age, color = Type)
          ) +
            geom_path() +
            facet_wrap(~LCov)
          
          ggsave(
            glue::glue("{dir_out}/summaryResults/flow_comparison_new.png"),
            plot = pl,
            width = 30,
            height = 22,
            units = "cm"
          )
          
          psaPixelFlow <- list(
            crds = init_crds %>% dplyr::select(-Landcover),
            flow = classFlow %>%
              as_tibble() %>%
              setNames(glue::glue("LCov_{c(-1000, psaFlow$Age)}"))
          )
          
          save(
            psaPixelFlow,
            file = glue::glue("{dir_out}/psaChange/psaPixelFlow_new.rda")
          )
        }
        
        #### 4.3 Output Maps
        if (!file.exists(glue::glue("{dir_out}/psaChange/LcovOut_new.tif"))) {
          load(glue::glue("{dir_out}/psaChange/psaPixelFlow_new.rda"))
          load(glue::glue("{rda_dir}/{name}/summaryResults/initRast.rda"))
          
          LcovOut <- lapply(1:ncol(psaPixelFlow$flow), function(x) {
            rast <- initRast[1] %>% setNames(names(psaPixelFlow$flow)[x])
            rast[[1]][psaPixelFlow$crds$pxlID] <- as.numeric(unlist(psaPixelFlow$flow[, x]))
            rast
          }) %>% 
            do.call("c", .) %>% 
            merge() %>% 
            setNames("LandcoverChange")
          
          save(LcovOut, 
               file = glue::glue("{dir_out}/psaChange/LcovOut_new.rda"))
          
          terra::writeRaster(
            as(LcovOut, "SpatRaster"),
            filename = glue::glue("{dir_out}/psaChange/LcovOut_new.tif")
          )
          
          #### Figures
          mergeTab <- class_defs %>% 
            mutate(
              lcov = ifelse(
                is.na(Merge_To_Class),
                Class_Code,
                Merge_To_Class
              )
            ) %>%
            dplyr::select(
              Class_Code, lcov, Class_Plotlabel,
              Predicted_In_Past, Color_Code, Transfer
            )
          
          dim <- as.numeric(dim(LcovOut)[3])
          
          for (i in 1:ceiling(dim / 25)) {
            maps <- ggplot() +
              geom_stars(
                data = LcovOut[, , , which(cut(1:dim, seq(0, 200, by = 25), labels = FALSE) == i)],
                mapping = aes(fill = as.factor(LandcoverChange))
              ) +
              facet_wrap(~attributes, ncol = 5, nrow = 5) +
              scale_fill_manual(
                values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel),
                name = "Landcover"
              ) +
              theme_void() +
              theme(legend.position = "top")
            
            ggsave(
              glue::glue("{project_dir}/{name}/summaryResults/changeMap_{i}_new.png"),
              plot = maps,
              width = 6 * 5,
              height = 4 + 6 * ceiling(length(which(cut(1:dim, seq(0, 200, by = 25), labels = FALSE) == i)) / 5),
              units = "cm",
              limitsize = FALSE,
              bg = "white"
            )
          }
        }
        
      }, error = function(e) e)
    )
  )
  
  if (any(class(runOut) %in% "error")) {
    png(glue::glue("{dir_out}/summaryResults/error_new.png"))
    plot(
      1:10, 1:10,
      type = "n",
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = ""
    )
    text(5, 5, glue::glue("Error Message: \n {runOut$message}"), cex = 2)
    dev.off()
  }
  
  ## end psaPixelFlow.rda
  if (
    file.exists(glue::glue("{dir_out}/psaChange/LcovOut_new.rda")) &
    !any(
      sapply(
        list.files(glue::glue("{dir_out}/summaryResults"), pattern = ".png"),
        function(x) any(c(grepl("error", x, fixed = TRUE), grepl("changeMap", x, fixed = TRUE)))
      )
    )
  ) {
    load(glue::glue("{dir_out}/psaChange/LcovOut_new.rda"))
    class_defs <- readxl::read_excel(
      glue::glue("{data_dir}/settings/landcover_class_definitions_mod.xlsx")
    )
    
    mergeTab <- class_defs %>% 
      mutate(lcov = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
      dplyr::select(
        Class_Code, lcov, Class_Plotlabel,
        Predicted_In_Past, Color_Code, Transfer
      )
    
    dim <- as.numeric(dim(LcovOut)[3])
    
    for (i in 1:ceiling(dim / 25)) {
      maps <- ggplot() +
        geom_stars(
          data = LcovOut[, , , which(cut(1:dim, seq(0, 200, by = 25), labels = FALSE) == i)],
          mapping = aes(fill = as.factor(LandcoverChange))
        ) +
        facet_wrap(~attributes, ncol = 5, nrow = 5) +
        scale_fill_manual(
          values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
          breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
          labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel),
          name = "Landcover"
        ) +
        theme_void() +
        theme(legend.position = "top")
      
      ggsave(
        glue::glue("{project_dir}/{name}/summaryResults/changeMap_{i}_new.png"),
        plot = maps,
        width = 6 * 5,
        height = 4 + 6 * ceiling(length(which(cut(1:dim, seq(0, 200, by = 25), labels = FALSE) == i)) / 5),
        units = "cm",
        limitsize = FALSE,
        bg = "white"
      )
    }
  }
}

### apply to all sites
library(parallel)

cl <- makeCluster(nodes)

clusterEvalQ(cl,
             suppressMessages({
               library(tidyr)
               library(dplyr)
               library(tibble)
               library(readr)
               library(readxl)
               library(sf); sf_use_s2(FALSE)
               library(stars)
               library(ggplot2)
               library(patchwork)
               library(transport)
               library(lpSolve)
               library(vegan)
               library(terra)
             }))

clusterExport(cl,
              ls())

clusterApplyLB(cl,
               1:nrow(rda_files),
               site_pixel_flow)
stopCluster(cl)


