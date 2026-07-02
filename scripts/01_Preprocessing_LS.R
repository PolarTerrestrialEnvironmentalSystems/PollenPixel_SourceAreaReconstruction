#################################
### Template: LandCoverChange ###
### 22.10.2023                ###
### authors: Simeon Lisovski   ###
###         Peter Ewald       ###
###         Thomas Böhmer     ###
#################################

#source functions scripts
source("scripts/LandCoverChange_ElasticNet.R")

#############
### Setup ###
#############

check_files <- TRUE

data_dir    <- "input"
project_dir <- "landcoverchange"

#### Google Earth Engine ####
{
  library(rgee)
  ee_Initialize()

}

### PSAs
metadata_path <- "input/data/psa_metadata" 
psa_metadata  <- read_csv(glue::glue("{metadata_path}/PSA_locations_northern_hemisphere.csv"), 
                          show_col_types = F)


#psa_selection  <- psa_metadata$Dataset_ID[1:5]
psa_selection <- readLines(glue::glue("{metadata_path}/ala_sites.txt")) |> 
  tail(-1) |> 
  as.numeric()

errorTab <- tibble(ID = 1:length(psa_selection), 
                   SA = psa_selection, 
                   time = NA, 
                   timeunit = NA, 
                   status = NA, 
                   errorMsg = NA)
i <- 0
for(pp in psa_selection) {

  rr <- which(psa_metadata$Dataset_ID==pp)
    
  cat("\n")
  print(paste0("PSA: ", 
               psa_metadata$Dataset_ID[rr],
               " (",rr,"/",nrow(psa_metadata),")"))
  
  ID <- psa_metadata$Dataset_ID[rr]

  i <- i + 1
  
  timetrack <- Sys.time()
  failsafe  <- tryCatch({
  
  
  ####################################
  ### Vegetation and Initial Setup ###
  ####################################
  
  #### 1.1 Directory Setup ####
  {
    current_psa  <- psa_metadata[rr,]
    psa_long     <- with(list(Long = current_psa$Longitude), glue::glue("{round(abs(Long),2)}{ifelse(Long>0, 'E', 'W')}"))
    psa_lat      <- with(list(Lat  = current_psa$Latitude),  glue::glue("{round(abs(Lat),2)}{ifelse(Lat>0, 'N', 'S')}"))   ### ifelse(Lat>0, 'S', 'N')

    dataFiles <- tibble(path = list.files(data_dir, recursive = T)) %>%
      mutate(folder = sapply(strsplit(gsub("data/", "", path), "/"), function(x)  x[1]),
             file = sapply(strsplit(gsub("data/", "", path), "/"), function(x)  x[2])) %>%
      dplyr::select(folder, file)

    dir_out <- glue::glue("{project_dir}/PSA_{current_psa$Dataset_ID}__{psa_long}_{psa_lat}/")
    if(!file.exists(dir_out)) {
      dir.create(dir_out)
      createFolders(dir_out)
      fs::dir_tree(path = dir_out, recurse = TRUE)
    }
  }
  
  #### 1.2 Spatial extent and resolution
  {
    ## PSA polygon and output resolution
    psa <- st_point(current_psa %>% dplyr::select(Longitude, Latitude) %>% unlist()) %>% st_sfc(crs = 4326)
    
    ## radius PSA [m]
    radius <- current_psa %>% pull(Pollen_Source_Radius)
    
    ## Trim PSAs if larger PSAmax [km]
    PSAmax <- 250 
    
    ## buffer to include PSAs for calibration [km]
    calibration_buffer  <- 1000            
    
    ## Native resolution = 309.2208 meter, possible options: 300, 3000, 5000, 10000 
    resolution   <- 300    
    
    ## Pixel per class for RDA
    pxl_extract_PSA        <- 800    # max. points for RDA to search within the PSA   
    pxl_extract_PSA_buffer <- 50     # min. points for RDA to search in buffer rings if they are not found in the PSA 
    
    ## Maximal spatial buffer around PSA
    max_buffer             <- 1000        # maximum buffer for search
    
    ## Flow cost parameters
    envs  <- c(0, 100) ## no scaling NULL
    dists <- c(0, 100) ## no scaling NULL
    distF <- 0.05       ## no factor 1
    
    ## Human classes will be transfered with (TRUE) or without (FALSE) scaling in relation to landcover distribution in PSA
    scaleDensity = T
  }
  
  #### 1.3 Create vegetation object for PSA
  {
    # load data, filter pollen source areas in region and/or calibration subset
    if(!file.exists(glue::glue("{dir_out}/summaryResults/psaVeg.rda")) & check_files) {
      psaVeg <- make_psaVeg(ID, psa, radius, calibration_buffer, resolution, PSAmax, proj = "laea",
                            veg_folder = file.path(data_dir, "data", "vegetation_cover"), 
                            fc = 1/500, dt = 500, k = 5, nrCores = 9)
      save(psaVeg, file = glue::glue("{dir_out}/summaryResults/psaVeg.rda"))
      
      plot(psaVeg)
      ggsave(glue::glue("{dir_out}/summaryResults/Map.png"), width = 10, height = 10, units = "cm")
      
    } else load(glue::glue("{dir_out}/summaryResults/psaVeg.rda"))
    
  }
  
  
  ###################################
  ### Landcover                   ###
  ###################################
  
  #### 2.1 Modern LandCover share
  {
    if(!file.exists(glue::glue("{dir_out}/psaEnv/modernVegetation.rda"))) {
      ### modern time threshold for calibration
      modern_age_threshold <- 1000 
      
      modernVegetation <- psaVeg@vegetation %>% 
        filter((Interp & Age <= modern_age_threshold) | Modern) %>%
        dplyr::select(-c(Interp, Age, Modern)) %>%
        group_by(Dataset_ID) %>%
        summarise(across(everything(), function(x) mean(x, na.rm = TRUE))) %>% filter(!apply(., 1, function(x) all(is.nan(x[-1]))))
      save(modernVegetation, file = glue::glue("{dir_out}/psaEnv/modernVegetation.rda"))
    } else load(glue::glue("{dir_out}/psaEnv/modernVegetation.rda"))
    
    class_defs <- readxl::read_excel(glue::glue("{data_dir}/settings/landcover_class_definitions_mod.xlsx")) 
  } ## end 2.1
  
  #### 2.2 CCI LandCover map of PSA
  {
    roi_sf  <- psaVeg@regionMap %>% st_transform(4326)
    roi     <- roi_sf %>% st_shift_longitude() %>% sf_as_ee()
    dataset <- ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))
    
    if(!file.exists(glue::glue("{dir_out}/summaryResults/dataset_rast.tiff")) & check_files) {
      if(diff(range(st_coordinates(roi_sf)[,1]))<100) {
        dataset_rast <- read_stars(glue::glue("{data_dir}/data/global_maps/C3S-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1.tif")) %>%
        st_crop(roi_sf) %>% st_as_stars() %>% suppressWarnings() %>% suppressMessages()
      } else {
        dataset_rast <- getLCCfromGEE(ee_dataset = dataset, sf_roi = psaVeg@regionMap, y_max = 500, buffer = 0.1)
      }
      write_stars(dataset_rast, glue::glue("{dir_out}/summaryResults/dataset_rast.tiff"))
    } else dataset_rast <- read_stars(glue::glue("{dir_out}/summaryResults/dataset_rast.tiff"))

    # extract number of pixel per class
    class_areas  = ee$Image$pixelArea()$addBands(dataset)$reduceRegions(
      collection = roi,
      reducer    = ee$Reducer$count()$group(groupField = 1),
      scale      = dataset$projection()$nominalScale()$getInfo()
    )$getInfo()
    
    regLCov = lapply(class_areas$features[[1]][[4]][[1]], function(x) {
      tibble(Class_Code = x$group, count = x$count)
    }) %>% Reduce("rbind",.) %>% 
      left_join(class_defs %>% 
                  dplyr::select(Class_Code, Variable_Name, Class_Label), 
                by = join_by(Class_Code)) %>%
      mutate(Percent = round((count/sum(count))*100,2)) %>% 
      relocate(Percent, .after = Variable_Name) %>%
      relocate(count, .after = Variable_Name) %>% 
      rename(Count = count) %>% 
      arrange(desc(Count))
    
  } ## end 2.2
  
  #### 2.3 Modern LandCover shares in PSA
  { 
    
    ### Classes present in buffer to use in calibration/prediction
    roi_maxbuffer = sf_as_ee(psaVeg@psas %>% filter(Dataset_ID==ID) %>% st_buffer(.$'Pollen_Source_Radius [m]' + (max_buffer*1000)) %>% 
                             st_transform(4326) %>% st_shift_longitude())    
    
    classes_buffer  = ee$Image$pixelArea()$addBands(dataset)$reduceRegions(
      collection = roi_maxbuffer,
      reducer    = ee$Reducer$count()$frequencyHistogram(),
      scale      = dataset$projection()$nominalScale()$getInfo()
    )$getInfo()
    
    filterLCover = tibble(Class_Code = as.numeric(names(classes_buffer$features[[1]][[4]][[6]])),
                          count = unlist(classes_buffer$features[[1]][[4]][[6]])) %>%
      left_join(class_defs %>% 
                  mutate(Merge_To_Class = ifelse(!is.na(Merge_To_Class), Merge_To_Class, Class_Code)) %>%
                  filter(Include_Class_In_Calibration) %>% 
                  dplyr::select(Class_Code, Merge_To_Class), 
                by = 'Class_Code') %>%
      filter(!is.na(Merge_To_Class)) %>% 
      group_by(Merge_To_Class) %>% 
      summarise(count = sum(count)) %>% 
      filter(count>=pxl_extract_PSA_buffer) %>%
      rename(Lcov = Merge_To_Class)
    
    if(!file.exists(glue::glue("{dir_out}/summaryResults/lcovShares.rda"))) {
      
      lcovShares <- modernLC(psaVeg, 
                             psaIDs = unique(modernVegetation$Dataset_ID), 
                             class_defs = class_defs, 
                             class_keep = filterLCover, 
                             resolution = resolution)
      
      save(lcovShares, file = glue::glue("{dir_out}/summaryResults/lcovShares.rda"))
    } else load(glue::glue("{dir_out}/summaryResults/lcovShares.rda"))
    
    
  } ## end 2.3
  
  #### 2.4 Translate Vegetation to Landcover
  {
    if(!file.exists(glue::glue("{dir_out}/summaryResults/veg_lc_model_new.rda")) & check_files) {
      Sys.unsetenv("LD_LIBRARY_PATH")
      cl <- parallel::makeCluster(20)
      doParallel::registerDoParallel(cl)
      veg_lc_model <- regressionModel(psaVeg, 
                                      lcovShares, 
                                      vegetation_modern_years = 1000, 
                                      prefilter_taxa = T) %>% 
        suppressWarnings()
      stopCluster(cl)
      
      save(veg_lc_model, 
           file = glue::glue("{dir_out}/summaryResults/veg_lc_model_new.rda"))
      
    } else load(glue::glue("{dir_out}/summaryResults/veg_lc_model_new.rda"))
    
    if(!file.exists(glue::glue("{dir_out}/summaryResults/psaLCover_new.rda")) & check_files) {
      psaLCover <- make_psaLCover(ID = ID, psaVeg, lcovShares, veg_lc_model)
      save(psaLCover, file = glue::glue("{dir_out}/summaryResults/psaLCover_new.rda"))
    } else load(glue::glue("{dir_out}/summaryResults/psaLCover_new.rda"))
    
  } ## end 2.4
  
  
  ###################################
  ### RDAs                        ###
  ###################################
  
  #### 3.1 Pixel selection
  {
    if(!file.exists(glue::glue("{dir_out}/psaCrds/psaCrds.rda")) & check_files) {
      
      ## filterLCover
      psa      <- psaVeg@psas %>% 
        filter(Dataset_ID==ID) %>% 
        st_buffer(.$'Pollen_Source_Radius [m]')
      bufferS  <- seq(0, max_buffer*1000, length = 10)
      
      nPxl <- pixelsBuffer(psa, 
                           psaLCover@resolution, 
                           class_defs, 
                           bufferS, 
                           cutoff_PSA = pxl_extract_PSA, 
                           cutoff_buffer = pxl_extract_PSA_buffer) %>% 
        filter(lcov%in%filterLCover$Lcov) %>% 
        group_split(lcov) %>% 
        lapply(., function(x) {
        if(x %>% filter(radius==0) %>% pull(count) >  pxl_extract_PSA) {
          x %>% filter(radius==0) %>% mutate(cumsum = count)
        } else {
          tmp <- x %>% mutate(cumsum = cumsum(count),
                              over   = c(FALSE, (cumsum>pxl_extract_PSA_buffer)[-1]))
          if(sum(tmp$over)>1) {
            tmp[1:min(which(tmp$over)),] %>% 
              filter(count>0) %>% 
              dplyr::select(-over)
            } else {
              tmp %>% 
                filter(count>0) %>% 
                dplyr::select(-over) 
            }
        }
      }) %>% Reduce("rbind",.) 
                      
      psaCrds <- lapply(unique(nPxl$lcov), function(class) {
        
        classTab <- nPxl %>% filter(lcov==class)
        
        lapply(1:nrow(classTab), function(rad) {
          
          if(classTab$radius[rad]==0) {
            bbox_ee  <- psa %>% st_geometry() %>% sf_as_ee()
          } else {
            bbox_ee <- psa %>% st_buffer(classTab$radius[rad]) %>% 
              st_difference(psa %>% st_buffer(bufferS[which(bufferS==classTab$radius[rad])-1])) %>%
              st_geometry() %>% sf_as_ee() %>% suppressWarnings()
          }
          
          dataset  <- ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))$clip(bbox_ee)
            
          fromList = class_defs %>% mutate(Merge_To_Class = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
              filter(Predicted_In_Past, Merge_To_Class == class) %>% pull(Class_Code) %>% as.list()
            
          toList = rep(1, length(fromList)) %>% as.list()
            
          classImage = dataset$remap(
            from = fromList,
            to = toList,
            defaultValue = -1,
            bandName = 'b1')
            
          dts <- dataset$updateMask(classImage$gt(0))
 
          table = dts$stratifiedSample(
            numPoints = min(c(classTab$count[rad], pxl_extract_PSA)),
            region = bbox_ee,
            scale = (dts$projection()$nominalScale()$getInfo()),
            geometries = TRUE
          )$getInfo()
            
          lapply(table[[4]], function(f) {
                tibble(Dataset_ID = current_psa$Dataset_ID,
                         lcov       = as.numeric(as.character(class)),
                         dist       = NA,
                         lon        = f$geometry$coordinates[1],
                         lat        = f$geometry$coordinates[2])
                }) %>% Reduce("rbind",.) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
            mutate(dist = (st_distance(., psa %>% st_centroid() %>% st_geometry() %>% st_transform(4326))/1000)[,1]) %>%
            suppressWarnings()
          
        }) %>% Reduce("rbind", .)
        
      }) %>% Reduce("rbind",.)  
      
      psaCrds_summary <- psaCrds %>% st_drop_geometry() %>% group_by(lcov) %>% summarise(pixel_found = n())
      
      # save pixel summary
      openxlsx::write.xlsx(psaCrds_summary, file = glue::glue("{dir_out}/summaryResults/psaCrds_summary.xlsx"), asTable = TRUE)
  
      ### save results
      save(psaCrds, file = glue::glue("{dir_out}/psaCrds/psaCrds.rda"))
        
    } 
  
  }## end 3.1
  
  #### 3.2 Add climate variables & z-transformation parameters
  {
    threshold <- 0.15 ### 15 percentile of sample size across samples in one PSA
      
    load(glue::glue("{dir_out}/psaCrds/psaCrds.rda"))
      
      if(!file.exists(glue::glue("{dir_out}/psaEnv/psaEnv.rda")) & check_files) {
        
        ## Climate maps
        clim_variables <- read_xlsx(glue::glue("{data_dir}/settings/bioclim_vars_definitions.xlsx"))
        clim_dir       <- list.files(glue::glue("{data_dir}/data/global_maps"), pattern = "*.tif", full.names = T)
        
        clim_stars <- lapply(clim_variables$Var_Name, function(x) {
          ind <- which(sapply(clim_dir, function(y) {grepl(glue::glue("_{x}_"), y, fixed = T)}))
          read_stars(clim_dir[ind]) %>% 
            setNames(x) %>% 
            st_crop(psaVeg@regionMap %>% 
                      st_buffer(500) %>% 
                      st_transform(4326)) %>% 
            st_as_stars() %>% 
            suppressMessages()
        })
        
        elev_stars <- read_stars(clim_dir[grepl("ETOPO", clim_dir)]) %>% st_set_crs(4326)  %>%
          st_warp(., clim_stars[[1]], use_gdal = F, method = "near") %>%
          setNames("elev") %>% st_crop(psaVeg@regionMap %>% st_buffer(500) %>% st_transform(4326)) %>% suppressMessages() %>% suppressWarnings()
        
        ### Environment for selected crds
        psaEnv_init <-  lapply(clim_variables$Var_Name, function(x)
                          {
                          read_stars(clim_dir[which(sapply(clim_dir, function(y) grepl(glue::glue("_{x}_"), y, fixed = T)))]) %>% 
                            setNames(x) %>% suppressMessages() %>% st_extract(psaCrds) %>% 
                            pull(1) %>% unlist()}) %>% Reduce("cbind", .) %>% as_tibble() %>% setNames(clim_variables$Var_Name) %>% 
                          mutate(elev = st_extract(
                            read_stars(clim_dir[grepl("ETOPO", clim_dir)]) %>% st_set_crs(4326) %>% setNames("elev"),
                            psaCrds) %>% pull(1))
        
        ### Environment for selected crds and complete psa
        allEnv_init <- lapply(1:ncol(psaEnv_init), function(x) {
          if(x<ncol(psaEnv_init)) {
            tibble(var = c(clim_stars[[x]][[1]], psaEnv_init %>% pull(x))) %>% filter(!is.na(var)) %>% pull(1)
          } else tibble(var = c(elev_stars[[1]], psaEnv_init %>% pull(x))) %>% filter(!is.na(var)) %>% pull(1)
        })
          
        ## z transformation table
        z_tranTab <- lapply(1:ncol(psaEnv_init), function(x) {
          tibble(Var_Name = names(psaEnv_init)[x], 
                 mean = mean(allEnv_init[[x]], na.rm = T),
                 sd   = sd(allEnv_init[[x]], na.rm = T))
        }) %>% Reduce("rbind", .)
        
        ## apply z transformation to psaEnv
        psaEnv <- lapply(1:ncol(psaEnv_init), function(x) {
          dat = psaEnv_init %>% pull(x)
          dat[is.na(dat)] <- median(dat, na.rm = T)
          tibble(c(scale(dat, z_tranTab[x,2], z_tranTab[x,3])))
        }) %>% Reduce("cbind",.) %>% setNames(names(psaEnv_init))
        
        ## apply z transformation to environment
        env_stars <- lapply(1:nrow(z_tranTab), function(x) {
          append(clim_stars, list(elev_stars))[[x]] %>% setNames("var") %>%
            mutate(var = scale(var, z_tranTab[x,2], z_tranTab[x,3]),
                   var = ifelse(is.na(var), 0, var))
        }) %>% do.call("c", .) %>% setNames(z_tranTab$Var_Name) %>% merge()
  
        save(env_stars, file = glue::glue("{dir_out}/psaEnv/env_stars.rda"))
        save(psaEnv,    file = glue::glue("{dir_out}/psaEnv/psaEnv.rda"))
        save(z_tranTab, file = glue::glue("{dir_out}/psaEnv/z_tranTab.rda"))
      }
    } ## end 3.2
    
  #### 3.3 RDA models   
  {
  
    ## env_stars: z-transformed environmental variables
    ## psaCrds:   coordinates of selected pixels per class
    ## psaEnv:    environmental variable (z-transformed)
    
    if(!file.exists(glue::glue("{dir_out}/psaOvlp/psaOvlp.rda")) & check_files) {
        
        load(glue::glue("{dir_out}/psaEnv/env_stars.rda")) ## env_stars
        load(glue::glue("{dir_out}/psaCrds/psaCrds.rda"))  ## psaCrds
        load(glue::glue("{dir_out}/psaEnv/psaEnv.rda"))    ## psaEnv
        
        datst <- psaCrds %>% st_drop_geometry() %>% dplyr::select(-Dataset_ID) %>%
          bind_cols(psaEnv) %>% mutate(lcov = as.factor(lcov))
        
        sampleSize  <- datst %>% group_by(lcov) %>% summarise(sample = n())
        thrsh       <- quantile(sampleSize$sample, probs = threshold)
        
        dataset_sub <- sampleSize %>% mutate(sample = ifelse(sample>thrsh, thrsh, sample)) %>% 
          group_split(lcov) %>% lapply(., function(x) {
            datst %>% filter(lcov == x$lcov) %>% arrange(dist) %>% dplyr::select(-dist) %>% slice(1:x$sample)
          }) %>% Reduce("rbind", .) %>% mutate(source = "calibration")
      
        rdaOut <- rdaMod(dataset_sub %>% dplyr::select(-source), ID = ID)
        
        save(rdaOut, file = glue::glue("{dir_out}/summaryResults/rdaOut.rda"))
        
        ### Overlap and centers
        #### ovlp range
        zTransEnvRDA <- as.data.frame(
            predict(
              rdaOut$rda_model,
              env_stars %>% st_as_sf() %>% st_drop_geometry() %>% as.data.frame(),
              type = "wa",
              scaling = 1
            )) %>% dplyr::select("RDA1", "RDA2")
        
        psaOvlp <- make_psaOvlp(rdaOut, range_axes = apply(zTransEnvRDA, 2, range))
        
        save(psaOvlp, file = glue::glue("{dir_out}/psaOvlp/psaOvlp.rda"))
      }
  } ## end 3.3
    
  ### 3.4 plot RDA output
  {
    if(!file.exists(glue::glue("{dir_out}/summaryResults/rda_summary.png")) & check_files) {
      load(glue::glue("{dir_out}/psaOvlp/psaOvlp.rda"))
      plotRDAsummary(ID, dir_out, class_defs)
      ggsave(glue::glue("{dir_out}/summaryResults/rda_summary.png"), width = 30, height = 22, units = "cm")
    }
  }
  
  ###################################
  ### Pixel flow                  ###
  ###################################
  
  #### 4.1 Init Maps
  {

   ## load(z-trans variables)
   load(glue::glue("{dir_out}/psaEnv/env_stars.rda"))
   
   ## LandCover Init Map 
   if(!file.exists(glue::glue("{dir_out}/summaryResults/map_init.png"))  & check_files) {
      dataset_rast <- read_stars(glue::glue("{dir_out}/summaryResults/dataset_rast.tiff"))
        
      mergeTab    <- class_defs %>% 
        mutate(lcov = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
          dplyr::select(Class_Code, lcov, Class_Plotlabel, Predicted_In_Past, Color_Code, Transfer)
        
      initRast    <- dataset_rast %>% 
        st_as_stars() %>% 
        setNames("Landcover") %>%
        mutate(Landcover = mergeTab$lcov[match(Landcover, mergeTab$Class_Code)]) %>%
        mutate(predict = mergeTab$Predicted_In_Past[match(Landcover, mergeTab$lcov)])
    
      
      ### Class transfer
      {
        trans <- class_defs %>% 
          mutate(Merge_To_Class = ifelse(!is.na(Merge_To_Class), 
                                         Merge_To_Class, 
                                         Class_Code)) %>% 
          filter(Transfer) %>% 
          pull(Merge_To_Class) %>% 
          unique()
        
        if(any(unique(c(initRast[[1]]))%in%trans)){
          
          pl_human <- ggplot() +
            geom_stars(data = initRast[1], mapping = aes(fill = as.factor(Landcover))) +
            scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                              breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                              labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
            xlab("") + ylab("") +
            theme_light()
          
          ggsave(plot = pl_human, file = glue::glue("{dir_out}/summaryResults/map_init_withHuman.png"), width = 20, height = 12, units = "cm")
          
          transPxl <- initRast[1] %>% 
            st_as_sf(na.rm = FALSE) %>% 
            st_centroid() %>% 
            mutate(index = 1:nrow(.))%>% 
            suppressWarnings()
          transEnv <- st_extract(env_stars, 
                                 transPxl %>% 
                                   filter(Landcover %in% trans))%>% 
            st_as_sf() %>% 
            st_drop_geometry()  
          
          if(scaleDensity) {
            ovl_scaled <- lapply(1:dim(psaOvlp@densities)[3], function(x) {
              (split(psaOvlp@densities)[x,] %>% 
                 setNames("densities"))*as.numeric((lcovShares %>% filter(Dataset_ID==ID))[x+1])
            }) %>% 
              do.call("c", .) %>% 
              setNames(names(split(psaOvlp@densities))) %>% 
              merge()
          } else ovl_scaled <- psaOvlp@densities
            
          transRDA <- st_extract(
            ovl_scaled,
              predict(
                psaOvlp@rda_mod[[1]],
                newdata = transEnv %>% as.data.frame(),
                type = "wa",
                scaling = 1
              )[,1:2] %>% 
              as.data.frame() %>% 
              st_as_sf(coords = c("RDA1", "RDA2"))
          ) %>% 
            st_as_sf()
          
          toClass <- as.numeric(gsub("LC_", 
                                     "", 
                                     names(transRDA %>% st_drop_geometry())))
          
          transDists <- (lapply(1:nrow(psaOvlp@centers), function(x) 
            (st_point(as.numeric(psaOvlp@centers[x,])) %>% 
               st_sfc() %>% 
               st_distance(., transRDA %>% 
                             st_geometry()))[1,]) %>%
            Reduce("cbind", .) %>% 
              as.matrix() %>% 
              scales::rescale(., dists))*distF
          
          transCls <- toClass[abind::abind(envs[2]-(transRDA %>% 
                                                      st_drop_geometry() %>% 
                                                      as.matrix() %>% 
                                                      scales::rescale(., envs)),
                                           transDists, along = 3) %>% 
                                apply(., 1:2, function(x) sum(x)) %>% 
                                apply(., 1, function(x) which.min(x))]
          
          transPxl$Landcover[as.numeric(transPxl %>% filter(Landcover %in% trans) %>% pull(index))] <- transCls
          initRast[1] <- initRast[1] %>% 
            mutate(Landcover = transPxl %>% 
                     pull(Landcover))
        }
      }
          
      ggplot() +
        geom_stars(data = initRast[1], mapping = aes(fill = as.factor(Landcover))) +
        # geom_sf(psaVeg@psas %>% filter(Dataset_ID %in% ID) %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="cyan") + 
        scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                          breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                          labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
        xlab("") + ylab("") +
        theme_light()
      
      save(initRast, file = glue::glue("{dir_out}/summaryResults/initRast.rda"))
      ggsave(glue::glue("{dir_out}/summaryResults/map_init.png"), 
             width = 20, 
             height = 12, 
             units = "cm")
   }
  
  }
    
  }, error = function(e) e)
  
  # -------------------------------------------------------------------------------------------------
  
  if(any(class(failsafe)=="error")){ 
    
    errorTab[which(errorTab$ID == rr),5] <- "failed"
    errorTab[which(errorTab$ID == rr),6] <- failsafe[1]
    
  }else{errorTab[which(errorTab$ID == rr),5] <- "ok"}
  
  errorTab[which(errorTab$ID == rr),3] <- round(as.numeric(difftime(Sys.time(), timetrack, units = "mins")), digits =2) 
  errorTab[which(errorTab$ID == rr),4] <- "min"
  
  # -------------------------------------------------------------------------------------------------
  
  ### save overview tables every 50 sites:
  
  if(i %% 2 == 0 | i == length(psa_selection)) {
    
    write.table(errorTab, 
                file = glue::glue("{project_dir}2024_02_13_logfile_until_pixelflow_part1.csv"), 
                sep="\t", 
                row.names=FALSE)
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
}