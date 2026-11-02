setClass(
  "psaVeg",
   slots = c(regionMap      = "sfc_MULTIPOLYGON",
             calibrationMap = "sf",
             resolution     = "numeric",
             psas           = "sf",
             interpolation  = "list",
             vegetation     = "data.frame")
)

setClass(
  "psaLCover",
  slots = c(map           = "sfc_MULTIPOLYGON",
            resolution    = "numeric",
            psas          = "sf",
            landcover_ts  = "data.frame")
)

setClass(
  "psaOvlp",
  slots = c(Dataset_ID = "numeric",
            rda_mod = "list",
            transl_lc_rds = "data.frame",
            scores = "data.frame",
            lcov = "character",
            densities = "stars",
            overlapp = "array",
            centers = "data.frame",
            dists = "matrix"))


setMethod("plot", "psaVeg", function(x) {
  
  ps_poly <- psaVeg@psas %>% left_join(psaVeg@vegetation %>% group_by(Dataset_ID) %>% summarise(oldestAge = max(Age)), by = 'Dataset_ID')
  
  ggplot() +
    geom_sf(psaVeg@calibrationMap, mapping = aes(geometry = geometry), fill = "grey80") +
    geom_sf(ps_poly %>% st_buffer(ps_poly$`Pollen_Source_Radius [m]`), mapping = aes(geometry = geometry), alpha = 0.2) +
    geom_sf(psaVeg@regionMap , mapping = aes(geometry = geometry), fill = NA, color = "red") +  #linewidth = 0.6, 
    theme_light()
  
})


# NODES (set in main.R) controls the PSOCK apply calls below.
.psock_lapply <- function(X, FUN, ..., packages = character()) {
  if (!is.numeric(NODES) || length(NODES) != 1L || is.na(NODES) ||
      !is.finite(NODES) || NODES < 1 || NODES != floor(NODES)) {
    stop("NODES must be a positive whole number.")
  }
  if (length(X) == 0L) return(list())

  cl <- parallel::makePSOCKcluster(as.integer(min(length(X), NODES)))
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterCall(cl, function(packages) {
    for (pkg in packages) library(pkg, character.only = TRUE)
    NULL
  }, packages)

  # Workers receive dependencies as arguments, without serializing caller frames.
  environment(FUN) <- .GlobalEnv
  parallel::parLapply(cl, X, FUN, ...)
}

make_psaVeg <- function(ID, psa, radius, buffer, resolution, PSAmax, proj = "laea",
                        veg_folder, 
                        interpolate = TRUE, 
                        fc = 1/500, dt = 500, k = 5, nrCores = NULL, silent = TRUE) {
  
  if(sf_use_s2()) sf_use_s2(F)
  
  psa_sf <- psa %>%
    st_transform(sprintf("+proj=%s +lon_0=%f +lat_0=%f", proj, st_coordinates(.)[1], st_coordinates(.)[2])) %>%
    st_buffer(radius)
  
  
  #### Load vegetation and meta data ###
  ######################################
  if(!silent) cat("\nLoad meta and vegetation data\n")
  
  # read vegetation cover and meta data
  file_names <- list.files(veg_folder)
  
  # nrCores is retained for caller compatibility; NODES controls workers.
  veg_cover <- .psock_lapply(file_names, function(x, veg_folder, silent) {
    dat <- data.table::fread(file.path(veg_folder, x), verbose = !silent,
                            showProgress = !silent, nThread = 1L)
    dplyr::mutate(dat, dataset = x)
  }, veg_folder = veg_folder, silent = silent)
  names(veg_cover) <- file_names
  
  # extract meta data
  meta_cols_veg <- c("Dataset_ID", "Longitude", "Latitude", "Pollen_Source_Radius [m]")
  meta_data     <- unique(do.call(rbind,
                                  lapply(1:length(veg_cover),
                                         function(x){
                                           veg_cover[[x]] %>% dplyr::select(all_of(meta_cols_veg)) %>% 
                                             cbind(File = names(veg_cover)[x])
                                         }))) %>% mutate(Modern = grepl("modern", File))
                  
  # filter PSAs in buffer region
  psas = meta_data %>% sf::st_as_sf(coords = c("Longitude", "Latitude")) %>% sf::st_set_crs(4326) %>% 
    st_transform(st_crs(psa_sf)) %>% st_intersection(psa_sf %>% st_buffer(buffer*1000)) %>% suppressMessages() %>% suppressWarnings()
  psas$`Pollen_Source_Radius [m]`[psas$Dataset_ID==ID] <- min(c(PSAmax*1000, psas$`Pollen_Source_Radius [m]`[psas$Dataset_ID==ID] ))
  
  
  # filter vegetation data based on selected pollen source areas, drop taxa not present, and normalize rows to 1
  {
    
    if(!silent) cat("\nFilter vegetation data within maps")
    
    # extract vegetation cover (Dataset_ID, Age, and taxa percentages)
    veg_cover_filtered = lapply(veg_cover,
                                function(x){
                                  veg_cov_out = x %>% filter(Dataset_ID %in% psas$Dataset_ID)
                                  if(nrow(veg_cov_out)>0) {
                                    veg_cov_out <- veg_cov_out %>% dplyr::select(all_of(c("Dataset_ID",
                                                                                          grep("mean", colnames(x), value = T))))
                                    colnames(veg_cov_out) = gsub(" [mean of cover in %]", "",
                                                                 colnames(veg_cov_out), fixed = T)
                                    colnames(veg_cov_out) = gsub("_mean [yrs BP]", "",
                                                                 colnames(veg_cov_out), fixed = T)
                                    veg_cov_out
                                  } else NULL
                                })
    
    suppressMessages({
      veg_cover_filtered <- veg_cover_filtered[!sapply(veg_cover_filtered, is.null)] %>% 
        Reduce("full_join", .)  %>%
        filter(!is.na(Age)) %>% 
        mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>% as.data.frame() %>%
        left_join(meta_data %>% dplyr::select(Dataset_ID, Modern), by = "Dataset_ID") %>%
        relocate("Modern", .after = "Age")
    }) 
    
    
    
    # reorder columns
    veg_cover_filtered = veg_cover_filtered[,c(meta_cols_veg[meta_cols_veg %in% colnames(veg_cover_filtered)], 
                                               colnames(veg_cover_filtered)[!(colnames(veg_cover_filtered) %in% meta_cols_veg)])] %>%
      dplyr::select(-(which(apply(veg_cover_filtered[,-c(1:3)], 2, sum)==0)+3))
    veg_cover_filtered[,-c(1,2,3)] <- t(apply(veg_cover_filtered[,-c(1,2, 3)],  1, function(u) u/sum(u)))
    
  }
  
  
  ## interpolation
  veg_cover_filtered_interp <- veg_cover_filtered %>% 
    filter(Modern) %>% 
    mutate(Interp = FALSE, .after = Age) %>%
    bind_rows(interpVegTS(veg_cover_filtered %>% 
                            filter(!Modern), 
                          nrCores = nrCores, 
                          fc = fc, 
                          dt = dt, 
                          k = 5, 
                          silent = silent))

  ### Maps ###
  ############
  map <- suppressWarnings(
    rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
      dplyr::select("sovereignt") %>% st_transform(psa_sf %>% st_crs()) %>% st_make_valid()
  )
  region_map <- psa %>% st_transform(st_crs(map)) %>% st_buffer(psas$`Pollen_Source_Radius [m]`[psas$Dataset_ID==ID]) %>% 
    st_geometry() %>% st_cast("MULTIPOLYGON")
  calibration_map <- map %>% st_intersection(psa_sf %>% st_buffer(buffer*1000) %>% st_bbox() %>% st_as_sfc(crs = st_crs(region_map))) %>% suppressWarnings()

  
  new(
    "psaVeg",
    regionMap      = region_map,
    calibrationMap = calibration_map,
    resolution     = resolution,
    psas           = psas,
    interpolation  = list(fc = fc, dt = dt, k = 5),
    vegetation     = veg_cover_filtered_interp
  )
  
}



interpVegTS = function(veg_cover, 
                       nrCores = NULL, 
                       fc = 1/500, 
                       dt = 500, 
                       k = 5, 
                       silent = TRUE){
  
  if(!silent) cat("\n\nInterpolate past vegetation time series")

    suppressMessages({
      # nrCores is retained for caller compatibility; NODES controls workers.
      (.psock_lapply(unique(veg_cover$Dataset_ID),
                          function(x, veg_cover, fc, dt, k, interpolate_fun){
                            
                            veg_site = veg_cover %>% 
                              filter(Dataset_ID == x) %>% 
                              arrange(Age)
                            
                            veg_site %>% 
                              mutate(Interp = FALSE) %>% 
                              relocate(Interp, .after = Age) %>%
                              bind_rows((suppressWarnings({interpolate_fun(df = veg_site,
                                                                                     meta_cols = c("Dataset_ID", "Age"),
                                                                                     fc = fc, 
                                                                                     dt = dt, 
                                                                                     k = k)}) %>% 
                                           mutate(Interp = TRUE) %>% 
                                           relocate(Interp, .after = Age))) %>% 
                              arrange(Age, Interp)
                            
                          }, veg_cover = veg_cover, fc = fc, dt = dt, k = k,
                          interpolate_fun = interpolate_ages_per_site,
                          packages = c("dplyr", "zoo", "corit"))) %>% Reduce("full_join", .)
      })
}

interpolate_ages_per_site = function(df, 
                                     #ID, 
                                     meta_cols,
                                     fc = 1/500, 
                                     dt = 500,
                                     k = 5, 
                                     int.method = "linear",
                                     appliedFilter = "gauss"){
  require(corit)
  
  taxa_cols = colnames(df)[!(colnames(df) %in% meta_cols)] 
  df_taxa = df[, taxa_cols]
  
  # sort out all NA columns if existing
  df_taxa = df_taxa %>% select(which(colSums(is.na(.)) == 0))
  taxa_cols = colnames(df_taxa)[!(colnames(df_taxa) %in% meta_cols)] 
  
  # in case all sample are within one time intervall average the vegetation cover to age dt
  # to keep 1 time slice at least. from interpolation there would results 0 time slices.
  if (df$Age[nrow(df)] <= dt){
    df = data.frame(Dataset_ID = unique(df$Dataset_ID), 
                    Age = dt, 
                    t(apply(df %>% select(-Dataset_ID, -Age), 2, mean)),
                    stringsAsFactors = F)
  } else{
    
    # interpolate to 500 years grid
    df_taxa_interp = do.call(cbind, 
                             lapply(
                               taxa_cols,
                               function(u){
                                 col = zoo(df_taxa[,u], order.by = df$Age)
                                 col_interp = InterpolationMethod(X = col, fc=fc, dt=dt,
                                                                  timser.length = tail(df$Age, n=1), 
                                                                  int.method = int.method, 
                                                                  appliedFilter = appliedFilter, 
                                                                  k=k)
                                 ages = index(col_interp)
                                 df_ret = data.frame(Age = ages, col_interp)
                                 colnames(df_ret) = c("Age", u); row.names(df_ret) = NULL
                                 
                                 return(df_ret)
                                 
                               }))
    df_taxa_interp$Dataset_ID = unique(df$Dataset_ID)
    df = df_taxa_interp[, c("Dataset_ID", "Age", taxa_cols)]
  }
  
  
  return(df)
}

plotVegSite = function(psaVeg, 
                       ID = unique(psaVeg@vegetation$Dataset_ID)[1], 
                       meta_cols = c("Dataset_ID", "Age", "Interp"),
                       n_taxa = NULL, ...) {
  
  df_plot = psaVeg@vegetation %>% filter(Dataset_ID == ID)
  
  # select n_taxa taxa
  subtitle <- ifelse(is.null(n_taxa), "All taxa available", glue::glue("Selection of {n_taxa}")) 
  n_taxa   <- ifelse(is.null(n_taxa), ncol(df_plot)-length(meta_cols), n_taxa)
  df_plot  <- df_plot %>% dplyr::select(c(meta_cols, 
                names(sort(apply((df_plot %>% dplyr::select(-meta_cols)), 2, mean, na.rm = T), decreasing = T))[1:n_taxa])) %>%
                pivot_longer(cols = -c(meta_cols), values_to = "Percentage")
  
  ggplot(df_plot %>% filter(!is.na(Percentage))) +
    geom_path(aes(x = Percentage, y = Age, group = as.factor(Interp), color = as.factor(Interp))) +
    scale_colour_manual(values = c("#E69F00", "grey60"), name = "Interpolated") +
    facet_wrap(~name, ...) +
    labs(title = glue::glue("Pollen source area ID: {ID}"),
         subtitle = subtitle) +
    ylim(min(df_plot$Age[]), max(df_plot$Age)) +
    theme_light() +
    theme(legend.position = "top")

}

modernLC <- function(psaVeg, 
                     psaIDs = NULL, 
                     class_defs, 
                     class_keep, 
                     resolution = 309.2208, 
                     percent = TRUE) {
  
  if(is.null(psaIDs)) {
    psaIDs <- psaVeg@psas$Dataset_ID
  }
  
  if(length(psaIDs)>1000) psaIDs_sub <- split(psaIDs, c(rep(1, (length(psaIDs)/2)+1), rep(2, length(psaIDs) - (length(psaIDs)/2)))) else psaIDs_sub <- list(psaIDs)
      
      tt <- lapply(1:length(psaIDs_sub), function(x) {
        
        roi     = sf_as_ee(psaVeg@psas %>% st_buffer(psaVeg@psas$`Pollen_Source_Radius [m]`) %>% filter(Dataset_ID %in% psaIDs_sub[[x]]) %>% st_transform(4326) %>% st_shift_longitude())
        dataset = ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))
        
        # extract number of pixel per class
        class_areas  = ee$Image$pixelArea()$addBands(dataset)$reduceRegions(
          collection = roi,
          reducer    = ee$Reducer$count()$group(groupField = 1),
          scale      = dataset$projection()$nominalScale()$getInfo()
        )$getInfo()
        
        regLCov = lapply(class_areas$features, function(x) {
          lapply(x[[4]][[5]], function(y) tibble(id   = x[[4]][[1]], 
                                       lcov  = as.factor(as.numeric(y$group)), 
                                       count = y$count)) %>% bind_rows() 
        }) %>% do.call("rbind",.)
        
        ## merge and filter
        out <- left_join(regLCov, 
                         tibble(lcov      = as.factor(class_defs$Class_Code),
                                out_class = as.factor(ifelse(is.na(class_defs$Merge_To_Class),
                                                             class_defs$Class_Code, 
                                                             class_defs$Merge_To_Class))), 
                         by = "lcov") %>%
          group_by(id, out_class) %>% 
          summarise(count = sum(count)) %>% 
          rename(Dataset_ID = id, lcov = out_class) %>%
          group_by(Dataset_ID) %>% 
          mutate(perc = (count/sum(count, na.rm = T))*100) %>% 
          unnest(cols = "Dataset_ID") %>% 
          filter(lcov %in% class_keep$Lcov) %>% 
          suppressMessages()
        
        if(percent) {
        out <- out %>% 
          dplyr::select(-count) %>%
          pivot_wider(names_from = lcov, 
                      values_from = perc, 
                      values_fill = 0) %>% 
          ungroup()
        } else {
          out <- out %>% 
            dplyr::select(-perc) %>%
            pivot_wider(names_from = lcov, 
                        values_from = count, 
                        values_fill = 0) %>% 
            ungroup()
        }
        
        out[,c(1, order(as.numeric(names(out)[-1]), decreasing = F)+1)] %>% 
          pivot_longer(cols = -Dataset_ID)
        
        
    }) %>% Reduce("rbind",.) %>% pivot_wider(names_from = name, values_fill = 0)
}

getPxls <- function(geom, scale = 300) {
  
  bbox <- st_bbox(psaVeg@regionMap)
  x_delta <- (as.numeric(bbox[3]) - as.numeric(bbox[1]))/scale
  y_delta <- (as.numeric(bbox[4]) - as.numeric(bbox[2]))/scale
  
  x_delta*y_delta
  
}

getLCCfromGEE <- function(ee_dataset, sf_roi, y_max = 500, buffer = 0.5) {
  
  pxl <- max(c(1, round(round(getPxls(sf_roi)/250000,0)/2,0)))
  
  psaGrid <- sf_roi %>% st_transform(4326) %>% st_shift_longitude() %>%
    st_make_grid(n=c(pxl,pxl)) %>% st_buffer(buffer) %>% suppressWarnings() %>% suppressMessages()
  
  out_list <- list()
  for(i in 1:length(psaGrid)) {
    out_list[[i]] <- ee_as_stars(ee_dataset$clip(psaGrid[i] %>% sf_as_ee()), via = 'drive', quiet = TRUE, maxPixels = 1e12) %>% 
      st_as_stars() %>% mutate(LandCover_CCI_300.tif = ifelse(LandCover_CCI_300.tif<1, NA, LandCover_CCI_300.tif))
  }
  
  do.call("st_mosaic", out_list)[sf_roi %>% st_transform(st_crs(out_list[[1]])),,,] %>% suppressMessages()
}

getClasses <- function(psaVeg, psaIDs = NULL, radius = 0) {
  
  roi     = sf_as_ee(psaVeg@psas  %>% filter(Dataset_ID %in% psaIDs) %>% st_buffer(radius) %>% st_transform(4326) %>% st_shift_longitude())
  dataset = ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))
  
  # extract number of pixel per class
  class_areas  = ee$Image$pixelArea()$addBands(dataset)$reduceRegions(
    collection = roi,
    reducer    = ee$Reducer$count()$group(groupField = 1),
    scale      = dataset$projection()$nominalScale()$getInfo()
  )$getInfo()
  
  regLCov = lapply(class_areas$features, function(x) {
    lapply(x[[4]][[5]], function(y) tibble(id   = x[[4]][[1]], 
                                           lcov  = as.factor(as.numeric(y$group)), 
                                           count = y$count)) %>% bind_rows() 
  }) %>% do.call("rbind",.)
  
  regLCov %>% filter(lcov %in% class_defs$Class_Code[class_defs$Include_Class_In_Calibration]) %>% pull(lcov) %>%
    as.character() %>% as.numeric()
}

pixelsBuffer <- function(psa, resolution, class_defs, bufferS, cutoff_PSA = 25000, cutoff_buffer = 100) {
  
  if(ee_check(quiet = TRUE) %>% suppressMessages()) {
    
    dataset <- ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))
    buffer  <- lapply(1:length(bufferS), function(x) {
      if(x==1) psa else {
        psa %>% st_buffer(bufferS[x]) %>% st_difference(psa %>% st_buffer(bufferS[x-1])) %>% suppressWarnings()
      }
    })
    
    mergeTab <- class_defs %>% dplyr::select(Class_Code, Merge_To_Class, Predicted_In_Past) %>%
      mutate(Merge_To_Class = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
      filter(Predicted_In_Past, Merge_To_Class<1000) %>% mutate(Class_Code = as.factor(Class_Code))
    
    
    for(r in 1:length(buffer)) {
      
      roi <- buffer[[r]] %>% st_geometry() %>% sf_as_ee()
      
      # extract number of pixel per class
      class_areas  = ee$Image$pixelArea()$addBands(dataset)$reduceRegions(
        collection = roi,
        reducer    = ee$Reducer$count()$group(groupField = 1),
        scale      = resolution
      )$getInfo()
      
      regLCov = lapply(class_areas$features[[1]][[4]]$groups, function(x) {
        tibble(lcov  = as.factor(as.numeric(x$group)), 
               count = x$count)
      }) %>% do.call("rbind",.)
      
      mrg <- regLCov %>%
        right_join(mergeTab, by = join_by("lcov" == "Class_Code")) %>%
        mutate(count = ifelse(is.na(count), 0, count)) %>%
        group_by(Merge_To_Class) %>% summarise(count = sum(count)) %>% rename(lcov = Merge_To_Class) %>% mutate(radius = bufferS[r], .after = lcov)
      
      if(r == 1) {
        out <- mrg
      } else {
        out <- out %>% bind_rows(mrg)
      }
      
      if(all(out %>% filter(radius==0) %>% pull(count)>=cutoff_PSA) | 
         (any(out$radius>0) & all(out %>% filter(radius>0) %>% group_by(lcov) %>% summarise(count = sum(count)) %>% pull(count)>=cutoff_buffer))) break
      
    }
    
    out %>% arrange(lcov, radius)
    
  } else stop("rgee not (correctly) configured. Fix or use GEE = FALSE instead.")
}


regressionModel <- function(psaVeg,
                            lcovShares,
                            vegetation_modern_years = 1000,
                            k_folds = 10,
                            prefilter_taxa = TRUE,
                            pseudocount = 1e-6,
                            alpha_grid = seq(0, 1, by = 0.25),
                            seed_int = 1234,
                            n_workers = 10,
                            inner_nfolds = 3) {
  
  `%>%` <- dplyr::`%>%`
  
  if (!requireNamespace("foreach", quietly = TRUE)) stop("Package 'foreach' is required.")
  if (!requireNamespace("doParallel", quietly = TRUE)) stop("Package 'doParallel' is required.")
  
  normalize_rows <- function(x) {
    x <- as.matrix(x)
    rs <- rowSums(x, na.rm = TRUE)
    if (any(!is.finite(rs)) || any(rs <= 0)) {
      stop("At least one row has non-positive or non-finite sum; cannot normalize rows.")
    }
    sweep(x, 1, rs, "/")
  }
  
  closure <- function(x) {
    x <- as.matrix(x)
    x[!is.finite(x)] <- 0
    x[x < 0] <- 0
    rs <- rowSums(x)
    bad <- !is.finite(rs) | rs <= 0
    if (any(bad)) {
      x[bad, ] <- 1 / ncol(x)
      rs <- rowSums(x)
    }
    sweep(x, 1, rs, "/")
  }
  
  clr_transform <- function(x, eps = 1e-6) {
    x <- as.matrix(x)
    x[x <= 0 | !is.finite(x)] <- eps
    x <- normalize_rows(x)
    gm <- exp(rowMeans(log(x)))
    log(x / gm)
  }
  
  inv_clr_to_comp <- function(z) {
    z <- as.matrix(z)
    closure(exp(z))
  }
  
  aitchison_distance <- function(obs, pred, eps = 1e-6) {
    obs_clr <- clr_transform(obs, eps = eps)
    pred_clr <- clr_transform(pred, eps = eps)
    sqrt(rowSums((obs_clr - pred_clr)^2))
  }
  
  rmse_vec <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
  mae_vec  <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
  r2_vec   <- function(obs, pred) {
    sst <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
    if (!is.finite(sst) || sst == 0) return(NA_real_)
    1 - sum((obs - pred)^2, na.rm = TRUE) / sst
  }
  
  make_folds <- function(n, k, seed = 1234) {
    set.seed(seed)
    sample(rep(seq_len(k), length.out = n))
  }
  
  predict_enet_composition <- function(model_sub, new_data, return_percent = TRUE) {
    taxa_cols <- model_sub$taxa_cols
    lc_cols <- model_sub$lc_cols
    pseudocount <- model_sub$pseudocount
    
    if (!all(taxa_cols %in% names(new_data))) {
      missing_cols <- setdiff(taxa_cols, names(new_data))
      stop("New data is missing required taxa columns: ",
           paste(missing_cols, collapse = ", "))
    }
    
    X_new <- as.data.frame(new_data[, taxa_cols, drop = FALSE])
    for (nm in taxa_cols) X_new[[nm]] <- suppressWarnings(as.numeric(X_new[[nm]]))
    
    X_new <- as.matrix(X_new)
    if (any(!is.finite(X_new))) stop("New data contains NA/NaN/Inf after numeric conversion.")
    if (any(X_new < 0)) stop("New data contains negative predictor values.")
    if (any(rowSums(X_new) <= 0)) stop("At least one new-data row has non-positive taxa sum.")
    
    X_comp <- normalize_rows(X_new + pseudocount)
    X_clr <- clr_transform(X_comp, eps = pseudocount)
    colnames(X_clr) <- taxa_cols
    
    pred_clr <- predict(model_sub$model_full_set, newx = X_clr, s = model_sub$final_lambda)
    pred_clr <- pred_clr[, , 1, drop = FALSE]
    dim(pred_clr) <- c(nrow(X_clr), length(lc_cols))
    
    pred_comp <- inv_clr_to_comp(pred_clr)
    colnames(pred_comp) <- lc_cols
    
    if (return_percent) pred_comp <- pred_comp * 100
    tibble::as_tibble(pred_comp)
  }
  
  # --------------------------------------------------
  # Prepare modern vegetation table(s)
  # --------------------------------------------------
  veg_modern <- psaVeg@vegetation %>%
    dplyr::filter((Interp & Age <= vegetation_modern_years) | Modern) %>%
    dplyr::select(-c(Age, Interp, Modern)) %>%
    dplyr::select(names(.)[apply(., 2, function(x) (sum(is.na(x)) / length(x)) < 0.5)]) %>%
    dplyr::group_by(Dataset_ID) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(apply(., 1, function(x) !all(is.nan(x[-1])))) %>%
    dplyr::mutate(calibrationMap = 1, .after = Dataset_ID) %>%
    dplyr::left_join(lcovShares, by = "Dataset_ID") %>%
    dplyr::group_split(calibrationMap) %>%
    lapply(as.data.frame)
  
  # --------------------------------------------------
  # Light taxa filtering only: remove taxa with no mass
  # --------------------------------------------------
  cols_model <- lapply(veg_modern, function(data) {
    cols_veg <- names(data)[-c(
      which(names(data) %in% c("Dataset_ID", "calibrationMap")),
      which(!is.na(suppressWarnings(as.numeric(names(data)))))
    )]
    cols_lc <- names(data)[which(!is.na(suppressWarnings(as.numeric(names(data)))))]
    
    if (prefilter_taxa) {
      taxa_sums <- sapply(data[, cols_veg, drop = FALSE], function(x) sum(x, na.rm = TRUE))
      cols_veg <- names(taxa_sums)[taxa_sums > 0]
    }
    
    list(cols_veg = cols_veg, cols_lc = cols_lc)
  })
  
  # --------------------------------------------------
  # Parallel backend
  # --------------------------------------------------
  n_workers <- min(n_workers, parallel::detectCores(logical = TRUE))
  cl <- parallel::makeCluster(n_workers)
  doParallel::registerDoParallel(cl)
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  # --------------------------------------------------
  # Fit per calibration map
  # --------------------------------------------------
  model_out <- lapply(seq_along(veg_modern), function(x) {
    
    data_map <- veg_modern[[x]]
    taxa_cols <- cols_model[[x]]$cols_veg
    lc_cols <- cols_model[[x]]$cols_lc
    
    keep <- complete.cases(data_map[, c(taxa_cols, lc_cols), drop = FALSE])
    data_map <- data_map[keep, , drop = FALSE]
    
    if (nrow(data_map) < 3) stop("Too few complete rows for calibration map ", x, ".")
    if (length(taxa_cols) < 2) stop("Too few taxa predictors after filtering for calibration map ", x, ".")
    if (length(lc_cols) < 2) stop("Need at least 2 land-cover classes for calibration map ", x, ".")
    
    X_raw <- as.matrix(data_map[, taxa_cols, drop = FALSE])
    Y_raw <- as.matrix(data_map[, lc_cols, drop = FALSE])
    
    if (any(X_raw < 0, na.rm = TRUE)) stop("Negative pollen values found.")
    if (any(Y_raw < 0, na.rm = TRUE)) stop("Negative land-cover values found.")
    
    X_comp <- normalize_rows(X_raw + pseudocount)
    Y_prop <- Y_raw / 100
    Y_comp <- normalize_rows(Y_prop)
    colnames(Y_comp) <- lc_cols
    
    X_clr <- clr_transform(X_comp, eps = pseudocount)
    colnames(X_clr) <- taxa_cols
    Y_clr <- clr_transform(Y_comp, eps = pseudocount)
    colnames(Y_clr) <- lc_cols
    
    n <- nrow(data_map)
    k_use <- if (is.null(k_folds) || k_folds > n) n else k_folds
    fold_id <- make_folds(n, k_use, seed = seed_int)
    
    # -------- parallel outer CV over folds --------
    #local copies
    alpha_grid_local <- alpha_grid
    X_clr_local <- X_clr
    Y_clr_local <- Y_clr
    fold_id_local <- fold_id
    inner_nfolds_local <- inner_nfolds
    
    fold_results <- foreach::foreach(
      fold = seq_len(k_use),
      .packages = c("glmnet"),
      .export = c(
        "inv_clr_to_comp",
        "alpha_grid_local",
        "X_clr_local",
        "Y_clr_local",
        "fold_id_local",
        "inner_nfolds_local"
      )
    ) %dopar% {
      
      test_idx <- which(fold_id_local == fold)
      train_idx <- setdiff(seq_len(nrow(X_clr_local)), test_idx)
      
      X_train <- X_clr_local[train_idx, , drop = FALSE]
      X_test  <- X_clr_local[test_idx, , drop = FALSE]
      Y_train_clr <- Y_clr_local[train_idx, , drop = FALSE]
      Y_test_clr  <- Y_clr_local[test_idx, , drop = FALSE]
      
      alpha_fold_scores <- rep(NA_real_, length(alpha_grid_local))
      alpha_fold_lambdas <- rep(NA_real_, length(alpha_grid_local))
      alpha_fold_preds <- vector("list", length(alpha_grid_local))
      
      for (i in seq_along(alpha_grid_local)) {
        a <- alpha_grid_local[i]
        
        cvfit <- glmnet::cv.glmnet(
          x = X_train,
          y = Y_train_clr,
          family = "mgaussian",
          alpha = a,
          standardize = TRUE,
          nfolds = min(inner_nfolds_local, length(train_idx)),
          type.measure = "mse",
          parallel = FALSE
        )
        
        lam <- cvfit$lambda.min
        pred_clr <- predict(cvfit, newx = X_test, s = "lambda.min")
        pred_clr <- pred_clr[, , 1, drop = FALSE]
        dim(pred_clr) <- c(length(test_idx), ncol(Y_clr_local))
        
        score <- mean(sqrt(rowSums((Y_test_clr - pred_clr)^2)), na.rm = TRUE)
        
        alpha_fold_scores[i] <- score
        alpha_fold_lambdas[i] <- lam
        alpha_fold_preds[[i]] <- inv_clr_to_comp(pred_clr)
      }
      
      best_i <- which.min(alpha_fold_scores)
      
      list(
        fold = fold,
        test_idx = test_idx,
        best_alpha = alpha_grid_local[best_i],
        best_lambda = alpha_fold_lambdas[best_i],
        best_pred = alpha_fold_preds[[best_i]],
        alpha_scores = alpha_fold_scores
      )
    }
    
    pred_final <- matrix(NA_real_, nrow = n, ncol = ncol(Y_comp), dimnames = list(NULL, lc_cols))
    alpha_by_fold <- rep(NA_real_, k_use)
    lambda_by_fold <- rep(NA_real_, k_use)
    
    alpha_perf <- lapply(alpha_grid, function(a) numeric(0))
    names(alpha_perf) <- as.character(alpha_grid)
    
    for (res in fold_results) {
      pred_final[res$test_idx, ] <- res$best_pred
      alpha_by_fold[res$fold] <- res$best_alpha
      lambda_by_fold[res$fold] <- res$best_lambda
      
      for (i in seq_along(alpha_grid)) {
        alpha_perf[[as.character(alpha_grid[i])]] <- c(
          alpha_perf[[as.character(alpha_grid[i])]],
          res$alpha_scores[i]
        )
      }
    }
    
    mean_Aitchison_by_alpha <- data.frame(
      alpha = alpha_grid,
      mean_Aitchison = vapply(alpha_perf, function(v) mean(v, na.rm = TRUE), numeric(1))
    )
    
    final_alpha <- mean_Aitchison_by_alpha$alpha[which.min(mean_Aitchison_by_alpha$mean_Aitchison)]
    
    final_cvfit <- glmnet::cv.glmnet(
      x = X_clr,
      y = Y_clr,
      family = "mgaussian",
      alpha = final_alpha,
      standardize = TRUE,
      nfolds = min(inner_nfolds, nrow(X_clr)),
      type.measure = "mse",
      parallel = FALSE
    )
    final_lambda <- final_cvfit$lambda.min
    
    coef_list <- glmnet::coef.glmnet(final_cvfit, s = final_lambda)
    intercept <- vapply(seq_along(coef_list), function(j) as.numeric(coef_list[[j]][1, 1]), numeric(1))
    coef_matrix <- do.call(cbind, lapply(seq_along(coef_list), function(j) {
      as.numeric(coef_list[[j]][-1, 1])
    }))
    rownames(coef_matrix) <- taxa_cols
    colnames(coef_matrix) <- lc_cols
    names(intercept) <- lc_cols
    
    per_class_metrics <- dplyr::bind_rows(lapply(seq_len(ncol(Y_comp)), function(j) {
      tibble::tibble(
        class = lc_cols[j],
        RMSE = rmse_vec(Y_comp[, j], pred_final[, j]),
        MAE  = mae_vec(Y_comp[, j], pred_final[, j]),
        R2   = r2_vec(Y_comp[, j], pred_final[, j])
      )
    }))
    
    sample_rmse <- sqrt(rowMeans((Y_comp - pred_final)^2, na.rm = TRUE))
    sample_ait <- aitchison_distance(Y_comp, pred_final, eps = pseudocount)
    
    overall_metrics <- tibble::tibble(
      method = "ElasticNet_CLR",
      mean_RMSE = mean(per_class_metrics$RMSE, na.rm = TRUE),
      mean_MAE = mean(per_class_metrics$MAE, na.rm = TRUE),
      mean_R2 = mean(per_class_metrics$R2, na.rm = TRUE),
      mean_sample_RMSE = mean(sample_rmse, na.rm = TRUE),
      mean_Aitchison = mean(sample_ait, na.rm = TRUE)
    )
    
    samplewise_metrics <- tibble::tibble(
      row_id = seq_len(nrow(Y_comp)),
      sample_RMSE = sample_rmse,
      sample_Aitchison = sample_ait,
      fold = fold_id
    )
    
    model_sub <- list(
      model_type = "ElasticNet_CLR",
      model_full_set = final_cvfit,
      final_alpha = final_alpha,
      final_lambda = final_lambda,
      alpha_grid = alpha_grid,
      alpha_by_fold = alpha_by_fold,
      lambda_by_fold = lambda_by_fold,
      mean_Aitchison_by_alpha = mean_Aitchison_by_alpha,
      coefficients = coef_matrix,
      intercept = intercept,
      taxa_cols = taxa_cols,
      lc_cols = lc_cols,
      pseudocount = pseudocount,
      overall_metrics = overall_metrics,
      per_class_metrics = per_class_metrics,
      samplewise_metrics = samplewise_metrics,
      observed = Y_comp,
      predicted = pred_final,
      calibration_data = data_map[, c("Dataset_ID", "calibrationMap"), drop = FALSE],
      predict_fun = predict_enet_composition
    )
    
    list(
      model_full_set = final_cvfit,
      final_alpha = final_alpha,
      final_lambda = final_lambda,
      alpha_by_fold = alpha_by_fold,
      lambda_by_fold = lambda_by_fold,
      mean_Aitchison_by_alpha = mean_Aitchison_by_alpha,
      coefficients = coef_matrix,
      intercept = intercept,
      mean_metric = overall_metrics$mean_Aitchison,
      sd_metric = stats::sd(sample_ait, na.rm = TRUE),
      mae_per_class = stats::setNames(per_class_metrics$MAE, per_class_metrics$class),
      sd_per_class = stats::setNames(rep(NA_real_, nrow(per_class_metrics)), per_class_metrics$class),
      overall_metrics = overall_metrics,
      per_class_metrics = per_class_metrics,
      samplewise_metrics = samplewise_metrics,
      taxa_cols = taxa_cols,
      lc_cols = lc_cols,
      pseudocount = pseudocount,
      predict_fun = predict_enet_composition,
      model_object = model_sub
    )
  })
  
  names(model_out) <- paste0("calibMap_", seq_along(model_out))
  model_out
}


make_psaLCover <- function(ID, psaVeg, lcovShares, veg_lc_model, model_index = 1) {
  
  
  `%>%` <- dplyr::`%>%`
  
  model_sub <- veg_lc_model[[model_index]]$model_object
  pred_fun <- veg_lc_model[[model_index]]$predict_fun
  
  veg_pred_input <- psaVeg@vegetation %>%
    dplyr::filter(Interp, Dataset_ID == ID) %>%
    dplyr::select(dplyr::all_of(c("Dataset_ID", "Age", model_sub$taxa_cols)))
  
  if (nrow(veg_pred_input) == 0) {
    stop("No interpolated vegetation rows found for Dataset_ID = ", ID)
  }
  
  pred_vals <- pred_fun(
    model_sub,
    veg_pred_input %>% dplyr::select(dplyr::all_of(model_sub$taxa_cols)),
    return_percent = TRUE
  )
  
  pred <- veg_pred_input %>%
    dplyr::select(Dataset_ID, Age) %>%
    dplyr::bind_cols(pred_vals) %>%
    dplyr::filter(rowSums(dplyr::across(dplyr::all_of(model_sub$lc_cols))) != 0)
  
  # add modern share from lcovShares if present
  modern_share <- lcovShares %>%
    dplyr::filter(Dataset_ID %in% pred$Dataset_ID)
  
  out <- pred
  if (nrow(modern_share) > 0) {
    modern_share <- modern_share %>%
      dplyr::mutate(Age = -1000, .after = 1)
    out <- pred %>%
      dplyr::filter(Dataset_ID %in% modern_share$Dataset_ID) %>%
      dplyr::bind_rows(modern_share) %>%
      tibble::as_tibble() %>%
      dplyr::arrange(Dataset_ID, Age)
  } else {
    out <- pred %>%
      tibble::as_tibble() %>%
      dplyr::arrange(Dataset_ID, Age)
  }
  
  psas <- psaVeg@psas %>%
    dplyr::filter(Dataset_ID %in% out$Dataset_ID) %>%
    dplyr::mutate(Region = names(veg_lc_model)[model_index]) %>%
    dplyr::select(Dataset_ID, Region)
  
  new(
    "psaLCover",
    map = psaVeg@regionMap,
    resolution = psaVeg@resolution,
    psas = psas,
    landcover_ts = out
  )
}


plotRDAsummary <- function(ID, dir_out, class_defs) {
  
  require(patchwork)
  
  load(glue::glue("{dir_out}/psaCrds/psaCrds.rda"))
  load(glue::glue("{dir_out}/psaOvlp/psaOvlp.rda"))
  load(glue::glue("{dir_out}/summaryResults/psaVeg.rda"))
  
  colorTab <- class_defs %>% mutate(Merge_To_Class = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>% 
    filter(!is.na(Merge_To_Class), !duplicated(Merge_To_Class), !is.na(Color_Code)) %>% filter(Merge_To_Class%in%as.numeric(gsub("LC_", "", psaOvlp@lcov)))
  
  pl1 <- ggplot() +
    geom_sf(data = psaVeg@calibrationMap$geometry %>% st_intersection(psaVeg@psas %>% filter(Dataset_ID==ID) %>% st_buffer(max(psaCrds$dist)*1000*2))) +
    geom_sf(data = psaVeg@psas %>% filter(Dataset_ID==ID) %>% st_buffer(max(psaCrds$dist)*1000), fill = NA) +
    geom_sf(data = psaVeg@psas %>% filter(Dataset_ID==ID) %>% st_buffer(psaVeg@psas %>% filter(Dataset_ID==ID) %>% pull(2)), mapping = aes(geometry = geometry), fill = NA, color = "red", linewidth = 1, linetype = 2) +
    geom_sf(data = psaVeg@psas %>% filter(Dataset_ID==ID) %>% st_geometry(), mapping = aes(geometry = geometry), fill = NA, color = "red") +  #, linewidth = 1
    geom_sf(data = psaCrds, mapping = aes(geometry = geometry, color = as.factor(lcov)), shape = 16, size = 0.5, show.legend = FALSE) +
    scale_color_manual(values = colorTab$Color_Code, breaks = colorTab$Merge_To_Class, name = "Landcover class", labels = colorTab$Class_Plotlabel) +
    theme_light()
  
  barTab <- psaCrds %>% st_drop_geometry() %>% group_by(lcov) %>% summarise(n = n()) %>% mutate(lcov = factor(lcov, levels = unique(lcov), ordered = F))
  pl2 <- ggplot(barTab, aes(x = lcov, y = n, fill = lcov)) +
    geom_bar(stat="identity", show.legend = FALSE) +
    scale_fill_manual(values = colorTab$Color_Code[match(colorTab$Merge_To_Class, barTab$lcov)]) +
    geom_hline(mapping = aes(yintercept = psaOvlp@scores %>% group_split(lcov) %>% lapply(., nrow) %>% unlist() %>% max()), linetype = 2) +
    xlab("Landcover class") + ylab("Number of pixels") +
    theme_light()
  
  ## Overlap
  ovlTab <- as_tibble(psaOvlp@overlapp[,,1]*100) %>% setNames(psaOvlp@lcov) %>% mutate(lcov = psaOvlp@lcov) %>%
    pivot_longer(-lcov) %>% mutate(lcov = factor(lcov, levels = psaOvlp@lcov), name = factor(name, levels = psaOvlp@lcov))
  
  pl3 <- ggplot(ovlTab) +
    geom_raster(mapping = aes(y = lcov, x = name, fill = value), show.legend = FALSE) +
    scale_fill_viridis_c(na.value = NA, direction = -1, name = "Overlap [%]") +
    xlab("") + ylab("") +
    theme(axis.text = element_text(size = 10))
  
  ## RDA
  colorRDA <- colorTab %>% left_join(psaOvlp@scores %>% group_by(lcov) %>% summarise(sample = n()) %>% 
                                       mutate(Merge_To_Class = as.numeric(gsub("LC_", "", lcov))) %>% dplyr::select(Merge_To_Class, sample), by = join_by(Merge_To_Class))
  
  pl4 <- psaOvlp@scores %>% mutate(lcov_code = gsub("LC_", "", lcov)) %>%
    ggplot(aes(x = RDA1, y = RDA2, color = lcov_code)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = colorRDA$Color_Code, breaks = colorRDA$Merge_To_Class, name = "", 
                       # labels = glue::glue("{colorRDA$Variable_Name} - {colorRDA$Class_Plotlabel} (n = {colorRDA$sample})")) +
                       labels = glue::glue("LC_{colorRDA$Merge_To_Class} - {colorRDA$Class_Plotlabel} (n = {colorRDA$sample})")) +
    theme_light() +
    theme(legend.position = "bottom") +
    guides(color=guide_legend(nrow=4, byrow=TRUE))
  
  
  
  print((pl1 + pl2) / (pl3 + pl4) / guide_area() + plot_layout(guides = "collect"))
  
  
}


k_fold_cv = function(..., data, k, x_names, y_names, analysis_function, 
                     set_seed = F, seed_int = NULL, prnt = F, 
                     tune = F, average="num", verb = 1){
  
  dots = list(...)
  
  if (k > nrow(data)){
    stop("k-fold CV: k must be smaller than the number of rows in data.")
  }
  # build folds
  {
    cv_results = list()
    folds = list() 
    fold.size = nrow(data)/k
    remain = 1:nrow(data)
    for (ii in 1:k){
      
      if(set_seed){
        if (!is.null(seed_int)){
          set.seed(seed_int)
        } else {
          stop("k-fold CV: provide seed_int when set_seed = T")
        }
      }
      
      subsample = sample(remain, fold.size, replace = F)
      
      folds[[ii]] = subsample # store indices
      if (ii == k){
        folds[[ii]] = remain
      }
      
      remain = setdiff(remain, subsample)
    }
    cv_results[["folds"]] = list()
  }
  
  for (ii in 1:k){
    if (prnt){ cat("\n\nFold:", ii) }
    inds = folds[[ii]] 
    train = data[-inds, ] 
    test = data[inds, ]
    
    # analysis is a function that fits, tunes, whatever a model on training set,
    # evaluates it on the test set and computes a metric for the current fold
    # it must return at least the model and the metric
    analysis_results = analysis_function(..., df_train = train, df_test = test, 
                                         x_names = x_names, y_names = y_names,
                                         prnt = prnt, tune = tune, set_seed = set_seed,
                                         verb = verb) 
    if (prnt){ cat("\n\n\tMean absolute error on training set of fold", 
                   ii, ":", analysis_results$metric) }
    
    cv_results$folds[[ii]] = analysis_results
  }
  
  # extract average best parameter over folds
  {
    if (tune){
      best_parameters = do.call(rbind,
                                lapply(cv_results$folds,
                                       function(u){
                                         u$tuning$best$parameter
                                       }))
      if (average == "num"){ 
        tuned_param = mean(best_parameters, na.rm = T) 
      } else if (average == "int"){
        tuned_param = as.integer(round(mean(best_parameters, na.rm = T),0))
      } else {
        stop("Tuning in k-fold CV: Choose average from ('num','int')")
      }
      
      # here hard coded which parameter to change in what way. could be generalized.
      top_taxa = (sort(apply(data[,x_names], 2, sd), decreasing = T))
      keep_taxa = names(top_taxa)
      perc_sd_all = sum(top_taxa)
      x_names = keep_taxa[c(1:tuned_param)]
      perc_sd = sum(top_taxa[1:tuned_param], na.rm = T)
      perc_sd_rel = perc_sd / perc_sd_all * 100
      
      cv_results$tuning = list()
      cv_results$tuning[["x_names"]] = x_names
      cv_results$tuning[["perc_sd_rel"]] = perc_sd_rel
    } 
  }
  
  # extract mean and sd of metrics per fold 
  {
    metrics = data.frame(fold = 1:k,
                         metric = do.call(rbind, lapply(cv_results$folds, 
                                                        function(u){ u$metric })),
                         mean_mae = do.call(rbind, lapply(cv_results$folds, 
                                                          function(u){ u$mae })),
                         do.call(rbind, lapply(cv_results$folds, 
                                               function(u){ u$mae_per_class })))
    mean_metric = mean(metrics$metric, na.rm = T)
    cv_results[["mean_metric"]] = mean_metric
    sd_metric = sd(metrics$metric, na.rm = T)
    cv_results[["sd_metric"]] = sd_metric
    if (prnt){ cat("\n\nMean Metric +/- st. deviation over all ", k, " folds:", 
                   mean_metric, "+/-", sd_metric) }
    
    mean_mae_per_class = apply(metrics[,-c(1:3)], 2, mean, na.rm = T)
    sd_mae_per_class = apply(metrics[,-c(1:3)], 2, sd, na.rm = T)
    cv_results[["mean_mae_per_class"]] = mean_mae_per_class
    cv_results[["sd_mae_per_class"]] = sd_mae_per_class
    
    mean_mae = mean(mean_mae_per_class)
    sd_mae = mean(sd_mae_per_class)
    cv_results[["mean_mae"]] = mean_mae
    cv_results[["sd_mae"]] = sd_mae
  }   
  
  # train and predict final model on full set and pass with mean metric from folds
  analysis_results_full = analysis_function(..., df_train = data, df_test = NULL, 
                                            x_names = x_names, y_names = y_names,
                                            prnt = prnt, tune = F, verb = verb)
  cv_results[["model_full_set"]] = analysis_results_full$model
  cv_results[["results_full_set"]] = analysis_results_full
  cv_results[["x_names"]] = analysis_results_full$x_names
  return(cv_results)
}


multi_reg_analysis = function(..., df_train, df_test, 
                              x_names, y_names, prnt = F, tune = F, verb = 1){
  
  ## initialize ====
  {
    return_list = list()
    
    dots = list(...)
    if ("prnt" %in% names(dots)){ prnt = dots$prnt }
    if ("tune" %in% names(dots)){ tune = dots$tune }
    if ("set_seed" %in% names(dots)){ set_seed = dots$set_seed }
    if ("seed_int" %in% names(dots)){ seed_int = dots$seed_int }
  }
  
  # tune nr. of taxa to include (sorted after highest variation) ====
  
  if (tune){
    if (verb >= 1){ cat("\n\tTuning optimal predictors") }
    return_list[["tuning"]] = list()
    x_names_all = x_names
    
    # do simple train-validation split from df_train
    {
      if (set_seed){
        set.seed(seed_int)
      }
      sub_sample = sample(1:nrow(df_train), ceiling(nrow(df_train)*0.8), 
                          replace = F)
      train_inner = df_train[sub_sample,]
      validation_inner = df_train[-sub_sample,]
    }
    
    for (nr_taxa in 2:max_nr_predictors){
      nr_taxa_name = paste0(nr_taxa, "_taxa")
      return_list$tuning[[nr_taxa_name]] = list()
      x_names = x_names_all[c(1:nr_taxa)]
      
      # how much of mass, variation, and correlation with y do we keep in advance
      {
        keep_sd_n = sum(top_taxa_sd[x_names], na.rm = T)
        keep_sd_rel_n = keep_sd_n / sum(keep_sd, na.rm = T) * 100
        return_list$tuning[[nr_taxa_name]][["perc_sd"]] = keep_sd_rel_n
        
        keep_mean_n = sum(top_taxa_mean[x_names], na.rm = T)
        keep_mean_rel_n = keep_mean_n / sum(keep_mean, na.rm = T) * 100
        return_list$tuning[[nr_taxa_name]][["perc_mean"]] = keep_mean_rel_n
        
        keep_cor_n = sum(top_taxa_cor[x_names,], na.rm = T)
        keep_cor_rel_n = keep_cor_n / sum(keep_cor, na.rm = T) * 100
        return_list$tuning[[nr_taxa_name]][["perc_cor"]] = keep_cor_rel_n
        
        if (verb == 2){
          cat("\n\t\t\tTry ", nr_taxa, " of (total ", length(x_names_all), ") taxa as predictors:",
              "\n\t\t\t\twith relative continental st. deviation of:", round(keep_sd_rel_n, 4),
              "\n\t\t\t\twith relative continental means of:", round(keep_mean_rel_n, 4),
              "\n\t\t\t\trelative continental abs. correlation with landcover:", round(keep_cor_rel_n, 4),
              "\n\t\t\t\t\t% compared to including all predictors.")
        }
      }
      
      fit_results = compute_multi_reg(df = train_inner, 
                                      x_names = x_names, y_names = y_names)
      model = fit_results$model
      return_list$tuning[[nr_taxa_name]][["model"]] = model
      
      # predict inner test data ====
      {
        pred_test = predict(model, newdata = validation_inner[, x_names])  
        df_pred_test = validation_inner; df_pred_test[,y_names] = pred_test
        df_actual_test = validation_inner
        
        # only y columns for metrics evaluation
        predicted = pred_test
        actual = df_actual_test[,y_names]  
        
        # shifting to > 0 % and rescale such that row sums are 100% (part of the model!)
        {
          # add lowest value per row to get rid of negative probabilites
          predicted = t(apply(predicted, 1, 
                              function(u){
                                if (any(u < 0)){
                                  u = u - min(u, na.rm = T)     
                                } else{
                                  u
                                }
                              }))
          row_sums = apply(predicted, 1, sum)
          row_sums[row_sums == 0] = 1
          predicted = predicted / row_sums 
          
          # update df's with shifted and rescaled probabilities
          df_pred_test[,y_names] = predicted
          predicted_save = cbind("Dataset_ID" = validation_inner$Dataset_ID,
                                 as.data.frame(predicted))
        }
        
        return_list$tuning[[nr_taxa_name]][["df_actual"]] = df_actual_test
        return_list$tuning[[nr_taxa_name]][["df_pred"]] = df_pred_test
      }
      
      # compute metrics for inner test data ====
      {
        mae = multi_reg_metrics(predicted = predicted, actual = actual,
                                per_class = F, metric_name = "MAE")
        mae_per_class = multi_reg_metrics(predicted = predicted, actual = actual,
                                          per_class = T, metric_name = "MAE")
        
        return_list$tuning[[nr_taxa_name]][["mae_test"]] = mae
        return_list$tuning[[nr_taxa_name]][["mae_per_class_test"]] = mae_per_class
        if (verb == 2){ cat("\n\t\t\t\tMAE: ", mae) }
      }
      
      # predict inner training data ====
      {
        pred_train = predict(model, newdata = train_inner[, x_names])  
        df_pred_train = train_inner; df_pred_train[,y_names] = pred_train
        df_actual_train = train_inner
        
        # only y columns for metrics evaluation
        predicted = pred_train
        actual = df_actual_train[,y_names]  
        
        # shifting to > 0 % and rescale such that row sums are 100% (part of the model!)
        {
          # add lowest value per row to get rid of negative probabilites
          predicted = t(apply(predicted, 1, 
                              function(u){
                                if (any(u < 0)){
                                  u = u - min(u, na.rm = T)     
                                } else{
                                  u
                                }
                              }))
          row_sums = apply(predicted, 1, sum)
          row_sums[row_sums == 0] = 1
          predicted = predicted / row_sums 
          
          # update df's with shifted and rescaled probabilities
          df_pred_train[,y_names] = predicted
          predicted_save = cbind("Dataset_ID" = train_inner$Dataset_ID,
                                 as.data.frame(predicted))
        }
        
        return_list$tuning[[nr_taxa_name]][["df_actual"]] = df_actual_train
        return_list$tuning[[nr_taxa_name]][["df_pred"]] = df_pred_train
      }
      
      # compute metrics for inner training data ====
      {
        mae = multi_reg_metrics(predicted = predicted, actual = actual,
                                per_class = F, metric_name = "MAE")
        mae_per_class = multi_reg_metrics(predicted = predicted, actual = actual,
                                          per_class = T, metric_name = "MAE")
        
        # additionally all other metrics
        return_list$tuning[[nr_taxa_name]][["mae_train"]] = mae
        return_list$tuning[[nr_taxa_name]][["mae_per_class_train"]] = mae_per_class
      }
    }
    
    # find model / nr predictors with minimum metric
    {
      # extract and plot training vs test metrics
      {
        metrics = data.frame(mae_test = do.call(rbind, 
                                                lapply(return_list$tuning, 
                                                       function(u){ u$mae_test*100})),
                             mae_train = do.call(rbind, 
                                                 lapply(return_list$tuning, 
                                                        function(u){ u$mae_train*100})))
        metrics$nr_taxa = 2:(nrow(metrics)+1) 
        colnames(metrics)[which(colnames(metrics) == "mae_train")] = "Training Set"
        colnames(metrics)[which(colnames(metrics) == "mae_test")] = "Validation Set"
        return_list$tuning[["metrics"]] = metrics
      }
      
      # find optimal model and extract model etc.
      {
        min_metric = min(metrics[["Validation Set"]], na.rm = T)
        min_ind = which.min(metrics[["Validation Set"]])
        metrics[["Validation Set"]][min_ind] == min_metric
        
        return_list$tuning[["best_model"]] = list()
        return_list$tuning$best_model[["model_ind"]] = min_ind
        return_list$tuning$best_model[["min_metric"]] = min_metric
        best_model = return_list$tuning[[min_ind]]$model
        return_list$tuning$best_model[["model"]] = best_model
        
        best_nr_taxa = metrics$nr_taxa[min_ind]
        return_list$tuning$best_model[["parameter"]] = best_nr_taxa
        best_x_names = keep_taxa[1:best_nr_taxa]
        return_list$tuning$best_model[["best_x_names"]] = best_x_names
        
        return_list$tuning$best_model[["perc_mean"]] = return_list$tuning[[min_ind]]$perc_mean
        return_list$tuning$best_model[["perc_sd"]] = return_list$tuning[[min_ind]]$perc_sd
        return_list$tuning$best_model[["perc_cor"]] = return_list$tuning[[min_ind]]$perc_cor
        
        if (verb >= 1){
          cat("\n\n\t\t Minimum MAE for current fold: ", 
              return_list$tuning$best_model$min_metric, "%", 
              "\n\t\t for including:", best_nr_taxa, "taxa", 
              "\n\t\t top variation of predictors with relative st. dev. of:", 
              round(return_list$tuning$best_model$perc_sd,4), " %",
              "\n\t\t relative prob. mass of:", 
              round(return_list$tuning$best_model$perc_mean, 4), " %", 
              "\n\t\t relative abs correlation with y:", 
              round(return_list$tuning$best_model$perc_cor, 4), "%",
              "\n\t\t\t compared to the sums with all predictors (", 
              "passed to the regression function).")
        }
      }
    }
    
    # set x_names now to the optimal predictor set to train on the full training set
    x_names = best_x_names 
  }
  
  # fit on full set ====
  {
    fit_results = compute_multi_reg(..., df = df_train, 
                                    x_names = x_names, y_names = y_names)
    model = fit_results$model
    return_list[["model"]] = model
    return_list[["formula"]] = fit_results$formula
    A = fit_results$A
    return_list[["A"]] = A
    return_list[["lc_classes"]] = fit_results$lc_classes
    return_list[["lc_class_inds"]] = fit_results$lc_class_inds
  }
  
  # predict values ====
  {
    if (!is.null(df_test)){
      
      # if test data is provided, return metric, actual, and predicted for this test data
      pred_test = predict(model, newdata = df_test[, x_names])  
      df_pred_test = df_test; df_pred_test[,y_names] = pred_test
      df_actual_test = df_test
      
      df_pred_train = NULL; df_actual_train = NULL
      
      # only y columns for metrics evaluation
      predicted = pred_test
      actual = df_actual_test[,y_names]  
      
    } else {
      
      # if no test data provided, return metric, actual, and predicted for training data
      pred_train = predict(model, newdata = df_train)
      df_pred_train = df_train; df_pred_train[,y_names] = pred_train
      df_actual_train = df_train
      
      df_pred_test = NULL; df_actual_test = NULL
      
      predicted = pred_train
      actual = df_actual_train[,y_names]
    }
    
    # shifting to > 0 % and rescale such that row sums are 100% (part of the model!)
    {
      return_list[["rescaled"]] = list()
      return_list[["rescaled"]][["max_row_sums_before"]] = max(as.numeric(apply(predicted, 1, sum)))*100
      return_list[["rescaled"]][["min_row_sums_before"]] = min(as.numeric(apply(predicted, 1, sum)))*100
      return_list[["rescaled"]][["max_pred_value_before"]] = max(predicted)*100
      return_list[["rescaled"]][["min_pred_value_before"]] = min(predicted)*100
      
      # add lowest value per row to get rid of negative probabilites
      predicted = t(apply(predicted, 1, 
                          function(u){
                            if (any(u < 0)){
                              u = u - min(u, na.rm = T)     
                            } else{
                              u
                            }
                          }))
      row_sums = apply(predicted, 1, sum)
      row_sums[row_sums == 0] = 1
      predicted = predicted / row_sums 
      
      return_list[["rescaled"]][["max_row_sums_after"]] = max(as.numeric(apply(predicted, 1, sum)))*100
      return_list[["rescaled"]][["min_row_sums_after"]] = min(as.numeric(apply(predicted, 1, sum)))*100
      return_list[["rescaled"]][["max_pred_value_after"]] = max(predicted)*100
      return_list[["rescaled"]][["min_pred_value_after"]] = min(predicted)*100
      
      # update df's with shifted and rescaled probabilities
      if (!is.null(df_test)){
        df_pred_test[,y_names] = predicted
        predicted_save = cbind("Dataset_ID" = df_test$Dataset_ID,
                               as.data.frame(predicted))
      } else {
        df_pred_train[,y_names] = predicted
        predicted_save = cbind("Dataset_ID"=df_train$Dataset_ID,
                               as.data.frame(predicted))
      }
      
      
    }
    
    if (!is.null(df_test)){
      return_list[["df_actual"]] = df_actual_test
      return_list[["df_pred"]] = df_pred_test
    } else {
      return_list[["df_actual"]] = df_actual_train
      return_list[["df_pred"]] = df_pred_train
    }
  }
  
  # compute metrics ====
  {
    mae = multi_reg_metrics(..., predicted = predicted, actual = actual,
                            per_class = F, metric_name = "MAE")
    mae_per_class = multi_reg_metrics(..., predicted = predicted, actual = actual,
                                      per_class = T, metric_name = "MAE")
    
    # over all metric for cross-fold validation
    return_list[["metric"]] = mae 
    
    # additionally other metrics
    return_list[["mae"]] = mae
    return_list[["mae_per_class"]] = mae_per_class
  }
  
  
  return(return_list)
}


multi_reg_metrics = function(..., predicted, actual,
                             per_class = F, metric_name = "MAE"){
  
  dots = list(...)
  if (('per_class' %in% names(dots))){ per_class = dots$per_class }
  if (('metric_name' %in% names(dots))){ metric_name = dots$metric_name }
  
  if(metric_name == "MAE"){
    resids = predicted - actual
    mae_per_class = apply(abs(resids), 2, 
                          mean, na.rm=T)
    mae = mean(mae_per_class, na.rm=T)
    if (per_class){ 
      metric = mae_per_class
    } else {
      metric = mae
    }
    
  } else if (metric_name == "KLD"){
    keras_available = require(keras)
    if (keras_available){
      kld_per_sample = as.array(keras::metric_kullback_leibler_divergence(actual, predicted))
      kld = mean(kld_per_sample, na.rm=T)
      metric = kld            
    } else { stop("Install Keras for metric KLD")}
  } else if (metric_name == "RMSE"){
    
    # default RMSE (root mean squared error) 
    resids = predicted - actual
    resids_sq = resids**2
    rmse_per_class = apply(resids_sq, 2, 
                           function(u){
                             sqrt(mean(u, na.rm = T))    
                           })
    rmse = mean(rmse_per_class, na.rm = T)
    
    if (per_class){ 
      metric = rmse_per_class
    } else {
      metric = rmse
    }
  }
  
  return(metric)
}


compute_multi_reg = function(..., df, x_names, y_names,
                             with_ic = T){
  
  dots = list(...)
  return_list = list()
  
  # seperate data
  {
    x = as.matrix(df[,x_names])
    y = as.matrix(df[,y_names])
    df_model = df[, c(x_names)]
  }
  
  # do regression
  {
    form = as.formula(paste0("y ~ .", c(" + 0", "")[(as.numeric(with_ic)+1)]))   
    reg = lm(form , data = df_model)
    return_list[["model"]] = reg
    return_list[["formula"]] = form
    
    A = reg$coefficients
    lc_classes = colnames(A) 
    lc_class_inds = 1:ncol(A)
    return_list[["A"]] = A
    return_list[["lc_classes"]] = lc_classes
    return_list[["lc_class_inds"]] = lc_class_inds
  }
  
  return(return_list)
}

rdaMod <- function(df, ID, significance_tests = FALSE) {
  
  # df <- dataset_sub_z
  
  ## one-hot-encode land cover type
  lc_type <- df %>% dplyr::select(lcov) %>% rownames_to_column() %>% 
    pivot_wider(names_from = lcov, values_from = lcov) %>%
    reframe(
      across(
        .cols  = c(-rowname),
        .fns   = function(x) ifelse(is.na(x), 0, 1),
        .names = "LC_{col}"
      ))
  
  
  ## ‘lc’: orthogonal linear combinations of the explanatory variable
  ## ‘wa’: more robust to noise in the environmental variables but are a step between constrained towards unconstrained.
  ##     Most of the time, the default is “wa”. Because these are the most robust to noise in the data.
  ## Scaling for species and site scores:
  ##     Either site (1)
  ##     or     species (2)
  ##     scores are scaled by eigenvalues, and the other set of scores is left unscaled,
  ##     or with 3 both are scaled symmetrically by square root of eigenvalues.
  
  dt <- as.data.frame(df %>% dplyr::select(-lcov))
  
  formula = as.formula(paste0("dt ~ ", paste0(names(lc_type), collapse = " + ")))
  rda_model = rda(formula = formula, 
                  data   = as.data.frame(lc_type),
                  scale = FALSE)
  
  ## extract model results
  {
    R2_adj = RsquareAdj(rda_model)$adj.r.squared
    eigen_vals_all = eigenvals(rda_model, "all")
    eigen_vals_constr = eigenvals(rda_model, "constrained")
    eigen_vals_unconstr = eigenvals(rda_model, "unconstrained")
    variance_df = data.frame("Axis" = c(names(eigen_vals_constr), names(eigen_vals_unconstr)), 
                             Variance = c(eigen_vals_constr, eigen_vals_unconstr),
                             Relative_Variance = c(eigen_vals_constr / sum(eigen_vals_all), 
                                                   eigen_vals_unconstr / sum(eigen_vals_all))); row.names(variance_df) = NULL
    
    # see how prediction is done and how the different scalings change euclidic distances in env space
    #vars_transf_rds_wa = as.data.frame(rda_model$CCA$wa) ==  predict(rda_model, new_data = df[, env_cols], type = "wa")
    vars_transf_rds_wa_sc_1 = as.data.frame(scores(rda_model, choices = c(1,2), display = "wa", scaling = 1)) %>% as_tibble() %>%
      bind_cols(tibble(lcov = glue::glue("LC_{df$lcov}"))) %>% relocate(lcov, .before = RDA1)
    #vars_transf_rds_wa_sc_2 = as.data.frame(scores(rda_model, choices = c(1,2), display = "wa", scaling = 2))
    
    # check how scaling types differ in scaling over RDA axes
    #apply(round(vars_transf_rds_wa_sc_1 / vars_transf_rds_wa[,1:2], 6), 2, unique)
    #apply(round(vars_transf_rds_wa_sc_2 / vars_transf_rds_wa[,1:2], 6), 2, unique)
    
    transl_lc_rds = as.data.frame(coef(rda_model))
    row.names(transl_lc_rds) = gsub("df[, y]", "LC", row.names(transl_lc_rds), fixed = T)
    
    # matrix that translates environmental data to RDA axes
    transl_env_rds = as.data.frame(predict(rda_model, new_data = df[, env_cols], type = "sp"))
  }
  
  ## check collinearity and whether to add all predictors
  {
    lc_type_vifs = sqrt(vif.cca(rda_model))
    # cat("\n\t", all(sqrt(vif.cca(rda_model)) < 2, na.rm = T)) # rule of thumb for collinearity too high
  }
  
  ## fitted and residuals on env. variables 
  {
    residuals = as.data.frame(vegan:::residuals.cca(rda_model))
    predicted = as.data.frame(predict(rda_model, new_data = df %>% dplyr::select(-lcov), 
                                      type = "response")) 
    actual = df %>% dplyr::select(-lcov)
    #cat("\n\t", all(round(actual - (residuals + predicted), 6) == 0)) 
    
    # predicted values can also be directly assessed from the model via fitted function
    #cat("\n\t", all(round(fitted(rda_model, type = "response", scaling = 1) - predicted, 6) == 0))
  }
  
  ## significance tests
  if (significance_tests) {
    
    # perform significance tests
    {
      sig_global = anova.cca(rda_model)
      sig_axes = anova.cca(rda_model, by = "axis")
      sig_axes[which(sig_axes$`Pr(>F)` <= p_val),]
      sig_terms = anova.cca(rda_model, by = 'terms', step = 100)
    }
    
    # build data frame with axes' significance
    {
      sig_df = (as.data.frame(sig_terms) %>% select(all_of("Pr(>F)")))
      sig_df$P_LC_Axes = row.names(sig_df)
      colnames(sig_df) = c("P", "Axes"); row.names(sig_df) = NULL
      sig_df = sig_df[-nrow(sig_df), ]; sig_df = sig_df[,c("Axes", "P")]
      
      sig_df = rbind(sig_df, 
                     as.data.frame(sig_axes) %>% 
                       mutate("Axes" = row.names(as.data.frame(sig_axes)),
                              P = !!as.name("Pr(>F)"),
                              .keep = "none"))
      sig_df = sig_df[-nrow(sig_df),]; row.names(sig_df) = NULL
      
      sig_df = rbind(data.frame(Axes = "Model",
                                P = sig_global$`Pr(>F)`[1]), sig_df)
      sig_df = sig_df %>% mutate(Significant = (P <= p_val))
    } 
  } else{
    sig_df = NULL
  }
  
  ## Output
  list(ID = ID,
       rda_model = rda_model, 
       R2_adj = R2_adj, 
       transl_lc_rds = transl_lc_rds, 
       vars_transf_rds_wa_sc_1 = vars_transf_rds_wa_sc_1,
       variance_df = variance_df, 
       significance_df = sig_df)
}

make_psaOvlp <- function(rdaOut, range_axes = NULL) {
  
  if(is.null(range_axes)) range_axes <- apply(rdaOut$vars_transf_rds_wa_sc_1[,2:3], 2, range)
  classes    <- rownames(rdaOut$transl_lc_rds)
  
  ## Density per class
  dens <- lapply(classes, function(cls) {
    scores <- rdaOut$vars_transf_rds_wa_sc_1 %>% filter(lcov == cls)
    if(nrow(scores)<2) {
      sd <- apply((rdaOut$vars_transf_rds_wa_sc_1 %>% group_by(lcov) %>% summarise(sd1 = sd(RDA1), sd2 = sd(RDA2)))[,2:3], 2, median, na.rm = T)
      scores <- tibble(lcov = scores$lcov, RDA1 = rnorm(5, scores$RDA1, sd[1]), RDA2 = rnorm(5, scores$RDA2, sd[2]))
    }
    MASS::kde2d(x = scores %>% pull(RDA1), y =  scores %>% pull(RDA2), 
                h = apply(scores[,2:3], 2, MASS::bandwidth.nrd), 
                n = 100,
                lims = as.numeric(c(range_axes[,1], range_axes[,2])))
  })
  
  x_seq <- dens[[1]]$x
  y_seq <- dens[[2]]$y
  
  z <- merge(c(lapply(dens, function(x) {
    data.frame(x = rep(x_seq, length(y_seq)), 
               y = rep(y_seq, each = length(x_seq)), 
               vals = as.numeric(x[[3]])/sum(x[[3]])) %>%
      st_as_sf(coords = c("x", "y")) %>% 
      st_rasterize()}) %>% 
      Reduce("c",.)) %>% 
      setNames(classes))
  
  ovlp <- apply(expand_grid(a = 1:dim(z)[3], b = 1:dim(z)[3]), 1, function(x) {
    tibble(a = x[1], 
           b = x[2], 
           D = 1 - (0.5 * (sum(abs(z[,,,x[1]][[1]] - z[,,,x[2]][[1]])))),
           I = 1 - (0.5 * (sqrt(sum((sqrt(z[,,,x[1]][[1]]) - sqrt(z[,,,x[2]][[1]]))**2)))))
  }) %>% Reduce("rbind",.)
  
  ovlpMatrix <- abind::abind(matrix(ovlp$D, 
                                    ncol = max(ovlp$a), 
                                    nrow = max(ovlp$b)),
                             matrix(ovlp$I, 
                                    ncol = max(ovlp$a), 
                                    nrow = max(ovlp$b)), along = 3)
  
  ## Centroids
  centers <- lapply(1:dim(z)[3], function(x) {
    stars_expl <- z[,,,x]
    cont <- stars_expl %>% st_contour(contour_lines = FALSE, breaks = c(quantile(stars_expl[[1]], probs = c(0.975, 0.999, 1)))) %>% st_sf() %>% arrange(Min)
    cntr <- suppressWarnings({cont[nrow(cont),] %>% st_centroid() %>% st_coordinates() %>% as_tibble()})
    
    # ggplot() +
    #   geom_stars(data = stars_expl) +
    #   geom_sf(data = cont, mapping = aes(geometry = geometry), fill = NA, col = "white") +
    #   geom_point(data = cntr, mapping = aes(x = X, y = Y), col = "orange")
    # 
   cntr 
    
  }) %>% Reduce("rbind",.) %>% setNames(c("RDA1", "RDA2"))
  
  ## center dists
  center_dists <- as.matrix(dist(centers))
  
  # opar <- par(mfrow = c(2,1))
  # plot(raster::raster(scales::rescale(1- ovlpMatrix[,,1], c(0,1))))
  # plot(raster::raster(scales::rescale(center_dists, c(0,1))))
  # par(opar)
  
  new(
    "psaOvlp",
    Dataset_ID = rdaOut$ID,
    rda_mod = list(rdaOut$rda_model),
    transl_lc_rds = rdaOut$transl_lc_rds,
    scores = rdaOut$vars_transf_rds_wa_sc_1,
    lcov = classes,
    densities = z %>% setNames("densities"),
    overlapp = ovlpMatrix,
    centers = as.data.frame(centers),
    dists = center_dists)
  
  
}


mergeMaps <- function(psaLCover, initRast, psaChangeTab, path = glue::glue("{dir_out}/psaChange"), 
                      method = "dist", selection = "exact") {
  
  emptyRast <- initRast; emptyRast[[1]][] <- NA
  
  psaCenter <- psaLCover@psas %>% st_centroid() %>% suppressWarnings() %>% filter(Dataset_ID %in% psaChangeTab$Dataset_ID)
  rastPts   <- initRast %>% st_as_sf() %>% st_centroid() %>% suppressWarnings() %>% st_geometry()
  
  distRast  <- lapply(1:nrow(psaChangeTab), function(x) {
    outR <- initRast[1]
    outR[[1]] <- as.numeric((st_distance(psaCenter[x,], rastPts %>% st_transform(st_crs(psaCenter))))/1000)
    outR
  }) %>% Reduce("c",.) %>% setNames(psaChangeTab$Dataset_ID) %>% merge()
  
  w <- t(distRast %>% st_as_sf() %>% st_drop_geometry() %>% 
           apply(., 1, function(z) scales::rescale(z, c(0,1), range(z)))) %>% as_tibble() %>% setNames(glue::glue("dist_{1:ncol(.)}"))
  
  
  ages <- unique(psaLCover@landcover_ts %>% filter(Dataset_ID %in% psaChangeTab$Dataset_ID) %>% pull(Age))
  ages <- ages[-which(ages <= 500)]
  ages[1] <- "Init"
  
  mergeRast <- .psock_lapply(ages, function(y, psaChangeTab, path, emptyRast, w) {
    
    selects <- lapply(psaChangeTab$Path, function(p) {
      load(glue::glue("{path}/{p}"))
      
      if(any(st_get_dimension_values(psaChange, 3)%in%glue::glue("year_{y}"))) {
        tibble(val = c(psaChange[,,,glue::glue("year_{y}")][[1]]))
      } else {
        tibble(val = c(emptyRast[[1]][]))
      }
    }) %>% Reduce("bind_cols", .) %>%
      bind_cols(select = w) %>%
      apply(., 1, function(x) {
        vals <- x[grepl("val", names(x), fixed=TRUE)][order(as.numeric(x[grepl("dist", names(x), fixed=TRUE)]))] 
        ifelse(any(!is.na(vals)), vals[!is.na(vals)][1], NA)
      }) %>% suppressMessages()
    
    outRast <- emptyRast; outRast[[1]][] <- selects
    outRast[1] %>% setNames(glue::glue("year_{y}"))
    
  }, psaChangeTab = psaChangeTab, path = path, emptyRast = emptyRast, w = w,
  packages = c("dplyr", "stars")) %>% Reduce("c",.) %>% merge() %>% setNames("LandcoverChange")
  
  mergeRast
}

createFolders = function(folder_path, silent = TRUE, indent = 1) {

    x = list("psaChange" = list(),
             "psaOvlp" = list(),
             "psaEnv" = list(),
             "psaCrds" = list(),
             "summaryResults" = list())


    invisible(lapply(names(x), 
                     function(f) dir.create(glue::glue("{folder_path}{f}"))))

}


