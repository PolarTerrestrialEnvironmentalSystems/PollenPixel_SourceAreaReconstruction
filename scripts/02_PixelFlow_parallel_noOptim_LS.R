#################################
### PollenPixel: pixel flow    ###
#################################
rm(list = setdiff(ls(), c("NODES", "RETRY_FAILED_FLOW")))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(readxl)
  library(sf)
  library(stars)
  library(ggplot2)
  library(vegan)
  library(terra)
})
sf_use_s2(FALSE)

source("scripts/LandCoverChange_ElasticNet.R")
source("scripts/ice_cover.R")
invisible(.worker_count())
configure_spatial_runtime()

project_dir <- "landcoverchange"
rda_dir <- "landcoverchange"
data_dir <- "input"
ice_settings <- ice_cover_settings(data_dir)
ice_source <- prepare_ice_cover(ice_settings)

ala_sites <- readLines("input/data/psa_metadata/ala_sites.txt") |>
  tail(-1) |>
  as.numeric()

rda_files <- tibble(Name = list.files(rda_dir)) %>%
  filter(grepl("^PSA_", Name)) %>%
  mutate(
    ID = as.numeric(vapply(strsplit(Name, "_"), function(parts) parts[2], ""))
  ) %>%
  filter(ID %in% ala_sites)

flow_status_file <- file.path(project_dir, "ice_flow_status.csv")
previous_flow_status <- NULL
if (isTRUE(get0("RETRY_FAILED_FLOW", inherits = FALSE))) {
  rm(RETRY_FAILED_FLOW)
  previous_flow_status <- read.csv(flow_status_file, stringsAsFactors = FALSE)
  if (!all(c("ID", "status", "message") %in% names(previous_flow_status)) ||
      anyNA(previous_flow_status$ID) || anyDuplicated(previous_flow_status$ID)) {
    stop("Expected one pixel-flow status row per site in ", flow_status_file)
  }
  failed_ids <- previous_flow_status$ID[which(previous_flow_status$status == "failed")]
  rda_files <- filter(rda_files, ID %in% failed_ids)
  if (!nrow(rda_files)) stop("No failed sites available to retry in ", flow_status_file)
  message("Retrying pixel flow for ", nrow(rda_files), " sites from ", flow_status_file)
}

plot_landcover_flow <- function(LcovOut, class_defs, dir_out) {
  legend <- class_defs %>%
    mutate(lcov = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
    filter(!is.na(Color_Code), !duplicated(lcov))

  layer_count <- dim(LcovOut)[3]
  for (page in seq_len(ceiling(layer_count / 25))) {
    layers <- seq.int((page - 1) * 25 + 1, min(page * 25, layer_count))
    maps <- ggplot() +
      geom_stars(
        data = LcovOut[, , , layers],
        mapping = aes(fill = as.factor(LandcoverChange))
      ) +
      facet_wrap(~attributes, ncol = 5) +
      scale_fill_manual(
        values = setNames(legend$Color_Code, legend$lcov),
        breaks = legend$lcov,
        labels = legend$Class_Plotlabel,
        name = "Landcover"
      ) +
      theme_void() +
      theme(legend.position = "top")

    ggsave(
      file.path(dir_out, "summaryResults", paste0("changeMap_", page, "_new.png")),
      plot = maps,
      width = 30,
      height = 4 + 6 * ceiling(length(layers) / 5),
      units = "cm",
      limitsize = FALSE,
      bg = "white"
    )
  }
}

site_pixel_flow <- function(run) {
  ID <- rda_files$ID[run]
  name <- rda_files$Name[run]
  dir_out <- file.path(project_dir, name)
  input_dir <- file.path(rda_dir, name)
  summary_dir <- file.path(dir_out, "summaryResults")
  flow_dir <- file.path(dir_out, "psaChange")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(flow_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    class_defs <- readxl::read_excel(
      file.path(data_dir, "settings", "landcover_class_definitions_mod.xlsx")
    ) %>%
      ice_class_definitions()

    # Class 140 affects calibration and RDA, not just allocation. Old products
    # must be rebuilt by preprocessing before the new flow can use them.
    preprocessing_manifest <- file.path(
      input_dir, "summaryResults", "ice_preprocessing.rds"
    )
    expected_preprocessing <- ice_preprocessing_signature(data_dir)
    if (!file.exists(preprocessing_manifest) ||
        !identical(readRDS(preprocessing_manifest), expected_preprocessing)) {
      stop("Preprocessing is outdated; run scripts/01_Preprocessing_LS.R first.")
    }

    input_files <- c(
      file.path(input_dir, "summaryResults", c(
        "initRast.rda", "psaLCover_new.rda", "psaVeg.rda"
      )),
      file.path(input_dir, "psaOvlp", "psaOvlp.rda"),
      file.path(input_dir, "psaEnv", "env_stars.rda")
    )
    if (!all(file.exists(input_files))) stop("Missing site preprocessing files.")
    signature <- ice_cache_key(list(
      settings = ice_settings,
      source = ice_source$checksum,
      preprocessing = expected_preprocessing,
      inputs = tools::md5sum(input_files),
      flow_code = tools::md5sum(c(
        "scripts/02_PixelFlow_parallel_noOptim_LS.R",
        "scripts/ice_cover.R"
      ))
    ))
    manifest_file <- file.path(flow_dir, "ice_flow_manifest.rds")
    cached_manifest <- if (file.exists(manifest_file)) {
      readRDS(manifest_file)
    } else {
      NULL
    }
    expected_outputs <- c(
      file.path(flow_dir, c(
        "psaPixelFlow_new.rda", "LcovOut_new.rda", "LcovOut_new.tif"
      )),
      file.path(summary_dir, c(
        "ice_cover_log.csv", "vegetation_proportions.csv",
        "flow_comparison_new.png", "changeMap_1_new.png",
        "class_140_status.csv"
      ))
    )
    if (is.list(cached_manifest) &&
        identical(cached_manifest$signature, signature) &&
        all(file.exists(c(
          expected_outputs, cached_manifest$masks, cached_manifest$layers
        )))) {
      message("PSA ", ID, ": reusing compatible ice-aware pixel flow")
      return(data.frame(ID = ID, status = "cached", message = ""))
    }
    if (file.exists(manifest_file) && !file.remove(manifest_file)) {
      stop("Could not invalidate the incomplete flow cache manifest.")
    }
    stale_maps <- list.files(
      summary_dir,
      pattern = "^changeMap_[0-9]+_new[.]png$",
      full.names = TRUE
    )
    if (length(stale_maps) && any(!file.remove(stale_maps))) {
      stop("Could not invalidate outdated map figures.")
    }

    load(input_files[1]) # initRast
    load(input_files[2]) # psaLCover
    load(input_files[3]) # psaVeg
    load(input_files[4]) # psaOvlp
    load(input_files[5]) # env_stars

    lcovs <- psaOvlp@lcov
    class_ids <- as.numeric(sub("^LC_", "", lcovs))

    site_predictions <- filter(psaLCover@landcover_ts, Dataset_ID == ID)
    oldest_ice_age <- ice_settings$age_origin - ice_settings$age_interval
    psaFlow <- site_predictions %>%
      filter(Age > 500, Age <= oldest_ice_age) %>%
      arrange(Age)
    age_notes <- character()
    older_count <- sum(site_predictions$Age > oldest_ice_age, na.rm = TRUE)
    if (older_count > 0L) {
      age_notes <- paste0("Capped pixel flow at ", oldest_ice_age,
                          " years. Excluded ", older_count, " older target rows.")
    }
    if (!nrow(psaFlow)) {
      psaFlow <- filter(site_predictions, Age == 500)
      if (!nrow(psaFlow)) stop("No vegetation targets in the supported ice age range.")
      age_notes <- c(age_notes, "Used the 500-year prediction as the only historical target.")
    }
    if (anyDuplicated(psaFlow$Age)) {
      stop("Duplicate vegetation targets for a reconstruction age.")
    }
    target_columns <- names(psaFlow)[!names(psaFlow) %in% c("Dataset_ID", "Age")]
    if (!setequal(target_columns, as.character(class_ids))) {
      stop("Vegetation prediction and RDA class sets differ; rebuild preprocessing.")
    }
    targets <- as.matrix(psaFlow[, as.character(class_ids), drop = FALSE])
    if (any(!is.finite(targets)) || any(targets < 0) ||
        any(rowSums(targets) <= 0)) {
      stop("Invalid vegetation targets.")
    }
    targets <- 100 * targets / rowSums(targets)

    output_ages <- c(-1000, psaFlow$Age)
    ice_ages <- c(0, psaFlow$Age)
    layer_lookup <- ice_layer_lookup(ice_settings$file, ice_ages, ice_settings)

    template <- as(initRast[1], "SpatRaster")
    base_values <- terra::values(template, mat = FALSE)
    cell_ids <- vegetation_cell_ids(base_values, class_ids, ice_settings)
    if (!length(cell_ids)) stop("No reconstructable vegetation pixels.")

    # Terra cell IDs are retained explicitly through extraction, masks, and
    # export; do not assume sf/stars and terra use the same flattened ordering.
    coordinates <- terra::xyFromCell(template, cell_ids)
    init_crds <- st_as_sf(
      data.frame(
        pxlID = cell_ids,
        cellID = cell_ids,
        Landcover = base_values[cell_ids],
        x = coordinates[, 1],
        y = coordinates[, 2]
      ),
      coords = c("x", "y"),
      crs = terra::crs(template)
    )

    environment <- as(env_stars, "SpatRaster")
    init_env <- terra::extract(
      environment,
      terra::vect(init_crds),
      ID = FALSE
    )
    for (column in names(init_env)) {
      missing <- !is.finite(init_env[[column]])
      replacement <- median(init_env[[column]][!missing], na.rm = TRUE)
      if (!is.finite(replacement)) {
        stop("No valid modern environmental values for ", column)
      }
      init_env[[column]][missing] <- replacement
    }

    pixel_scores <- as.data.frame(predict(
      psaOvlp@rda_mod[[1]],
      newdata = init_env,
      type = "wa",
      scaling = 1
    ))[, c("RDA1", "RDA2"), drop = FALSE]

    # Modern ice pixels have no observed vegetation state. Seed once from the
    # closest stored class centroid, and freeze that state while ice-covered.
    initial_state <- init_crds$Landcover
    needs_seed <- !initial_state %in% class_ids
    initial_state[needs_seed] <- nearest_vegetation_class(
      pixel_scores[needs_seed, , drop = FALSE],
      psaOvlp@centers[, c("RDA1", "RDA2"), drop = FALSE],
      class_ids
    )

    lichen_file <- file.path(summary_dir, "class_140_status.csv")
    lichen_status <- read.csv(file.path(
      input_dir, "summaryResults", "class_140_status.csv"
    ))
    lichen_status$retained <- FALSE
    lichen_status$stage <- "class_policy"
    lichen_status$reason <-
      "Excluded by policy; modern lichens/mosses remain static beneath ice"
    lichen_status$converted_pixels <- 0L
    lichen_status$static_pixels <- sum(base_values == 140, na.rm = TRUE)
    write.csv(lichen_status, lichen_file, row.names = FALSE)

    # Preserve the existing environmental density and distance cost ingredients.
    densities <- as(psaOvlp@densities, "SpatRaster")
    density_order <- match(lcovs, names(densities))
    if (anyNA(density_order)) stop("RDA density layer identifiers do not align.")
    densities <- densities[[density_order]]
    density_values <- as.matrix(terra::extract(
      densities,
      as.matrix(pixel_scores)
    ))
    density_values[!is.finite(density_values)] <- 0
    environmental_cost <- 100 - scales::rescale(density_values, c(0, 100))

    density_peaks <- t(vapply(seq_along(class_ids), function(column) {
      peak <- which.max(terra::values(densities[[column]], mat = FALSE))
      as.numeric(terra::xyFromCell(densities, peak))
    }, numeric(2)))
    distance_values <- vapply(seq_along(class_ids), function(column) {
      sqrt(rowSums(sweep(as.matrix(pixel_scores), 2, density_peaks[column, ], "-")^2))
    }, numeric(nrow(pixel_scores)))
    pixel_costs <- environmental_cost +
      0.05 * scales::rescale(distance_values, c(0, 100))
    class_costs <- 1 - psaOvlp@overlapp[, , 2]

    mask_files <- vapply(ice_ages, function(age) {
      aligned_ice_mask(
        ice_settings,
        age,
        template,
        file.path(dir_out, "iceMasks"),
        ice_source$checksum
      )
    }, character(1))
    initial_ice <- ice_mask_values(mask_files[1])
    previous_ice <- initial_ice[cell_ids]
    previous_grid_ice <- initial_ice
    exposed_ever <- rep(FALSE, length(cell_ids))
    grid_exposed_ever <- rep(FALSE, length(base_values))

    source_point <- psaVeg@psas %>% filter(Dataset_ID == ID)
    source_polygon <- st_buffer(
      source_point,
      source_point[["Pollen_Source_Radius [m]"]]
    )
    source_geometry <- prepare_ice_source_geometry(
      template, source_point, source_polygon
    )

    classFlow <- matrix(initial_state, nrow = length(cell_ids), ncol = length(ice_ages))
    diagnostics <- vector("list", length(ice_ages))
    comparisons <- vector("list", length(ice_ages))
    output_files <- character(length(ice_ages))

    for (step in seq_along(ice_ages)) {
      # Every exported cell needs a known mask, including water and background.
      ice <- ice_mask_values(mask_files[step])
      current_ice <- ice[cell_ids]
      eligible <- !current_ice
      switches <- 0L
      newly_exposed <- rep(FALSE, length(cell_ids))
      newly_covered <- rep(FALSE, length(cell_ids))
      grid_newly_exposed <- rep(FALSE, length(base_values))
      grid_newly_covered <- rep(FALSE, length(base_values))

      if (step > 1) {
        newly_exposed <- previous_ice & !current_ice
        newly_covered <- !previous_ice & current_ice
        exposed_ever <- exposed_ever | newly_exposed
        grid_newly_exposed <- previous_grid_ice & !ice
        grid_newly_covered <- !previous_grid_ice & ice
        grid_exposed_ever <- grid_exposed_ever | grid_newly_exposed
        allocation <- ice_free_pixel_flow(
          previous = classFlow[, step - 1],
          eligible = eligible,
          class_ids = class_ids,
          target_shares = targets[step - 1, ],
          class_costs = class_costs,
          pixel_costs = pixel_costs
        )
        classFlow[, step] <- allocation$state
        switches <- allocation$switches
      }

      overlap <- ice_source_overlap(ice, source_geometry)
      eligible_count <- sum(eligible)
      counts <- tabulate(
        match(classFlow[eligible, step], class_ids),
        length(class_ids)
      )
      proportions <- if (eligible_count > 0) {
        100 * counts / eligible_count
      } else {
        rep(NA_real_, length(class_ids))
      }
      requested <- if (step > 1) targets[step - 1, ] else rep(NA_real_, length(class_ids))

      comparisons[[step]] <- data.frame(
        Dataset_ID = ID,
        Age = output_ages[step],
        ice_age = ice_ages[step],
        class = class_ids,
        requested_percent = requested,
        realized_percent = proportions,
        pixels = counts,
        eligible_pixels = eligible_count
      )
      diagnostics[[step]] <- data.frame(
        Dataset_ID = ID,
        Age = output_ages[step],
        ice_age = ice_ages[step],
        ice_layer = layer_lookup$layer_name[step],
        modeled_pixels = length(cell_ids),
        ice_free_pixels = eligible_count,
        ice_covered_pixels = sum(current_ice),
        unknown_pixels = sum(is.na(ice)),
        newly_covered_pixels = sum(newly_covered),
        newly_exposed_pixels = sum(newly_exposed),
        unique_exposed_pixels = sum(exposed_ever),
        grid_pixels = length(base_values),
        grid_ice_covered_pixels = sum(ice),
        grid_newly_covered_pixels = sum(grid_newly_covered),
        grid_newly_exposed_pixels = sum(grid_newly_exposed),
        grid_unique_exposed_pixels = sum(grid_exposed_ever),
        switches = switches,
        water_ice_overlap_pixels = sum(ice[which(base_values == 210)], na.rm = TRUE),
        lichen_ice_overlap_pixels = sum(ice[which(base_values == 140)]),
        source_point_ice = overlap$source_point_ice,
        source_area_ice_fraction = overlap$source_area_ice_fraction
      )
      previous_ice <- current_ice
      previous_grid_ice <- ice

      output <- template
      terra::values(output) <- ice_output_values(
        base_values, cell_ids, classFlow[, step], ice, ice_settings
      )
      names(output) <- paste0("LCov_", output_ages[step])
      output_files[step] <- file.path(flow_dir, paste0("LCov_", output_ages[step], ".tif"))
      terra::writeRaster(output, output_files[step], overwrite = TRUE)
    }

    diagnostics <- bind_rows(diagnostics)
    comparisons <- bind_rows(comparisons)
    write.csv(
      diagnostics,
      file.path(summary_dir, "ice_cover_log.csv"),
      row.names = FALSE
    )
    write.csv(
      comparisons,
      file.path(summary_dir, "vegetation_proportions.csv"),
      row.names = FALSE
    )

    comparison_plot <- comparisons %>%
      pivot_longer(
        c(requested_percent, realized_percent),
        names_to = "Type",
        values_to = "P"
      ) %>%
      ggplot(aes(x = P, y = Age, color = Type)) +
      geom_path(na.rm = TRUE) +
      facet_wrap(~class)
    ggsave(
      file.path(summary_dir, "flow_comparison_new.png"),
      plot = comparison_plot,
      width = 30,
      height = 22,
      units = "cm"
    )

    psaPixelFlow <- list(
      crds = init_crds %>% dplyr::select(-Landcover),
      flow = as_tibble(classFlow, .name_repair = "minimal") %>%
        setNames(paste0("LCov_", output_ages)),
      ice_masks = mask_files,
      ice_lookup = layer_lookup,
      initial_vegetation = initial_state,
      source = ice_source,
      settings = ice_settings,
      signature = signature
    )
    save(psaPixelFlow, file = file.path(flow_dir, "psaPixelFlow_new.rda"))

    output_stack <- terra::rast(output_files)
    terra::writeRaster(
      output_stack,
      file.path(flow_dir, "LcovOut_new.tif"),
      overwrite = TRUE
    )
    LcovOut <- stars::st_as_stars(output_stack) %>% setNames("LandcoverChange")
    save(LcovOut, file = file.path(flow_dir, "LcovOut_new.rda"))
    plot_landcover_flow(LcovOut, class_defs, dir_out)

    # Write last: a partially completed run must never look reusable.
    saveRDS(
      list(signature = signature, masks = mask_files, layers = output_files),
      manifest_file
    )
    message("PSA ", ID, ": completed ice-aware pixel flow")
    data.frame(ID = ID, status = "ok", message = paste(age_notes, collapse = " "))
  }, error = function(error) {
    message("PSA ", ID, ": ", conditionMessage(error))
    data.frame(ID = ID, status = "failed", message = conditionMessage(error))
  })
}

# Workers open the NetCDF and mask files themselves. Terra objects contain
# external pointers and must not be exported from the parent R process.
if (!nrow(rda_files)) stop("No sites available for pixel flow.")
cl <- .make_worker_cluster(nrow(rda_files))
results <- tryCatch({
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(sf)
    library(stars)
    library(ggplot2)
    library(vegan)
    library(terra)
  })
  sf::sf_use_s2(FALSE)
  source("scripts/LandCoverChange_ElasticNet.R")
  source("scripts/ice_cover.R")
  configure_spatial_runtime()
  NULL
})
parallel::clusterExport(
  cl,
  c(
    "project_dir", "rda_dir", "data_dir", "rda_files",
    "ice_settings", "ice_source", "site_pixel_flow", "plot_landcover_flow"
  )
)
parallel::clusterApplyLB(cl, seq_len(nrow(rda_files)), site_pixel_flow)
}, finally = parallel::stopCluster(cl))
results <- bind_rows(results)
if (!is.null(previous_flow_status)) {
  results <- bind_rows(
    previous_flow_status[!previous_flow_status$ID %in% results$ID, , drop = FALSE],
    results
  ) %>% arrange(ID)
}
status_tmp <- tempfile("ice-flow-status-", tmpdir = dirname(flow_status_file))
write.csv(results, status_tmp, row.names = FALSE)
if (!file.rename(status_tmp, flow_status_file)) {
  unlink(status_tmp)
  stop("Could not publish pixel-flow status: ", flow_status_file)
}
if (any(results$status == "failed")) {
  stop("Some sites failed; see landcoverchange/ice_flow_status.csv.")
}
