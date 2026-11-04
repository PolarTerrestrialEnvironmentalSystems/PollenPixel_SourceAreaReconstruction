# Ice masks prescribe the available land. Vegetation state is stored separately
# so that pixels exposed by retreating ice can rejoin vegetation reconstruction.

configure_spatial_runtime <- function() {
  if (.Platform$OS.type != "windows") return(invisible(NULL))

  # Other applications (for example PostGIS) can set machine-wide paths to an
  # incompatible proj.db. Prefer the data bundled with this terra installation.
  # Changes affect only this R process and workers it starts, not Windows settings.
  proj_directory <- system.file("proj", package = "terra")
  gdal_directory <- system.file("gdal", package = "terra")

  if (file.exists(file.path(proj_directory, "proj.db"))) {
    Sys.setenv(PROJ_DATA = proj_directory, PROJ_LIB = proj_directory)
    terra::projPaths(unique(c(proj_directory, terra::projPaths())))
  }
  if (nzchar(gdal_directory) && dir.exists(gdal_directory)) {
    Sys.setenv(GDAL_DATA = gdal_directory)
    terra::setGDALconfig("GDAL_DATA", gdal_directory)
  }
  if (!terra::proj_ok()) {
    stop("PROJ cannot read its coordinate-system database in this R process.")
  }

  invisible(NULL)
}

ice_cover_settings <- function(data_dir = "input") {
  filename <- "TOPicemsk.GLACD26kN9894GE90227A6005GGrBgic.nc"

  list(
    version = "ice-cover-v2-static-140",
    file = file.path(data_dir, "data", "global_maps", "ice_cover", filename),
    url = paste0(
      "https://sharebox.lsce.ipsl.fr/index.php/s/yfuUw91ruuJXroC/",
      "download?path=%2F&files=", filename
    ),
    age_origin = 26100,
    age_interval = 100,
    ice_threshold = 0,
    ice_class = 220,
    static_classes = c(0, 140, 200, 201, 202, 205, 210, 9999),
    reconstructable_classes = 220
  )
}

ice_cache_key <- function(object) {
  temporary_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temporary_file), add = TRUE)
  saveRDS(object, temporary_file, version = 2)
  unname(tools::md5sum(temporary_file))
}

ice_class_definitions <- function(class_defs) {
  # Lichens/mosses stay at their modern positions. Static map classes and ice
  # are never vegetation targets, even if the workbook is edited otherwise.
  lichens <- class_defs$Class_Code == 140
  fixed <- class_defs$Class_Code %in% c(0, 140, 210, 220, 9999)

  class_defs$Merge_To_Class[lichens] <- NA_real_
  class_defs$Keep_Constant_in_Past[lichens] <- TRUE
  class_defs$Set_Zero_In_Past[lichens] <- FALSE
  class_defs$Include_Class_In_Calibration[fixed] <- FALSE
  class_defs$Predicted_In_Past[fixed] <- FALSE
  class_defs$Transfer[class_defs$Class_Code %in% c(140, 210, 220)] <- FALSE
  class_defs
}

ice_layer_lookup <- function(filename, ages, settings = ice_cover_settings()) {
  # The file contains three subdatasets. Select ICEM explicitly; the GDAL time
  # coordinate is not the chronology supplied for this simulation.
  ice <- terra::rast(filename, subds = "ICEM")
  source_index <- (settings$age_origin - ages) / settings$age_interval

  if (any(!is.finite(source_index)) ||
      any(abs(source_index - round(source_index)) > 1e-8)) {
    stop("Reconstruction ages must map to integer ICEM indices.")
  }

  layer_name <- paste0("ICEM_", as.integer(round(source_index)))
  layer_number <- match(layer_name, names(ice))

  if (anyNA(layer_number) || anyDuplicated(names(ice))) {
    stop(
      "Missing or ambiguous ice layers: ",
      paste(layer_name[is.na(layer_number)], collapse = ", ")
    )
  }

  data.frame(
    age = ages,
    source_index = as.integer(round(source_index)),
    layer_name = layer_name,
    layer_number = layer_number
  )
}

prepare_ice_cover <- function(settings = ice_cover_settings()) {
  dir.create(dirname(settings$file), recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(settings$file)) {
    partial_file <- paste0(settings$file, ".part")
    options(timeout = max(1e6, getOption("timeout")))
    download_status <- download.file(
      settings$url,
      destfile = partial_file,
      method = "libcurl",
      mode = "wb"
    )
    if (download_status != 0) stop("Ice download did not complete successfully.")

    # Force reads at both ends so an HTML response or truncated NetCDF fails
    # before it becomes the reusable input file.
    lookup <- ice_layer_lookup(partial_file, c(0, 26000), settings)
    for (layer in lookup$layer_number) {
      raster <- terra::rast(partial_file, subds = "ICEM", lyrs = layer)
      values <- terra::values(raster, mat = FALSE)
      if (!length(values) || all(is.na(values))) {
        stop("Downloaded ice layer contains no readable values.")
      }
    }

    if (!file.rename(partial_file, settings$file)) {
      stop("Could not promote the validated ice download: ", partial_file)
    }
  }

  ice_layer_lookup(settings$file, 0, settings)
  info <- file.info(settings$file)
  identity <- list(size = info$size, modified = as.numeric(info$mtime))
  manifest_file <- paste0(settings$file, ".rds")
  manifest <- if (file.exists(manifest_file)) readRDS(manifest_file) else NULL

  if (is.null(manifest) || !identical(manifest$identity, identity)) {
    endpoints <- ice_layer_lookup(settings$file, c(0, 26000), settings)
    for (layer in endpoints$layer_number) {
      values <- terra::values(
        terra::rast(settings$file, subds = "ICEM", lyrs = layer),
        mat = FALSE
      )
      if (!length(values) || all(is.na(values))) {
        stop("Ice input contains an unreadable layer; replace ", settings$file)
      }
    }
    manifest <- list(
      identity = identity,
      checksum = unname(tools::md5sum(settings$file)),
      source_url = settings$url
    )
    saveRDS(manifest, manifest_file)
  }

  manifest
}

aligned_ice_mask <- function(settings, age, template, cache_dir, source_checksum) {
  lookup <- ice_layer_lookup(settings$file, age, settings)
  grid <- list(
    crs = terra::crs(template),
    extent = as.vector(terra::ext(template)),
    resolution = terra::res(template),
    dimensions = dim(template)[1:2]
  )
  mask_configuration <- list(
    age_origin = settings$age_origin,
    age_interval = settings$age_interval,
    ice_threshold = settings$ice_threshold,
    alignment = "rotated-nearest-v1"
  )
  signature <- ice_cache_key(list(mask_configuration, source_checksum, lookup, grid))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, paste0("ice_", signature, ".tif"))

  if (!file.exists(cache_file)) {
    # Version 1 included vegetation policy in the mask key. Its binary masks
    # remain valid: the source, layer, threshold, and grid are unchanged.
    legacy_settings <- list(
      version = "ice-cover-v1",
      file = settings$file,
      url = settings$url,
      age_origin = settings$age_origin,
      age_interval = settings$age_interval,
      ice_threshold = settings$ice_threshold,
      ice_class = 220,
      constant_classes = c(0, 210, 9999),
      reconstructable_classes = c(140, 220)
    )
    legacy_key <- ice_cache_key(list(legacy_settings, source_checksum, lookup, grid))
    legacy_file <- file.path(cache_dir, paste0("ice_", legacy_key, ".tif"))
    if (file.exists(legacy_file) && terra::compareGeom(
      terra::rast(legacy_file), template, stopOnError = FALSE
    )) {
      cache_file <- legacy_file
    }
  }

  if (!file.exists(cache_file)) {
    ice <- terra::rast(
      settings$file,
      subds = "ICEM",
      lyrs = lookup$layer_number
    )

    if (!nzchar(terra::crs(ice))) {
      stop("Ice source has no coordinate reference system.")
    }
    if (terra::is.lonlat(ice) && terra::xmax(ice) > 180) {
      ice <- terra::rotate(ice)
    }

    # Threshold before resampling: this is a categorical ice-status raster.
    ice <- ice > settings$ice_threshold
    aligned <- terra::project(ice, template, method = "near")
    names(aligned) <- "ice"
    temporary_file <- tempfile(tmpdir = cache_dir, fileext = ".tif")
    on.exit(unlink(temporary_file), add = TRUE)
    terra::writeRaster(aligned, temporary_file, datatype = "INT1U")
    # Release GDAL's file handle before renaming on Windows.
    rm(aligned)
    gc(verbose = FALSE)

    # Each site has its own cache. A second worker may already have written it.
    if (!file.exists(cache_file) && !file.rename(temporary_file, cache_file)) {
      stop("Could not save aligned ice mask: ", cache_file)
    }
  }

  raster <- terra::rast(cache_file)
  if (!terra::compareGeom(raster, template, stopOnError = FALSE)) {
    stop("Cached ice mask is not aligned to the reconstruction grid.")
  }
  cache_file
}

ice_mask_values <- function(filename, required_cells = NULL) {
  values <- terra::values(terra::rast(filename), mat = FALSE)
  checked <- if (is.null(required_cells)) values else values[required_cells]

  if (anyNA(checked) || any(!checked %in% c(0, 1))) {
    stop("Unknown ice status in the reconstruction domain: ", filename)
  }
  values == 1
}

nearest_vegetation_class <- function(pixel_scores, centroids, class_ids) {
  pixel_scores <- as.matrix(pixel_scores)
  centroids <- as.matrix(centroids)
  if (!nrow(pixel_scores)) return(numeric())
  if (ncol(pixel_scores) != ncol(centroids) ||
      nrow(centroids) != length(class_ids) ||
      any(!is.finite(pixel_scores)) || any(!is.finite(centroids))) {
    stop("Cannot initialize vegetation from invalid RDA scores or centroids.")
  }

  # Sorting makes exact ties deterministic, independent of model column order.
  class_order <- order(class_ids)
  centroids <- centroids[class_order, , drop = FALSE]
  class_ids <- class_ids[class_order]
  distances <- vapply(seq_along(class_ids), function(column) {
    rowSums(sweep(pixel_scores, 2, centroids[column, ], "-")^2)
  }, numeric(nrow(pixel_scores)))
  distances <- matrix(distances, nrow = nrow(pixel_scores))
  class_ids[max.col(-distances, ties.method = "first")]
}

vegetation_target_counts <- function(shares, total, class_ids) {
  if (length(shares) != length(class_ids) ||
      any(!is.finite(shares)) || any(shares < 0) || sum(shares) <= 0) {
    stop("Invalid vegetation target proportions.")
  }

  exact <- total * shares / sum(shares)
  counts <- floor(exact)
  remainder <- as.integer(total - sum(counts))
  if (remainder > 0) {
    priority <- order(-(exact - counts), class_ids)
    counts[priority[seq_len(remainder)]] <-
      counts[priority[seq_len(remainder)]] + 1
  }
  as.integer(counts)
}

ice_free_pixel_flow <- function(previous, eligible, class_ids, target_shares,
                               class_costs, pixel_costs) {
  if (length(previous) != length(eligible) || anyNA(eligible) ||
      nrow(pixel_costs) != length(previous) ||
      ncol(pixel_costs) != length(class_ids) ||
      !all(previous %in% class_ids)) {
    stop("Vegetation states, eligibility, and allocation costs do not align.")
  }

  current <- previous
  total <- sum(eligible)
  target_counts <- vegetation_target_counts(target_shares, total, class_ids)
  old_counts <- tabulate(match(previous[eligible], class_ids), length(class_ids))

  if (total == 0 || identical(as.integer(old_counts), target_counts)) {
    return(list(state = current, targets = target_counts, switches = 0L))
  }
  if (any(!is.finite(class_costs)) || any(!is.finite(pixel_costs[eligible, ]))) {
    stop("Nonfinite vegetation transport or pixel allocation costs.")
  }

  # Integer supply/demand avoids floor(percent * pixel_count) residuals.
  flow <- transport::transport(
    old_counts,
    target_counts,
    class_costs,
    fullreturn = FALSE
  )
  if (any(abs(flow$mass - round(flow$mass)) > 1e-6)) {
    stop("Transport returned noninteger pixel counts.")
  }
  flow$mass <- as.integer(round(flow$mass))
  flow <- flow[flow$from != flow$to & flow$mass > 0, , drop = FALSE]

  for (donor in unique(flow$from)) {
    transfers <- flow[flow$from == donor, , drop = FALSE]
    candidates <- which(eligible & previous == class_ids[donor])

    # Pixels prefer their cheapest still-needed destination. If its quota is
    # filled, the remaining pixels can choose again among the remaining classes.
    while (nrow(transfers) > 0) {
      preferences <- max.col(
        -pixel_costs[candidates, transfers$to, drop = FALSE],
        ties.method = "first"
      )
      selected <- integer()

      for (destination in seq_len(nrow(transfers))) {
        available <- candidates[preferences == destination]
        ranking <- order(pixel_costs[available, transfers$to[destination]], available)
        take <- min(length(available), transfers$mass[destination])
        chosen <- available[ranking[seq_len(take)]]
        current[chosen] <- class_ids[transfers$to[destination]]
        transfers$mass[destination] <- transfers$mass[destination] - take
        selected <- c(selected, chosen)
      }

      if (!length(selected)) stop("Insufficient eligible transport donor pixels.")
      candidates <- setdiff(candidates, selected)
      transfers <- transfers[transfers$mass > 0, , drop = FALSE]
    }
  }

  realized <- tabulate(match(current[eligible], class_ids), length(class_ids))
  if (!identical(as.integer(realized), target_counts)) {
    stop("Pixel allocation did not conserve the requested vegetation counts.")
  }

  list(
    state = current,
    targets = target_counts,
    switches = sum(current != previous)
  )
}

vegetation_cell_ids <- function(base_values, class_ids,
                                settings = ice_cover_settings()) {
  if (anyNA(class_ids) || anyDuplicated(class_ids) ||
      any(class_ids %in% c(settings$static_classes, settings$ice_class))) {
    stop("Static classes or ice found in the vegetation class set; rebuild preprocessing.")
  }
  which(base_values %in% c(class_ids, settings$reconstructable_classes))
}

ice_output_values <- function(base_values, cell_ids, vegetation, ice,
                              settings = ice_cover_settings()) {
  if (length(ice) != length(base_values) || anyNA(ice) ||
      length(vegetation) != length(cell_ids)) {
    stop("Output requires an aligned, known ice mask across the complete grid.")
  }
  if (any(is.na(base_values[cell_ids])) ||
      any(base_values[cell_ids] %in% settings$static_classes)) {
    stop("Static pixels cannot receive reconstructed vegetation.")
  }
  base_values[cell_ids] <- vegetation
  # Ice is the last overlay and overrides water, static land, and missing data.
  base_values[ice] <- settings$ice_class
  base_values
}

ice_source_overlap <- function(mask_file, source_point, source_polygon) {
  mask <- terra::rast(mask_file)
  point <- terra::project(terra::vect(source_point), terra::crs(mask))
  polygon <- terra::project(terra::vect(source_polygon), terra::crs(mask))
  point_ice <- terra::extract(mask, point, ID = FALSE)[1, 1]

  # Weight the actual source polygon by cell area and boundary-cell fraction.
  cell_area <- terra::cellSize(mask, unit = "m")
  names(cell_area) <- "area"
  extracted <- terra::extract(c(mask, cell_area), polygon, exact = TRUE)
  if (anyNA(extracted$ice) || is.na(point_ice)) {
    stop("Unknown ice status at the pollen source point or source polygon.")
  }
  weights <- extracted$area * extracted$fraction
  list(
    source_point_ice = point_ice == 1,
    source_area_ice_fraction = sum(extracted$ice * weights) / sum(weights)
  )
}

ice_preprocessing_signature <- function(data_dir = "input") {
  ice_cache_key(list(
    version = ice_cover_settings(data_dir)$version,
    files = tools::md5sum(c(
      "scripts/01_Preprocessing_LS.R",
      "scripts/LandCoverChange_ElasticNet.R",
      file.path(data_dir, "settings", "landcover_class_definitions_mod.xlsx")
    )),
    # Class-definition changes require recalibration. Runtime path fixes and
    # mask/allocation helpers affect flow only, so they must not discard models.
    class_configuration = paste(deparse(body(ice_class_definitions)), collapse = "\n")
  ))
}

ice_preprocessing_cache <- function(dir_out, signature, complete = FALSE) {
  manifest_file <- file.path(dir_out, "summaryResults", "ice_preprocessing.rds")
  if (complete) {
    saveRDS(signature, manifest_file)
    return(invisible(TRUE))
  }
  previous <- if (file.exists(manifest_file)) readRDS(manifest_file) else NULL
  if (identical(previous, signature)) return(invisible(TRUE))

  # Only known derived files are invalidated. Raw data and psaVeg are retained.
  stale_files <- c(
    "summaryResults/lcovShares.rda",
    "summaryResults/veg_lc_model_new.rda",
    "summaryResults/psaLCover_new.rda",
    "summaryResults/rdaOut.rda",
    "summaryResults/initRast.rda",
    "summaryResults/map_init.png",
    "summaryResults/map_init_withHuman.png",
    "summaryResults/rda_summary.png",
    "psaCrds/psaCrds.rda",
    "psaEnv/psaEnv.rda",
    "psaEnv/env_stars.rda",
    "psaEnv/z_tranTab.rda",
    "psaOvlp/psaOvlp.rda"
  )
  existing <- file.path(dir_out, stale_files)
  existing <- existing[file.exists(existing)]
  if (length(existing)) {
    message("Rebuilding class-dependent preprocessing for ", basename(dir_out))
    if (any(!file.remove(existing))) stop("Could not invalidate preprocessing cache.")
  }
  if (file.exists(manifest_file)) file.remove(manifest_file)
  invisible(FALSE)
}
