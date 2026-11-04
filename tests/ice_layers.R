# Rscript tests/ice_layers.R -- naming compatibility and download restart tests.
# Uses synthetic labels and mocked downloads; never accesses the network.
source("scripts/ice_cover.R")
expect_error <- function(expression) {
  stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
}
ordinal <- paste0("ICEM_", 1:261)
coordinate <- paste0("ICEM_T122KP11=", sprintf("%.17g", seq(-26, 0, by = 0.1)))
coordinate[2] <- "ICEM_T122KP11=-25.89999999999999"
coordinate[260] <- "ICEM_T122KP11=-0.09999999999999432"
ages <- c(0, 2000, 20000, 26000)
expected <- ice_layer_lookup_names(ordinal, ages)
stopifnot(
  identical(expected, ice_layer_lookup_names(coordinate, ages)),
  identical(expected$layer_number, c(261L, 241L, 61L, 1L))
)
for (labels in list(ordinal[-261], c(ordinal[-1], "ICEM_2"), coordinate[-261],
                    rev(coordinate), replace(coordinate, 2, coordinate[1]),
                    replace(coordinate, 2, "ICEM_T122KP11=-25.85"),
                    replace(coordinate, 2, "ICEM_T122KP11=NaN"),
                    sub("ICEM", "OTHER", coordinate))) {
  expect_error(ice_layer_lookup_names(labels, ages))
}
expect_error(ice_layer_lookup_names(coordinate, 55))
expect_error(ice_layer_lookup_names(coordinate, 30000))

# Execute the actual download/promotion/manifest control flow with injected
# validation and downloader functions. Files contain fixture text, not NetCDF.
local({
  directory <- tempfile("ice-restart-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  for (initial in c("valid", "invalid", "absent")) {
    fixture <- new.env(parent = globalenv())
    fixture$prepare <- prepare_ice_cover
    environment(fixture$prepare) <- fixture
    fixture$downloads <- 0L
    fixture$validate_ice_cover_file <- function(filename, settings) {
      if (!identical(readLines(filename), "valid")) stop("Invalid fixture download")
      invisible(TRUE)
    }
    fixture$ice_layer_lookup <- function(...) data.frame(layer_number = 261L)
    fixture$download.file <- function(url, destfile, ...) {
      fixture$downloads <- fixture$downloads + 1L
      writeLines("valid", destfile)
      0L
    }
    settings <- list(file = file.path(directory, paste0(initial, ".nc")), url = "fixture")
    partial <- paste0(settings$file, ".part")
    if (initial != "absent") writeLines(initial, partial)
    manifest <- fixture$prepare(settings)
    stopifnot(
      fixture$downloads == if (initial == "valid") 0L else 1L,
      file.exists(settings$file), !file.exists(partial),
      identical(manifest$checksum, unname(tools::md5sum(settings$file)))
    )
    # A subsequent call reuses the promoted file and its manifest.
    stopifnot(identical(fixture$prepare(settings), manifest))
  }
})
cat("Ordinal/coordinate naming, invalid axes, and download restart checks passed.\n")
