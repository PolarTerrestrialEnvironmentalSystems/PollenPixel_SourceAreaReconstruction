# Read saved metrics only; never source the pipeline, predict, or refit models.
args <- commandArgs(trailingOnly = TRUE)
manifest <- read.csv(args[1], stringsAsFactors = FALSE)
rows <- list()
digest <- function(x) {
  if (is.null(x)) return("")
  path <- tempfile()
  on.exit(unlink(path))
  saveRDS(x, path, version = 2)
  unname(tools::md5sum(path))
}
for (i in seq_len(nrow(manifest))) {
  entry <- manifest[i, ]
  tryCatch({
    env <- new.env(parent = baseenv())
    load(entry$path, envir = env)
    if (!exists("veg_lc_model", env, inherits = FALSE)) stop("Missing veg_lc_model object")
    models <- env$veg_lc_model
    if (!is.list(models) || !length(models)) stop("Empty or unsupported model list")
    model_names <- names(models)
    if (is.null(model_names)) model_names <- paste0("calibMap_", seq_along(models))
    if (anyNA(model_names) || any(!nzchar(model_names)) || anyDuplicated(model_names)) {
      stop("Invalid or duplicated model names")
    }
    for (j in seq_along(models)) {
      m <- models[[j]]
      modern <- is.data.frame(m$overall_metrics) && is.data.frame(m$per_class_metrics)
      classes <- if (modern) as.character(m$per_class_metrics$class) else names(m$mean_mae_per_class)
      metadata <- list(
        method = if (modern) as.character(m$overall_metrics$method[1]) else "legacy",
        class_set = paste(sort(classes), collapse = ";"),
        samples = if (is.data.frame(m$samplewise_metrics)) nrow(m$samplewise_metrics) else NA_integer_,
        folds = if (is.data.frame(m$samplewise_metrics)) length(unique(m$samplewise_metrics$fold)) else length(m$folds),
        pseudocount = if (is.null(m$pseudocount)) NA_real_ else m$pseudocount,
        calibration_digest = digest(m$model_object$calibration_data),
        fold_digest = digest(m$samplewise_metrics$fold)
      )
      emit <- function(scope, class, metric, value, direction) {
        if (!is.numeric(value) || length(value) != 1) stop(paste("Invalid metric", metric))
        rows[[length(rows) + 1L]] <<- data.frame(
          version = entry$version, file = entry$file, model = model_names[j],
          scope = scope, class = class, metric = metric, value = value,
          direction = direction, status = "ok", message = "",
          as.data.frame(metadata, stringsAsFactors = FALSE), stringsAsFactors = FALSE
        )
      }
      if (modern) {
        if (nrow(m$overall_metrics) != 1L || anyDuplicated(classes)) stop("Invalid metric table")
        for (metric in c("mean_RMSE", "mean_MAE", "mean_R2", "mean_sample_RMSE", "mean_Aitchison")) {
          if (metric %in% names(m$overall_metrics)) {
            emit("overall", "", metric, m$overall_metrics[[metric]][1],
                 if (metric == "mean_R2") "higher" else "lower")
          }
        }
        for (k in seq_len(nrow(m$per_class_metrics))) {
          for (metric in intersect(c("RMSE", "MAE", "R2"), names(m$per_class_metrics))) {
            emit("class", classes[k], metric, m$per_class_metrics[[metric]][k],
                 if (metric == "R2") "higher" else "lower")
          }
        }
      } else {
        # Legacy metrics retain their own names: scales/evaluation differ from
        # the ElasticNet composition metrics and must not be silently equated.
        if (is.null(m$mean_mae_per_class)) stop("Unsupported saved metric structure")
        if (!is.null(m$mean_mae)) emit("overall", "", "legacy_mean_MAE", m$mean_mae, "lower")
        for (class in classes) emit("class", class, "legacy_MAE", m$mean_mae_per_class[[class]], "lower")
      }
    }
  }, error = function(error) {
    rows[[length(rows) + 1L]] <<- data.frame(
      version = entry$version, file = entry$file, model = "", scope = "", class = "",
      metric = "", value = NA_real_, direction = "", status = "error",
      message = conditionMessage(error), method = "", class_set = "", samples = NA,
      folds = NA, pseudocount = NA, calibration_digest = "", fold_digest = ""
    )
  })
}
if (!length(rows)) stop("No saved regression metrics found")
write.csv(do.call(rbind, rows), args[2], row.names = FALSE, na = "")
