# Run from the repository root: Rscript tests/workers.R
# Real PSOCK workers, without loading spatial inputs or overwriting outputs.
expressions <- parse("scripts/LandCoverChange_ElasticNet.R")
for (expression in expressions) {
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      as.character(expression[[2]]) %in% c(".worker_count", ".make_worker_cluster", ".psock_lapply")) {
    eval(expression)
  }
}
expect_error <- function(expression) {
  stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
}
for (invalid in list(NULL, NA_real_, Inf, NaN, 0, -1, 1.5, "2", c(1, 2), 2^32)) {
  NODES <- invalid
  expect_error(.worker_count())
}
rm(NODES)
expect_error(.worker_count())

for (limit in c(1L, 2L)) {
  NODES <- limit
  cluster <- .make_worker_cluster()
  tryCatch({
    stopifnot(length(cluster) == limit)
    pids <- unlist(parallel::clusterCall(cluster, Sys.getpid))
    stopifnot(length(unique(pids)) == limit, !Sys.getpid() %in% pids)
  }, finally = parallel::stopCluster(cluster))
}
NODES <- 2L
cluster <- .make_worker_cluster(1L)
tryCatch(stopifnot(length(cluster) == 1L), finally = parallel::stopCluster(cluster))
stopifnot(identical(.psock_lapply(list(), identity), list()))
stopifnot(identical(.psock_lapply(1:4, function(x, increment) x + increment,
                                increment = 3L), as.list(4:7)))

# Failures in dispatched work must close all connections to the worker pool.
before <- nrow(showConnections(all = TRUE))
expect_error(.psock_lapply(1:2, function(x) stop("test worker failure")))
stopifnot(nrow(showConnections(all = TRUE)) == before)

# Pixel-flow startup clears old objects while retaining the configured limit.
flow <- parse("scripts/02_PixelFlow_parallel_noOptim_LS.R")
stage <- new.env(parent = globalenv())
stage$NODES <- 2L
stage$old_object <- TRUE
eval(flow[[1]], envir = stage)
stopifnot(identical(stage$NODES, 2L), !exists("old_object", stage, inherits = FALSE))
cat("Worker configuration, task limits, dispatch, cleanup, and stage preservation passed.\n")
