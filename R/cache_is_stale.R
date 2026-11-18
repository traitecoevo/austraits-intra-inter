# TRUE when a cached result needs (re)computing: the cache file is missing, or
# the model fit it was derived from has been saved since (e.g. after a refit).
cache_is_stale <- function(cache_path, source_path) {
  !file.exists(cache_path) ||
    (file.exists(source_path) && file.mtime(source_path) > file.mtime(cache_path))
}
