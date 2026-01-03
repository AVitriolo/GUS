is_empty <- function(f) {
  !file.exists(f) || file.size(f) == 0
}