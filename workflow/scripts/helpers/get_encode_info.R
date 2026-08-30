
get_encode_info <- function(acc) {
  url <- paste0("https://www.encodeproject.org/files/", acc, "/?format=json")
  res <- GET(url, add_headers(`accept` = "application/json"))
  if(status_code(res) != 200) return(c(epigenome = NA, target = NA))

  data <- content(res, as = "text", encoding = "UTF-8")
  meta <- fromJSON(data)

  epigenome <- meta$biosample_summary
  target <- if (!is.null(meta$target)) meta$target$label else NA

  c(epigenome = epigenome, target = target)
}
