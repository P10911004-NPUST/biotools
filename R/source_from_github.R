suppressMessages({
    if (!require(devtools)) install.packages("devtools")
})

source_from_github <- function(url){
  url_0 <- paste0(url, "?raw=TRUE")
  devtools::source_url(url = url_0)
}