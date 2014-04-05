### Function to source code from github. The only thing it needs is a URL

source_github <- function(URL) {
  # load package
  require(RCurl)

  # read script lines from website and evaluate
  script <- getURL(URL, ssl.verifypeer = FALSE)
  eval(parse(text = script),envir=.GlobalEnv)
}  

