.onAttach <- function(...) {
  packageStartupMessage(paste("This is poppr version", utils::packageVersion("poppr")))
  if (!interactive() || stats::runif(1) > 0.1) return()

  tips <- c(
    "\nNeed help? Try the poppr mailing list: http://groups.google.com/group/poppr.\n",
    "\nUse suppressPackageStartupMessages to eliminate package startup messages.\n"
  )
  
  tip <- sample(tips, 1)
  packageStartupMessage(tip)
}

