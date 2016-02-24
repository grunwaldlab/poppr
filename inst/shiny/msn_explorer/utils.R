get_dist <- function(indist){
  indist <- switch(indist,
      Dissimilarity = "diss.dist",
      Bruvo         = "bruvo.dist",
      Nei           = "nei.dist",
      Rogers        = "rogers.dist",
      Edwards       = "edwards.dist",
      Provesti      = "provesti.dist",
      Reynolds      = "reynolds.dist"
  )
  return(indist)
}
#------------------------------------------------------------------------------#
# If the user does select custom, this function will make sure that it is 
# encapsulated in parentheses. This makes sure that the entire expression is
# evaluated.
#------------------------------------------------------------------------------#
parse_distfun <- function(x){
  if (grepl("function", x)){
    x <- paste0("(", x, ")")
  }
  return(x)
}
#------------------------------------------------------------------------------#
# This is a function that will print vectors of numerics or characters in a way
# that can be run directly from R. For example, if the user checks two boxes
# labeling the populations "pop1" and "pop2", that vector of populations gets 
# passed onto the popsub function. By default, it will print like so:
#
# [1] "pop1" "pop2"
#
# You cannot copy and paste this into R because it will throw an error, thus, 
# this function will print the vector like so:
#
# c("pop1", "pop2")
#
# This is then usable in the Command tab for the popsub function.
#------------------------------------------------------------------------------#
make_dput <- function(x){
  return(capture.output(dput(x)))
}

#------------------------------------------------------------------------------#
# A function to search the user's global environment and grab all of the useable
# objects. In this case, it's genind and genclone objects, but the objclass
# argument allows this to be extensible to any class. This is immensely useful
# so that the user does not have to save their objects as rda files, nor do they
# have to save them as text files for input.
#------------------------------------------------------------------------------#
get_globals <- function(objclass = c("genind", "genclone")){
  # Grab all object names in users R session.
  myobjs <- ls(envir = .GlobalEnv) 
  if (length(myobjs) == 0){
    return(myobjs)
  }
  # Go through each name and test if it is any of the classes in objclass.
  gens <- vapply(myobjs, FUN = is_usable, FUN.VALUE = logical(1), objclass)
  myobjs[gens]
}

#------------------------------------------------------------------------------#
# This function tests a single object name to see if it is of class objclass.
# The function is used in get_globals
#------------------------------------------------------------------------------#
is_usable <- function(object, objclass = c("genind", "genclone", "genlight", "snpclone")){
  inherits(get(object, .GlobalEnv), objclass)
}