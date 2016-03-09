#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#
# This software was authored by Zhian N. Kamvar and Javier F. Tabima, graduate 
# students at Oregon State University; Jonah C. Brooks, undergraduate student at
# Oregon State University; and Dr. Nik Grünwald, an employee of USDA-ARS.
#
# Permission to use, copy, modify, and distribute this software and its
# documentation for educational, research and non-profit purposes, without fee, 
# and without a written agreement is hereby granted, provided that the statement
# above is incorporated into the material, giving appropriate attribution to the
# authors.
#
# Permission to incorporate this software into commercial products may be
# obtained by contacting USDA ARS and OREGON STATE UNIVERSITY Office for 
# Commercialization and Corporate Development.
#
# The software program and documentation are supplied "as is", without any
# accompanying services from the USDA or the University. USDA ARS or the 
# University do not warrant that the operation of the program will be 
# uninterrupted or error-free. The end-user understands that the program was 
# developed for research purposes and is advised not to rely exclusively on the 
# program for any reason.
#
# IN NO EVENT SHALL USDA ARS OR OREGON STATE UNIVERSITY BE LIABLE TO ANY PARTY 
# FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
# LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, 
# EVEN IF THE OREGON STATE UNIVERSITY HAS BEEN ADVISED OF THE POSSIBILITY OF 
# SUCH DAMAGE. USDA ARS OR OREGON STATE UNIVERSITY SPECIFICALLY DISCLAIMS ANY 
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE AND ANY STATUTORY 
# WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
# BASIS, AND USDA ARS AND OREGON STATE UNIVERSITY HAVE NO OBLIGATIONS TO PROVIDE
# MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#==============================================================================#

#==============================================================================#
# Message to print after running the poppr function
#
# Public functions utilizing this function:
# # poppr poppr.all
#
# Internal functions utilizing this function:
# # none
#==============================================================================#

poppr_message <- function(){
  msg <- paste0(
  "-------------------------------------------------------------------------|\n",
  "Pop     = Population name (Total == Pooled)\n",
  "N       = Census population size\n",
  "MLG     = Number of unique multilocus genotypes (MLG) observed\n",
  "eMLG    = Number of expected MLG based on rarefaction at smallest N >= 10\n",
  "SE      = Standard error of rarefaction analysis\n",
  "H       = Shannon-Wiener Index of MLG diversity\n",
  "G       = Stoddart and Taylor's Index of MLG diversity\n",
  "lambda  = Simpson's index\n",
  "E.5     = Evenness\n",
  "Hexp    = Nei's 1978 expected heterozygosity\n",
  "Ia      = Index of association\n",
  "rbarD   = Standardized index of association\n",
  "-------------------------------------------------------------------------|\n"
  )
  message(msg)
}

#==============================================================================#
# A function that will quit the function if a level in the hierarchy is not
# present in the given data frame.
#
# Public functions utilizing this function:
# # setPop strata poppr.amova
#
# Internal functions utilizing this function:
# # make_hierarchy make_ade_df
#==============================================================================#
hier_incompatible_warning <- function(levs, df){
  msg <- paste("One or more levels in the given hierarchy is not present", 
               "in the data frame.",
               "\nHierarchy:\t", paste(levs, collapse = ", "), "\nData:\t\t", 
               paste(names(df), collapse = ", "))
  return(msg)
}

#==============================================================================#
# Warning message for when a distance matrix is non-euclidean and the user 
# did not specify an appropriate correction.
#
# Public functions utilizing this function:
# # poppr.amova
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
not_euclid_msg <- function(correction){
  msg <- paste0("\nThe distance matrix generated is non-euclidean and a ", 
  							"correction is needed.\n",
                "You supplied: correction = '", correction, "'\nPlease change",
                " it to one of the following:\n",
                "\t'cailliez'\t'quasieuclid'\t'lingoes'")
  return(msg)
}

#==============================================================================#
# Warning message for the function popsub.
# Public functions utilizing this function:
#
# # popsub
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
unmatched_pops_warning <- function(pops, sublist){
  msg <- paste("The sublist provided does not match any of the populations:\n",
               "\tsublist.......", paste(sublist, collapse = " "), "\n", 
               "\tPopulations...", paste(pops, collapse = " "))
  return(msg)
}

#==============================================================================#
# Warning messages for Bruvo's distance calculation.
# Public functions utilizing this function:
#
# # bruvo.dist bruvo.boot bruvo.msn
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
repeat_length_warning <- function(replen){
	msg <- paste("\n\nRepeat length vector for loci is not equal to the number",
							 "of loci represented.\nEstimating repeat lengths from data:\n",
							 paste0("c(", paste(replen, collapse = ", "),")"), "\n\n")
  return(msg)
}

non_ssr_data_warning <- function(){
	msg <- paste("\nThis dataset does not appear to be microsatellite data.",
				       "Bruvo's Distance can only be applied for true microsatellites.")
	return(msg)
}

#==============================================================================#
# Warning message for Neighbor-Joining trees.
# Public functions utilizing this function:
#
# # aboot bruvo.boot
#
# Internal functions utilizing this function:
# # none
#==============================================================================#
negative_branch_warning <- function(){
	msg <- paste("Some branch lengths of the tree are negative.", 
							 "Normalizing branches according to Kuhner and Felsenstein",
							 "(1994)")
	return(msg)
}

#==============================================================================#
# Warning message for mlg.crosspop with the flag mlgsub.
# Public functions utilizing this function:
#
# # mlg.crosspop
#
# Internal functions utilizing this function:
# # none
#==============================================================================#

mlg_sub_warning <- function(mlgs){
  msg <- paste0("The following multilocus genotypes are not defined in this ",
                "dataset: ", paste(mlgs, collapse = ", "))
  return(msg)
}


#==============================================================================#
# create a message about missing data
#
# param things a character vector of names that are removed
# param type a vector of length two giving the singular and plural of things
# param nremoved the number of removed items
# param cutoff the cutoff at which items were removed
# 
# Public functions utilizing this function:
# ## missingno
# 
# Internal functions utilizing this function:
# ## none
#==============================================================================#
missing_messenger <- function(things, type = c("locus", "loci"), nremoved = 1, 
                              cutoff = 0.05){
  type <- ifelse(length(things) == 1, type[1], type[2])
  cutoff <- cutoff*100
  msg <- paste0("\nFound ", nremoved, " missing values.\n\n",
                length(things), " ", type, 
                " contained missing values greater than ", cutoff, "%\n\n",
                "Removing ", length(things), " ", type, ":\n", 
                format_char_width(things, width = getOption("width") - 10))
  message(msg)
}
