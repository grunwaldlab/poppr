




###########################
# Method seploc for genind
###########################
setGeneric("seploc", function(x, ...) standardGeneric("seploc"))

setMethod("seploc", signature(x="genind"), function(x,truenames=TRUE,res.type=c("genind","matrix")){
    if(x@type=="PA"){
        msg <- paste("seploc is not implemented for presence/absence markers")
        cat("\n",msg,"\n")
        return(invisible())
    }


    if(!is.genind(x)) stop("x is not a valid genind object")
    res.type <- match.arg(res.type)
    if(res.type=="genind") { truenames <- TRUE }

    temp <- x@loc.fac
    nloc <- length(levels(temp))
    levels(temp) <- 1:nloc

    kX <- list()

    for(i in 1:nloc){
        kX[[i]] <- matrix(x@tab[,temp==i],ncol=x@loc.nall[i])

        if(!truenames){
            rownames(kX[[i]]) <- rownames(x@tab)
            colnames(kX[[i]]) <- paste(names(x@loc.names)[i],names(x@all.names[[i]]),sep=".")
        }else{
            rownames(kX[[i]]) <- x@ind.names
            colnames(kX[[i]]) <- paste(x@loc.names[i],x@all.names[[i]],sep=".")
        }
    }

    if(truenames) {
        names(kX) <- x@loc.names
    } else{
        names(kX) <- names(x@loc.names)
    }

    prevcall <- match.call()
    if(res.type=="genind"){
        ## ploidy bug fixed by Zhian N. Kamvar
        ##kX <- lapply(kX, genind, pop=x@pop, prevcall=prevcall)
        kX <- lapply(kX, genind, pop=x@pop, prevcall=prevcall, ploidy=x@ploidy, type=x@type)
        for(i in 1:length(kX)){
            kX[[i]]@other <- x@other
        }
    }

    return(kX)
})



###########################
# Method seploc for genpop
###########################
setMethod("seploc", signature(x="genpop"), function(x,truenames=TRUE,res.type=c("genpop","matrix")){
     if(x@type=="PA"){
         msg <- paste("seploc is not implemented for presence/absence markers")
         cat("\n",msg,"\n")
         return(invisible())
    }


    if(!is.genpop(x)) stop("x is not a valid genpop object")
    res.type <- match.arg(res.type)
    if(res.type=="genpop") { truenames <- TRUE }

    temp <- x@loc.fac
    nloc <- length(levels(temp))
    levels(temp) <- 1:nloc

    kX <- list()

    for(i in 1:nloc){
        kX[[i]] <- matrix(x@tab[,temp==i],ncol=x@loc.nall[i])

        if(!truenames){
            rownames(kX[[i]]) <- rownames(x@tab)
            colnames(kX[[i]]) <- paste(names(x@loc.names)[i],names(x@all.names[[i]]),sep=".")
        }else{
            rownames(kX[[i]]) <- x@pop.names
            colnames(kX[[i]]) <- paste(x@loc.names[i],x@all.names[[i]],sep=".")
        }
    }

    if(truenames) {
        names(kX) <- x@loc.names
    } else{
        names(kX) <- names(x@loc.names)
    }

    prevcall <- match.call()
    if(res.type=="genpop"){
        kX <- lapply(kX, genpop, prevcall=prevcall, ploidy=x@ploidy, type=x@type)
        for(i in 1:length(kX)){
            kX[[i]]@other <- x@other
        }
    }

    return(kX)
})





