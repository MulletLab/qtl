######################################################################
#
# est.rf.exHet.R
#
# copyright (c) 2001-2013, Karl W Broman; modified 2014 by Ryan McCormick
# last modified April, 2014
# first written Apr, 2001
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
# 
# Part of the R/qtl package
# Contains: est.rf.exHet, mapfromRF
#
######################################################################

######################################################################
#
# est.rf.exHet: Estimate sex-averaged recombination fractions between
#         all pairs of markers
#
# mapFromRF: Use pairwise recombination fractions of adjacent markers
#	to calculate a genetic map
#
######################################################################

est.rf.exHet <-
function(cross, maxit=10000, tol=1e-6, het=0.5) 
{
  print("This code is in development by RFM and SKT, and should not be used for general purposes")

  if(!any(class(cross) == "cross"))
    stop("Input should have class \"cross\".")

  n.chr <- nchr(cross)
  n.mar <- totmar(cross)
  n.ind <- nind(cross)
  mar.names <- unlist(lapply(cross$geno,function(a) colnames(a$data)))
  
  type <- class(cross)[1]
  chrtype <- sapply(cross$geno,class)

  is.bcsft <- (type == "bcsft")
  if(is.bcsft) {
    cross.scheme <- attr(cross, "scheme") ## c(s,t) for BC(s)F(t)
    is.bcsft <- cross.scheme[2] > 0 ## used for fixX only
  }
  
  xchrcol <- NULL
  fixX <- FALSE
  Geno <- NULL
  # create full genotype matrix
  for(i in 1:n.chr) {
    temp <- cross$geno[[i]]$data

    # treat X chromosome specially in an intercross or BCsFt with t>0.
    if((type=="f2" || is.bcsft) && chrtype[i]=="X") {
      fixX <- TRUE
      if(i != 1) xchrcol <- c(xchrcol,ncol(Geno)+(1:ncol(cross$geno[[i]]$data)))
      else xchrcol <- 1:ncol(cross$geno[[i]]$data)
      xchr <- temp
      xchr[is.na(xchr)] <- 0
      temp <- reviseXdata("f2","simple",getsex(cross),geno=temp,
                          cross.attr=attributes(cross))
    }
    Geno <- cbind(Geno,temp)
  }

  # which type of cross is this?
  if(type == "f2")
    cfunc <- "est_rf_f2"
  else if(type == "bc" || type=="risib" || type=="riself" || type=="dh" || type=="haploid") 
    cfunc <- "est_rf_bc"
  else if(type == "4way") 
    cfunc <- "est_rf_4way"
  else if(type=="ri8sib" || type=="ri8self" || type=="ri4sib" || type=="ri4self") {
    cfunc <- paste("est_rf_", type, sep="")
    if(any(chrtype == "X"))
      warning("est.rf not working properly for the X chromosome for 4- or 8-way RIL.")
  }
  else if(type == "bcsft")
    cfunc <- "est_rf_bcsft_exHet"
  else 
    stop("est.rf not available for cross type ", type, ".")

  Geno[is.na(Geno)] <- 0
  
  if(type=="bc" || type=="risib" || type=="riself" || type=="dh" || type=="haploid")
    z <- .C(cfunc,
            as.integer(n.ind),         # number of individuals
            as.integer(n.mar),         # number of markers
            as.integer(Geno),
            rf = as.double(rep(0,n.mar*n.mar)),
            PACKAGE="qtl")
  else {
    ## Hide cross scheme in genoprob to pass to routine. BY
    temp <- as.double(rep(0,n.mar*n.mar))
    if(type == "bcsft")
      temp[1:2] <- cross.scheme
    
    z <- .C(cfunc,
            as.integer(n.ind),         # number of individuals
            as.integer(n.mar),         # number of markers
            as.integer(Geno),
            rf = as.double(temp),
            as.integer(maxit),
            as.double(tol),
	    as.double(het),
            PACKAGE="qtl")
  }

  cross$rf <- matrix(z$rf,ncol=n.mar)
  dimnames(cross$rf) <- list(mar.names,mar.names)

  if(fixX) {
    temp <- as.double(rep(0, ncol(xchr) ^ 2))
    if(type == "bcsft") {
      temp[1] <- cross.scheme[1] + cross.scheme[2] - (cross.scheme[1] == 0)
    
      zz <- .C(cfunc,
               as.integer(n.ind),         # number of individuals
               as.integer(ncol(xchr)),   # number of markers on X chr.
               as.integer(xchr),
               rf = as.double(temp),
               as.integer(maxit),
               as.double(tol),
               PACKAGE="qtl")
    }
    else {
      zz <- .C("est_rf_bc",
               as.integer(n.ind),
               as.integer(ncol(xchr)),
               as.integer(xchr),
               rf=as.double(temp),
               PACKAGE="qtl")
    }
    zz <- matrix(zz$rf,ncol=ncol(xchr))
    cross$rf[xchrcol,xchrcol] <- zz
  }

  # check for alleles switches
  if(type == "risib" || type=="riself" || type=="f2" || type=="bc" || type=="dh" || type=="haploid") {
    out <- checkAlleles(cross, 5, FALSE)
    if(!is.null(out)) {
      out <- as.character(out[,1])
      warning("Alleles potentially switched at markers \n  ",
              paste(out, collapse=" "))
    }
  }

  cross
}

mapFromRF <- function(cross_cross, mapfunc=c("haldane", "kosambi", "c-f", "morgan"))
{
  i_numInd <- nind(cross_cross)
  i_numMar <- nmar(cross_cross)
  i_numChr <- nchr(cross_cross)
  newMap <- vector("list", i_numChr)
  names(newMap) <- names(cross_cross$geno)
  
  map.function <- match.arg(mapfunc)
  if(map.function=="kosambi") {
    mf <- mf.k; imf <- imf.k
  }
  else if(map.function=="c-f") {
    mf <- mf.cf; imf <- imf.cf
  }
  else if(map.function=="morgan") {
    mf <- mf.m; imf <- imf.m
  }
  else {
    mf <- mf.h; imf <- imf.h
  }
  
  for(j in 1:i_numChr) {
    
    rf_matrix <- pull.rf(cross_cross, chr=j)
    
    vector_rf <- vector(length=0)
    for (i in 2:(nrow(rf_matrix))) {
      vector_rf <- append(vector_rf, rf_matrix[i, i-1], after=length(vector_rf))
    }
    newMap[[j]] <- cumsum(c(0, imf(vector_rf)))
  }
  class(newMap) <- "map"
  newMap
}

# end of est.rf.exHet.R
