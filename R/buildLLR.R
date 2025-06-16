# Copyright (C) 2020  Jochen Weile, Roth Lab
# Edited by Cindy Zhang, Roth Lab (2024)
#
# This file is part of maveLLR
#
# maveLLR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# maveLLR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with maveLLR  If not, see <https://www.gnu.org/licenses/>.

#' Build LLR function using kernel density estimation
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param bw the bandwith to use for kernel density estimation (see kdensity package for more info).
#'   This can either be a numerical value or the name of the algorithm used to automatically choose one.
#' @param kernel the type of kernel to use (see kdensity package for more info)
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.kernel(posScores,negScores,bw=0.1,kernel="gaussian")
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.kernel <- function(posScores, negScores, bw=0.1, kernel="gaussian", outlierSuppression=0.0001) {
  
  library(kdensity)
  
  posDens <- kdensity(posScores,bw=bw,kernel=kernel)
  negDens <- kdensity(negScores,bw=bw,kernel=kernel)
  
  #generate a dummy outlier point far away from the rest of the distribution
  # so we can use it to measure it's hypothetical density
  pseudo <- min(posScores) - 10*ifelse(is.numeric(bw),bw,0.1)
  #measure the hypothetical density of an outlier point
  minDensPos <- kdensity(c(posScores[-1],pseudo),bw=bw,kernel=kernel)(pseudo)
  #do the same for the negative distribution
  pseudo <- min(negScores) - 10*ifelse(is.numeric(bw),bw,0.1)
  minDensNeg <- kdensity(c(negScores[-1],pseudo),bw=bw,kernel=kernel)(pseudo)
  #whichever outlier density is higher will serve as our uniform prior
  minDens <- outlierSuppression * max(minDensPos,minDensNeg)
  
  llrFun <- function(score) sapply(score,function(s)
    # log10( max(posDens(s),minDens) / max(negDens(s),minDens) )
    log10(max(posDens(s),minDens)) - log10(max(negDens(s),minDens))
  )
  
  return(list(llr=llrFun,posDens=posDens,negDens=negDens))
}

#' Build LLR function using kernel density estimation - experimental version
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param bw the bandwith to use for kernel density estimation (see kdensity package for more info).
#'   This can either be a numerical value or the name of the algorithm used to automatically choose one.
#' @param kernel the type of kernel to use (see kdensity package for more info)
#' 
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.kernel(posScores,negScores,bw=0.1,kernel="gaussian")
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.kernelExperimental <- function(posScores, negScores, bw=0.1, kernel="gaussian", outlierSuppression=0.0001) {
  
  library(kdensity)
  
  posDens <- kdensity(posScores,bw=bw,kernel=kernel)
  negDens <- kdensity(negScores,bw=bw,kernel=kernel)
  refDens <- kdensity(c(posScores,negScores), bw=bw, kernel=kernel)
  
  # generate a dummy outlier point far away from the rest of the distribution
  # so we can use it to measure it's hypothetical density
  pseudo <- min(c(posScores, negScores)) - 10*ifelse(is.numeric(bw),bw,0.1)
  # measure the hypothetical density of an outlier point
  pseudoDens <- kdensity(c(posScores[-1],negScores,pseudo),bw=bw,kernel=kernel)(pseudo)
  priorWeight = pseudoDens * outlierSuppression
  
  llrFun <- function(score) sapply(score,function(s) {
    rawLLR = log10(posDens(s)) - log10(negDens(s))
    weight = refDens(s) / (refDens(s)+priorWeight)
    finalLLR = rawLLR * weight
  })
  
  return(list(llr=llrFun,posDens=posDens,negDens=negDens))
}


#' Build LLR function using Gaussian densities.
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param spline logical; whether or not to apply spline monotonization. TRUE by default.
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#'
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.gauss(posScores,negScores)
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.gauss <- function(posScores, negScores, spline=TRUE) {
  
  mpos <- mean(posScores)
  spos <- sd(posScores)
  mneg <- mean(negScores)
  sneg <- sd(negScores)
  
  posDens <- function(x) dnorm(x,mpos,spos)
  negDens <- function(x) dnorm(x,mneg,sneg)
  
  llrFun <- function(score) log10(posDens(score)/negDens(score))
  
  if (spline) {
    minPoint <- optimize(llrFun,interval=c(mpos, qnorm(0.999,mneg,sneg)),maximum=FALSE)
    if (minPoint$minimum < mneg) {
      minPoint$minimum <- Inf
    }
    maxPoint <- optimize(llrFun,interval=c(qnorm(0.001,mpos,spos),mneg),maximum=TRUE)
    if (maxPoint$maximum > mpos) {
      maxPoint$maximum <- -Inf
    }
    
    llrSpline <- function(scores) {
      sapply(scores,function(score) {
        if (is.na(score)) {
          NA
        } else if (score > minPoint$minimum) {
          minPoint$objective
        } else if (score < maxPoint$maximum) {
          maxPoint$objective
        } else {
          llrFun(score)
        }
      })
    }
    
    return(list(llr=llrSpline,posDens=posDens,negDens=negDens))
    
  } else {
    return(list(llr=llrFun,posDens=posDens,negDens=negDens))
  }
  
}


#' Draw a plot that summarizes the densities and LLR calculated by other functions
#'
#' @param scores the numerical scores in the Mave map
#' @param llrFun the LLR function
#' @param posDens the density function of the positive reference set
#' @param negDens the density function of the negative reference set
#' @param posScores the numerical scores in the postive reference set
#' @param negScores the numerical scores in the negative reference set
#'
#' @return nothing
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.gauss(posScores,negScores)
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
drawDensityLLR <- function(scores, llrFun, posDens, negDens, posScores, negScores, prior=0.1) {
  
  llrTs <- llrThresholds(optiLLR(prior))
  
  opar <- par(mfrow=c(2,1))
  xlim <- range(scores,na.rm=TRUE,finite=TRUE)
  ymax <- max(c(
    posDens(seq(xlim[[1]],xlim[[2]],length.out = 100)),
    negDens(seq(xlim[[1]],xlim[[2]],length.out = 100))
  ))
  par(mar=c(.1,4,1,1))
  xs <- seq(xlim[[1]],xlim[[2]],length.out=200)
  ys <- llrFun(xs)
  ylim <- range(ys,na.rm=TRUE,finite=TRUE)
  plot(NA,type="n",xlim=xlim,ylim=ylim,
       axes=FALSE,xlab="",ylab="LLR"
  )
  drawThresh <- function(t,col,label) {
    if (t >= 0) {
      if (t < ylim[[2]]) {
        rect(xlim[[1]],t,xlim[[2]],ylim[[2]],col=col,border=NA)
        text(.1*xlim[[2]]+.9*xlim[[1]],t,label,pos=3,cex=0.8)
      }
    } else {
      if (t > ylim[[1]]) {
        rect(xlim[[1]],ylim[[1]],xlim[[2]],t,col=col,border=NA)
        text(.1*xlim[[2]]+.9*xlim[[1]],t,label,pos=1,cex=0.8)
      }
    }
  }
  drawThresh(llrTs[["patho.support"]],"lemonchiffon","Patho. support.")
  drawThresh(llrTs[["patho.moderate"]],"lightgoldenrod1","Patho. moderate")
  drawThresh(llrTs[["patho.strong"]],"goldenrod1","Patho. strong")
  drawThresh(llrTs[["patho.vstrong"]],"indianred1","Patho. very str.")
  drawThresh(llrTs[["benign.support"]],"lemonchiffon","Benign support.")
  drawThresh(llrTs[["benign.strong"]],"goldenrod1","Benign strong")
  
  lines(xs,ys)
  # plot(llrFun,from=xlim[[1]],to=xlim[[2]],xlim=xlim,axes=FALSE,xlab="",ylab="LLR")
  abline(h=0,col="gray",lty="dashed")
  axis(2)
  par(mar=c(5,4,.1,1))
  hist(scores,col="gray90",border=NA,freq=FALSE,main="",xlim=xlim,ylim=c(0,ymax),breaks=50)
  plot.function(posDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="firebrick3",lwd=2)
  plot.function(negDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="darkolivegreen3",lwd=2)
  abline(v=posScores,col="firebrick3")
  abline(v=negScores,col="darkolivegreen3")
  par(opar)
  
  return(invisible(NULL))
}

#' @export
optiLLR <- function(prior=0.1,posterior=0.9) {
  log10(posterior*(1-prior)/(prior*(1-posterior)))*4/3
}

#' @export
llrThresholds <- function(LLRpvst=optiLLR(0.1),X=2) {
  c(
    patho.vstrong=LLRpvst,
    patho.strong=LLRpvst/X,
    patho.moderate=LLRpvst/(X^2),
    patho.support=LLRpvst/(X^3),
    benign.support=-LLRpvst/(X^3),
    benign.strong=-LLRpvst/(X^1)
  )
}

#' Like drawDensityLLR, but with fixed y range for better comparison across different LLRs
#'
#' @param scores the numerical scores in the Mave map
#' @param llrFun the LLR function
#' @param posDens the density function of the positive reference set
#' @param negDens the density function of the negative reference set
#' @param posScores the numerical scores in the postive reference set
#' @param negScores the numerical scores in the negative reference set
#'
#' @return nothing
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.gauss(posScores,negScores)
#' drawDensityLLR_fixedRange(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
drawDensityLLR_fixedRange <- function(scores, llrFun, posDens, negDens, posScores, negScores, prior=0.1) {
  
  llrTs <- llrThresholds(optiLLR(prior))
  
  opar <- par(mfrow=c(2,1))
  xlim <- range(scores,na.rm=TRUE,finite=TRUE)
  ymax <- max(c(
    posDens(seq(xlim[[1]],xlim[[2]],length.out = 100)),
    negDens(seq(xlim[[1]],xlim[[2]],length.out = 100))
  ))
  par(mar=c(.1,4,1,1))
  xs <- seq(xlim[[1]],xlim[[2]],length.out=200)
  ys <- llrFun(xs)
  #ylim <- range(ys,na.rm=TRUE,finite=TRUE)
  ylim <- c(-3,4) # fixed range
  plot(NA,type="n",xlim=xlim,ylim=ylim,
       axes=FALSE,xlab="",ylab="LLR"
  )
  drawThresh <- function(t,col,label) {
    if (t >= 0) {
      if (t < ylim[[2]]) {
        rect(xlim[[1]],t,xlim[[2]],ylim[[2]],col=col,border=NA)
        text(.1*xlim[[2]]+.9*xlim[[1]],t,label,pos=3,cex=0.8)
      }
    } else {
      if (t > ylim[[1]]) {
        rect(xlim[[1]],ylim[[1]],xlim[[2]],t,col=col,border=NA)
        text(.1*xlim[[2]]+.9*xlim[[1]],t,label,pos=1,cex=0.8)
      }
    }
  }
  drawThresh(llrTs[["patho.support"]],"#FFC0CB","Patho. support.")
  drawThresh(llrTs[["patho.moderate"]],"#FF9999","Patho. moderate")
  drawThresh(llrTs[["patho.strong"]],"#FF6666","Patho. strong")
  drawThresh(llrTs[["patho.vstrong"]],"red","Patho. very str.")
  drawThresh(llrTs[["benign.support"]],"lightblue","Benign support.")
  drawThresh(llrTs[["benign.strong"]],"dodgerblue","Benign strong")
  
  lines(xs,ys)
  # plot(llrFun,from=xlim[[1]],to=xlim[[2]],xlim=xlim,axes=FALSE,xlab="",ylab="LLR")
  abline(h=0,col="gray",lty="dashed")
  axis(2)
  par(mar=c(5,4,.1,1))
  hist(scores,col="gray90",border=NA,freq=FALSE,main="",xlim=xlim,ylim=c(0,ymax),breaks=50)
  plot.function(posDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="firebrick3",lwd=2)
  plot.function(negDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="deepskyblue2",lwd=2)
  abline(v=posScores,col="firebrick3")
  abline(v=negScores,col="deepskyblue2")
  par(opar)
  
  return(invisible(NULL))
}

#' Find x-values at which the given LLR function crosses specified thresholds.
#'
#' @param llrFun the LLR Function
#' @param thresholds thresholds for each ACMG category
#' @param xlim range of x values i.e. scores
#' @param nPoints an integer specifying number of points to sample across xlim
#'
#' @return a table with threshold labels, values, and x when crossing
#' @export
#'
#' @examples
#' # Using a previously defined LLR function
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llrList <- buildLLR.kernel(posScores, negScores, bw=0.1, kernel="gaussian")
#' llrFun <- llrList$llr
#'
#' # Define thresholds of interest
#' llrTs <- llrThresholds(optiLLR(0.1))
#'
#' # Find the crossings in the range of observed scores
#' x_range <- range(c(posScores, negScores))
#' findLLRcrossings(llrFun, llrTs, xlim = x_range)
findLLRcrossings <- function(llrFun, thresholds, xlim, nPoints = 1000) {
  xs <- seq(xlim[1], xlim[2], length.out = nPoints)
  ys <- llrFun(xs)
  
  # Get indices
  patho_idx <- which(names(thresholds) == "patho.support")
  benign_support_idx <- which(names(thresholds) == "benign.support")
  benign_strong_val <- thresholds[which(names(thresholds) == "benign.strong")]
  
  # Insert 'none' as category between patho.support and benign.support
  benign_support_val <- thresholds[benign_support_idx]
  thresholds <- append(
    thresholds,
    setNames(benign_support_val, "none"),
    after = patho_idx
  )
  
  # Assign the lower value of none as the upper value of benign.support, but without changing display thresholds
  display_thresholds <- thresholds
  thresholds[which(names(thresholds) == "benign.support")] <- benign_strong_val
  display_thresholds["benign.support"] <- thresholds["none"]
  
  thresholdNames <- names(thresholds)
  
  # Interval range formatting for thresholds
  format_range <- function(name, i) {
    if (name == "none") return("[-0.32, 0.32)")
    val <- display_thresholds[i]
    if (is.na(val)) return("")
    
    if (grepl("^patho\\.", name)) {
      upper <- if (i > 1 && !is.na(display_thresholds[i - 1])) display_thresholds[i - 1] else Inf
      return(sprintf("[%.2f, %s)", val, ifelse(is.infinite(upper), "∞", sprintf("%.2f", upper))))
    } else if (grepl("^benign\\.", name)) {
      lower <- if (i < length(display_thresholds) && !is.na(display_thresholds[i + 1])) display_thresholds[i + 1] else -Inf
      bracket <- if (is.infinite(lower)) "(" else "["
      return(sprintf("%s%s, %.2f)", bracket, ifelse(is.infinite(lower), "-∞", sprintf("%.2f", lower)), val))
    }
    return("")
  }
  
  formatted_thresholds <- setNames(sapply(seq_along(thresholds), function(i) {
    format_range(names(thresholds)[i], i)
  }), names(thresholds))
  
  # Determine crossings and direction of crossing
  crossings <- list()
  for (i in seq_along(thresholds)) {
    thr <- thresholds[i]
    
    fvals <- ys - thr
    sign_changes <- which(fvals[-length(fvals)] * fvals[-1] < 0)
    
    for (j in sign_changes) {
      interval <- c(xs[j], xs[j + 1])
      if (!is.finite(fvals[j]) || !is.finite(fvals[j + 1])) next
      
      root <- tryCatch({
        uniroot(function(z) llrFun(z) - thr, interval = interval)$root
      }, error = function(e) NA)
      
      if (!is.na(root)) {
        slope <- ys[j + 1] - ys[j]
        cur_idx <- i
        
        if (slope > 0 && cur_idx + 1 <= length(thresholds)) {
          to <- thresholdNames[cur_idx]
          from <- thresholdNames[cur_idx + 1]
        } else if (slope < 0 && cur_idx + 1 <= length(thresholds)) {
          from <- thresholdNames[cur_idx]
          to <- thresholdNames[cur_idx + 1]
        } else {
          next
        }
        
        crossings[[length(crossings) + 1]] <- list(
          from = from,
          to = to,
          score = round(root, 2)
        )
      }
    }
  }
  
  # Build crossing data frame
  crossing_df <- do.call(rbind, lapply(crossings, as.data.frame))
  crossing_df$score <- as.numeric(as.character(crossing_df$score))
  
  # Point at where LLR = 0
  fvals <- ys
  zero_crossings <- which(fvals[-length(fvals)] * fvals[-1] < 0)
  for (j in zero_crossings) {
    interval <- c(xs[j], xs[j + 1])
    root <- tryCatch({
      uniroot(function(z) llrFun(z), interval = interval)$root
    }, error = function(e) NA)
    if (!is.na(root)) {
      crossing_df <- rbind(
        crossing_df,
        data.frame(
          from = "LLR = 0",
          to = "",
          score = round(root, 2),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  if (nrow(crossing_df) == 0) return(NULL)  # Return early if no crossings at all
  
  crossing_df <- crossing_df[order(crossing_df$score), ]
  
  # Initialize output with first 'from' row
  output <- list()
  
  first_from <- crossing_df$from[1]
  if (!is.na(first_from) && first_from != "LLR = 0") {
    output[[length(output) + 1]] <- data.frame(
      LLR_category = first_from,
      LLR_threshold = formatted_thresholds[[first_from]],
      score_crossing = ""
    )
  }
  
  # Add crossing rows and 'to' category rows
  for (i in seq_len(nrow(crossing_df))) {
    from <- crossing_df$from[i]
    to <- crossing_df$to[i]
    score <- crossing_df$score[i]
    
    if (from == "LLR = 0") {
      # Insert a blank row before "Crossing 0"
      #output[[length(output) + 1]] <- data.frame(
      # LLR_category = "",
      # LLR_threshold = "",
      #  score_crossing = ""
      # )
      
      output[[length(output) + 1]] <- data.frame(
        LLR_category = "LLR = 0",
        LLR_threshold = 0,
        score_crossing = score
      )
    } else {
      output[[length(output) + 1]] <- data.frame(
        LLR_category = "---------------------",
        LLR_threshold = "---------------------",
        score_crossing = score
      )
      
      output[[length(output) + 1]] <- data.frame(
        LLR_category = to,
        LLR_threshold = formatted_thresholds[[to]],
        score_crossing = ""
      )
    }
  }
  
  do.call(rbind, output)
}


#
# mthfr <- read.csv("~/projects/mthfr/folate_response_model5.csv")
# variants <- mthfr[mthfr$type=="substitution","hgvs"]
# scores <- mthfr[mthfr$type=="substitution","m25.score"]
#
# idb <- 337
# prsnrs <- read.csv("~/projects/mthfr/mthfr_prsNrs.csv")
#
# posref <- prsnrs[prsnrs$reference=="positive" & prsnrs$start < idb,"hgvs"]
# negref <- prsnrs[prsnrs$reference=="negative" & prsnrs$start < idb,"hgvs"]
# posScores <- na.omit(scores[variants %in% posref])
# negScores <- na.omit(scores[variants %in% negref])
# llr <- buildLLR.kernel(posScores, negScores,bw=0.1,kernel="gaussian")
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
# llr <- buildLLR.gauss(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#
# posref <- prsnrs[prsnrs$reference=="positive" & prsnrs$start > idb,"hgvs"]
# negref <- prsnrs[prsnrs$reference=="negative" & prsnrs$start > idb,"hgvs"]
# posScores <- na.omit(scores[variants %in% posref])
# negScores <- na.omit(scores[variants %in% negref])
# llr <- buildLLR.kernel(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
# llr <- buildLLR.gauss(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
