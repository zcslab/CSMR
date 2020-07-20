#' Function to Fit Mixture of Regressions
#'
#' The main function in this package.
#'
#' @author The mixtureReg package is developed by Tianxia Zhou on github. All right reserved by Tianxia Zhou.
#' @param regData data frame used in fitting model.
#' @param formulaList a list of the regression components that need to be estimated.
#' @param xName character; Name used to pick x variable from data.
#' @param yName character; Name used to pick y variable from data.
#' @param mixingProb character;
#' Specify how the mixing probabilities are estimated in the M step.
#' "Constant" specifies a constant mixing probabilities;
#' "loess" specifies predictor dependent mixing probabilities obtained by loess smoothing.
#' @param initialWList a list of weights guesses (provided by user).
#' Typically this is not used, unless the user has a good initial guess.
#' @param epsilon a small value that the function consider as zero.
#' The value is used in determine matrix sigularity and in determine convergence.
#' @param max_iter the maximum number of iterations.
#' @param max_restart the maximum number of restart before giving up.
#' @param min_lambda a value used to ensure estimated mixing probabilities (lambda's) are not too close to zero.
#' @param min_sigmaRatio a value used to prevent estimated variaces of any regression component from collapsing to zero.
#' @param silently a switch to turn off the screen printout.
#' @return A class 'mixtureReg' object.
#'
#' @export mixtureReg
mixtureReg <- function(regData, formulaList,
                       xName = NULL, yName = NULL,
                       mixingProb = c("Constant", "loess"),
                       initialWList = NULL,
                       epsilon = 1e-08, max_iter = 10000, max_restart = 15,
                       min_lambda = 0.01, min_sigmaRatio = 0.1,
                       silently = TRUE
) {
  if(is.null(yName)) {yName = all.vars(formulaList[[1]])[1]}
  if(is.null(xName)) {xName = all.vars(formulaList[[1]])[2]}

  # missing values in data or data with no variance will crash the algorithm
  stopifnot(
    all(!is.na(regData[, yName])),
    var(regData[ ,yName]) > 0
  )

  #require(dplyr)

  # E step: update WList, ll
  conditionalP <- function(res, lambda, sigma) {
    # res: vector of residuals from predictions
    # lambda: vector or scalar of prior probability
    # sigma: variance estimate
    # return: vector of posterior probability
    lambda * dnorm(x = res, mean = 0, sd = sigma)
  }

  EUpdate <- function(result) {
    lmList = result$lmList
    lambdaList = result$lambdaList
    resList = lapply(lmList, function(m) m$residuals)
    sigmaList = lapply(X = lmList, FUN = function(m) summary(m)$sigma)
    PList = mapply(conditionalP, resList, lambdaList, sigmaList, SIMPLIFY = FALSE)
    sumsP = rowSums(do.call(cbind, PList))
    WList = lapply(PList, function(p) p/sumsP)
    llList = lapply(lmList, logLik_mixtureReg)
    ll = sum(log(sumsP)) # Log-likelihood
    return(list(WList = WList, ll = ll))
  }

  # M step: update lmList, lambdaList
  MUpdate <- function(WList) {
    updateLambda <- function(WList) {
      standardizeW <- function(WList) {
        WSums <- WList %>%
          unlist %>%
          (function(x) (x > min_lambda)*x + (x <= min_lambda)*min_lambda) %>%
          matrix(., ncol = length(WList)) %>%
          #matrix(.data, ncol = length(WList)) %>%
          rowSums()
        lapply(X = WList, FUN = function(x) x/WSums)
      }

      estimateLambda <- function(WList) {
        # function to estimate prior weights
        if (mixingProb == "Constant") {
          lamList <- lapply(X = WList,
                            FUN = function(x) rep(mean(x), length(x)))
        } else if (mixingProb == "loess") {
          lamList <- lapply(
            X = WList,
            FUN = function(w) {
              predict(loess(formula = w ~ regData[ , xName],
                            degree = 0,
                            control=loess.control(surface="direct")))
            }
          )
        }
        return(lamList)
      }

      WList %>%
        standardizeW %>%
        estimateLambda
    }

    lambdaList = updateLambda(WList = WList)

    wlm <- function(fml, W) {
      # fml: formula
      # W: weights
      tempData <- regData
      tempData$Wt <- W/mean(W)  # rescale weights so that mean == 1
      #print(head(tempData))
      #print(W)
      #myweight <<- W/mean(W) #global variable
      reg = lm(formula = fml, weights = Wt, data = tempData)
      #reg = lm(formula = fml, weights = myweight, data = tempData) #interesting!!! Adding global variable on package-level. see global_variable.R for example.
      #print(reg$weights)
      return(reg)
    }

    lmList = mapply(FUN = wlm, formulaList, WList, SIMPLIFY = FALSE)

    return(list(lmList = lmList, lambdaList = lambdaList))
  }

  # restart step
  isSingular <- function(lambdaList) {
    return(any(unlist(lambdaList) < epsilon))
  }

  needRestart <- function(newResult, newLL, iter, restart) {
    if (newLL < ll) {
      errorMessage <- "smaller logLik"
      answer <- TRUE
    }

    if (isSingular(newResult$lambdaList)) {
      errorMessage <- "sigular lambda list"
      answer <- TRUE
    } else {
      if (is.na(newLL) || abs(newLL) == Inf) {
        errorMessage <- "abnormal logLik"
        answer <- TRUE
      } else {
        sigma_s = sapply(X = newResult$lmList, FUN = function(m) summary(m)$sigma)
        if (any(is.na(sigma_s))) {
          errorMessage <- "sigma is NA"
          answer <- TRUE
        } else {
          sigmaRatio_s = sigma_s/sigma_s[1]
          if (any(sigmaRatio_s < min_sigmaRatio) ||
              any(sigmaRatio_s > 1/min_sigmaRatio)) {
            errorMessage <- "sigma ratio constraint"
            answer <- TRUE
          } else {
            errorMessage <- NA
            answer <- FALSE
          }
        }
      }
    }
    return(list("error_message" = errorMessage, "answer" = answer))
  }

  randomWList <- function(n, k) {
    # n: number of obs
    # k: number of regression lines
    randomWeightMatrix = matrix(runif(n*k, min = 0.1, max = 0.9), nrow = n, ncol = k)
    WMatrix = randomWeightMatrix / (rowSums(randomWeightMatrix))
    WList = lapply(as.list(1:k), function(i) WMatrix[,i])
    return(WList)
  }

  # Main routine
  {
    # initialize
    {
      diff <- 1
      iter <- 0
      restart <- 0

      n = dim(regData)[1]
      k = length(formulaList)

      # initialize (E step)
      if (!is.null(initialWList)) { # get WList from initial guess
        WList = initialWList
      } else {
        WList = randomWList(n, k)
      }

      # first M step
      newResult = MUpdate(WList)
      # second E step
      newE = EUpdate(newResult)
      # update
      result <- newResult
      WList <- newE$WList
      ll <- newE$ll

      monitor <- data.frame(
        "diff" = diff, "iter" = iter, "restart" = restart, "logLik" = ll,
        "newLL" = NA,
        "sigma1" = summary(newResult$lmList[[1]])$sigma,
        "sigma2" = summary(newResult$lmList[[2]])$sigma,
        "ratio" = summary(newResult$lmList[[1]])$sigma/summary(newResult$lmList[[2]])$sigma,
        "lambda1" = mean(newResult$lambdaList[[1]]),
        "lambda2" = mean(newResult$lambdaList[[2]]),
        "error_message" = NA)
    }

    # while loop
    while (abs(diff) > epsilon && iter < max_iter && restart < max_restart) { # while not convegent
      # do an M step
      # do an E step
      # if sigular or no improvement, restart;
      # else update result and ll, and do another iteration.
      newResult = MUpdate(WList)
      newE = EUpdate(newResult)
      newLL = newE$ll

      judge <- needRestart(newResult, newLL, iter, restart)
      if (judge$"answer") {
        WList = randomWList(n, k)
        restart <- restart + 1
      } else {
        diff = newLL - ll
        result <- newResult
        WList <- newE$WList
        ll <- newLL
      }
      iter <- iter + 1
      monitor <-
        rbind(
          monitor,
          c(diff, iter, restart, ll,
            newLL,
            summary(newResult$lmList[[1]])$sigma,
            summary(newResult$lmList[[2]])$sigma,
            summary(newResult$lmList[[1]])$sigma/summary(newResult$lmList[[2]])$sigma,
            mean(newResult$lambdaList[[1]]),
            mean(newResult$lambdaList[[2]]),
            judge$"error_message")
        )
    }

    if (!silently) {
      cat("diff = ", diff, "\n")
      cat("iter = ", iter, "\n")
      cat("restart = ", restart, "\n")
      cat("log-likelihood = ", ll, "\n")
    }

    mixtureRegModel = list(
      "lmList" = result$lmList,
      "logLik" = ll,
      "monitor" = monitor,
      "regData" = regData,
      "prior" = result$lambdaList,
      "posterior" = WList)
    class(mixtureRegModel) = c("mixtureReg", class(mixtureRegModel))
    return(mixtureRegModel)
  }
}

#------

# methods for mixtureReg

#' Sort by X Coordinates and Add Line to a Plot
#'
#' Rearrange X and Y coordinates before calling "lines()" function.
#'
#' @param x X coordinate vectors of points to join.
#' @param y Y coordinate vectors of points to join.
#' @param ...	Further graphical parameters.
orderedLines <- function(x, y, ...) {
  # a helper function used in plotting
  xOrder <- order(x)
  lines(x = x[xOrder], y = y[xOrder], ...)
}


#' Plot Fit and Mixing Probability of a mixtureReg Object
#'
#' S3 plot method for class 'mixtureReg'.
#'
#' @param mixtureModel mixtureReg object, typically result from 'mixtureReg()'.
#' @param which numeric; choose which plot to display.
#' '1' gives a plot of fit; '2' gives a plot of mixing probability.
#' @param xName character; Name used to pick x variable from data.
#' @param yName character; Name used to pick y variable from data.
#' @param xlab character; label that should be put on the x axis.
#' @param ylab character; label that should be put on the y axis.
#' @param ...	Further graphical parameters.
#'
#' @export
#@S3method plot mixtureReg
plot_mixtureReg <- function(mixtureModel, which = 1:2,
                            xName = NULL, yName = NULL,
                            xlab = NULL, ylab = NULL,
                            ...) {
  # plot method for "mixtureReg" class

  if (is.null(yName)) {yName = all.vars(mixtureModel$lmList[[1]]$terms)[1]}
  if (is.null(xName)) {xName = all.vars(mixtureModel$lmList[[1]]$terms)[2]}
  if (is.null(xlab)) {xlab = xName}
  if (is.null(ylab)) {ylab = yName}

  XX = mixtureModel$regData[ , xName]
  YY = mixtureModel$regData[ , yName]
  YhatList = lapply(X = mixtureModel$lmList, FUN = function(x) predict(x))

  if (which == 1) {
    pt_color = rep('grey60', length(XX))
    for (i in 1:length(mixtureModel$lmList)) {
      ind_loca = which(mixtureModel$posterior[[i]] > mixtureModel$prior[[i]])
      pt_color[ind_loca] = c("tomato","dodgerblue","gold","greenyellow","darkorchid1","darkorange")[i]
    }
    plot(x = XX, y = YY, xlab = xlab, ylab = ylab, pch=19, col = alpha(pt_color, 0.9),main="Mixture Regression", ...)
    for (i in 1:length(mixtureModel$lmList)) {
      orderedLines(x = XX, y = YhatList[[i]], col = c("tomato","dodgerblue","gold","greenyellow","darkorchid1","darkorange")[i])
    }
  }

  if (which == 2) {
    for (i in 1:length(mixtureModel$lmList)) {
      plot(x = XX, y = mixtureModel$posterior[[i]],
           xlab = xlab, ylab = paste0("Weights_", i),
           ylim = c(-0.01,1.01),
           pch=19, col = 'grey', main='Posterior',
           ...)
      orderedLines(x = XX, y = mixtureModel$prior[[i]], col = c("tomato","dodgerblue","gold","greenyellow","darkorchid1","darkorange")[i])
    }
  }
}

#' Plot a List of mixtureReg Objects
#'
#' Feed in a list of mixtureReg models and get an overlayed plot.
#'
#' @param mixtureRegList a list of multiple mixtureReg objects.
#' @param xName character; Name used to pick x variable from data.
#' @param yName character; Name used to pick y variable from data.
#' @param ...	Further graphical parameters.
#'
#' @export
#@export plot.mixtureRegList
plot_mixtureRegList <- function(mixtureRegList,
                                xName = NULL, yName = NULL,
                                ...) {
  # plot overlayed plots for a list of 'mixtureReg' models

  getPlotData <- function(mReg) {
    XX = mReg$regData[ , xName]
    YY = mReg$regData[ , yName]
    YhatList = lapply(X = mReg$lmList, FUN = function(x) predict(x))
    return(list("XX" = XX,
                "YY" = YY,
                "YhatList" = YhatList))
  }

  plotDataList <- lapply(X = mixtureRegList,
                         FUN = getPlotData)

  plot(x = bind_rows(
    lapply(X = plotDataList,
           FUN = function(pd) {
             dd <- data_frame("Xs" = pd$"XX", "Ys" = pd$"YY")
             return(dd)
           }
    )
  ),
  xlab = xName, ylab = yName, type = 'n',
  ...)

  for (i in seq_along(plotDataList)) {
    points(x = plotDataList[[i]]$"XX",
           y = plotDataList[[i]]$"YY",
           col = i + 1, pch = i + 1)
    for (j in seq_along(plotDataList[[i]]$"YhatList")) {
      orderedLines(x = plotDataList[[i]]$"XX",
                   y = plotDataList[[i]]$"YhatList"[[j]],
                   col = i + 1)
    }
  }
}

#' Obtain Log-likelihood from a mixtureReg Object
#'
#' S3 method for class 'mixtureReg'.
#' However, it doesn't return a 'logLik' object.
#' For simlicity, it returns a 'numeric' value.
#'
#' @param mixtureModel mixtureReg object, typically result from 'mixtureReg()'.
#' @return Return a numeric value of log likelihood.
#'
#' @export
#@S3method logLik mixtureReg
logLik_mixtureReg <- function(mixtureModel) {
  return(mixtureModel$"logLik")
}


#------------------example by wnchang@iu.edu----------------
# #library(mixtools)
# data("CO2data")
# head(CO2data)
#
# #library(mixtureReg)
# mx1 <- mixtureReg(
#   regData = CO2data,
#   formulaList = list(formula(CO2 ~ GNP),
#                      formula(CO2 ~ GNP)),
#   mixingProb = "Constant"
# )
#
# plot_mixtureReg(mx1, which = 1)
#
# plot_mixtureReg(mx1, which = 2)
