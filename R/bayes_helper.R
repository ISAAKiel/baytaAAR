#' @title Diagnostic summary of MCMC samples
#'
#' @description
#' Summarising diagnostics partly derived from
#'  \href{https://jkkweb.sitehost.iu.edu/DoingBayesianDataAnalysis/}{Kruschke (2015)}.
#'  Enter a codaMCMClist and receive a dataframe.
#'  {Columns:} PSRF Point est.,
#'  PSRF Upper C.I., Mean, Median, Mode, ESS, MCSE, HDImass, HDIlow, HDIhigh.
#'
#'  In certain situations, invoking the Gelman diagnostics can lead to errors,
#'  for example when constants are included or if the data is too sparse. In
#'  such cases, set first \code{gelman_diag_multivariate} to \code{FALSE}, and,
#'  if this still produces an error, deactivate \code{gelman_diag} altogether
#'  by setting it to \code{FALSE}, too.
#'
#'
#' @param codaMCMClist List. List in codaMCMClist format.
#' @param HDImass numeric. Value within 0 and 1. Default = 0.95.
#' @param gelman_diag logical. If TRUE, the gelman diagnostics for computing
#' the PSRF is invoked. Default: TRUE
#' @param gelman_diag_multivariate logical. If TRUE, the gelman diagnostics for
#' computing the PSRF multivariately is invoked. Default: FALSE.
#'
#' @rdname diagnostic.summary
#' @export
diagnostic.summary <- function(
    codaMCMClist,
    HDImass = 0.95,
    gelman_diag = TRUE,
    gelman_diag_multivariate = TRUE
    ) {
  parameterNames = varnames(codaMCMClist)
  mcmcMat = as.matrix(codaMCMClist,chains=TRUE)
  summaryInfo = NULL
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName], credMass = HDImass ) )
    thisRowName = parName
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  summaryInfo_df <- as.data.frame(summaryInfo)
  if(gelman_diag == TRUE) {
    psrf_df <- as.data.frame((gelman.diag(codaMCMClist, multivariate = gelman_diag_multivariate))$psrf)
    colnames(psrf_df) <- c("PSRF Point est.", "PSRF Upper C.I.")
    diagnostic_summary <- cbind(psrf_df, summaryInfo_df)
  }  else {
    diagnostic_summary <- summaryInfo_df
  }
}

# simplified version of a similar function in Kruschke 2015
#' @rdname summarizePost
#' @export
#' @noRd
summarizePost = function(
    paramSampleVec,
    credMass=0.95
    ) {
  paramSampleVec <- na.omit(paramSampleVec)
  meanParam = mean( paramSampleVec )
  medianParam = median( paramSampleVec )
  dres = density( paramSampleVec )
  modeParam = dres$x[which.max(dres$y)]
  mcmcEffSz = tryCatch({
    es <- round(effectiveSize(paramSampleVec), 1)
    unname(es)
  }, error = function(e) NA)

  MCSE = if (!is.na(mcmcEffSz)) sd(paramSampleVec)/sqrt(mcmcEffSz) else NA

  hdiLim = tryCatch({
    HDIofMCMC(paramSampleVec, credMass = credMass)
  }, error = function(e) c(NA, NA))
  return( c( Mean=meanParam , Median=medianParam , Mode=modeParam ,
             ESS=mcmcEffSz , MCSE = MCSE,
             HDImass=credMass , HDIlow=hdiLim[1] , HDIhigh=hdiLim[2]) )
}

# simplified version of a similar function in Kruschke 2015
#' @rdname HDIofMCMC
#' @export
#' @noRd
HDIofMCMC = function(
    sampleVec,
    credMass=0.95
    ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}


# max. diagnostics values
bay.ta.diagnostics <- function(x) {
  result <- data.frame(PSRF_max = max(x$`PSRF Point est.`[which(x$MCSE > 0)]),
                       PSRF_upper_max = max(x$`PSRF Upper C.I.`[which(x$MCSE > 0)]),
                       ESS_min = min(x$ESS[which(x$MCSE > 0)]
                       ))
  return(result)
}



#' @title gomp.a0
#'
#' @description
#' Internal function for generating starting values for the Gompertz model if the starting age is not
#' 15 years. Not run if the minimum age is actually 15. The original forumula
#' derives from ##.
#'
#' @param sampling integer. Number of sampling steps. Default: 100000.
#' @param b_min numeric. Minimum of Gompertz beta parameter. Default: 0.02.
#' @param b_max numeric. Maximum of Gompertz beta parameter. Default: 0.1.
#' @param minimum_age numeric. Minimum age in years. Default: 15.
#'
#' @return vector with coefficients for generating alpha and beta parameters
#' for Gompertz function.
#'
#'
#' @export
#' @noRd
#' @examples
#' NULL
#'
#' gomp.a0(minimum_age = 12)
#'
gomp.a0 <- function(
    sampling = 100000,
    b_min = 0.02,
    b_max = 0.1,
    minimum_age = 15) {

  # we do not want too much overhead so no computation if the default age of 15 is true
  if (minimum_age == 15) {
    fit_coeff <- c(-66.77, -2.324914, 0.0823)
  } else {
    null_age <- minimum_age - 15

    ind_df <- data.frame(b = runif(n = sampling, min = b_min, max = b_max)) |>
      dplyr::mutate(a = exp(rnorm(dplyr::n(), (-66.77 * (b - 0.0718) - 7.119), sqrt(0.0823) ))) |>
      dplyr::mutate(a0 = a * exp(b * null_age))

    fit <- lm(log(a0) ~ b, data = ind_df)
    rse <- sum(fit$residuals**2)/fit$df.residual # without squaring
    fit_coeff <- c(fit$coefficients[2], fit$coefficients[1], rse )
    fit_coeff <- unname(fit_coeff)
  }
  return(fit_coeff)
}


#' @title Extract thresholds
#'
#'@description
#' Extract thresholds
#'
#' @param x output from function `diagnostic.summary`
#'
#' @param mean_value a character string of either "Mode", "Median" or "Mode".
#' Default: "Mode".
#'
#' @return a matrix with threshold values of traits
#' @export
#'
#' @examples
#'NULL
#'
threshold.matrix <- function(x,
                             mean_value = "Mode")
  {
  rn <- rownames(x)
  sel <- grepl("^thresh_age\\[", rn)

  vals <- x[sel, mean_value]

  idx <- do.call(rbind, regmatches(rn[sel], gregexpr("\\d+", rn[sel])))
  idx <- matrix(as.integer(idx), ncol = 2)

  nrow <- max(idx[, 1], na.rm = TRUE)
  ncol <- max(idx[, 2], na.rm = TRUE)

  mat <- matrix(NA, nrow = nrow, ncol = ncol)
  for (k in seq_along(vals)) {
    mat[idx[k, 1], idx[k, 2]] <- vals[k]
  }
  return(mat)
}


#' @title Extract correlation matrix from Cholesky factor
#'
#'@description
#' As the LKJ prior uses the Cholesky decomposition of the correlation matrix,
#' getting the correlation indices from the coda chains is less straightforward
#' than it seems. It involves taking the cross product from the resulting coda
#' estimates.
#'
#' @param x matrix. Output from coda chains
#'
#' @return a matrix with correlations between traits
#' @export
#'
#' @examples
#'NULL
#'
extract.corr <- function(x) {
  samples_Ustar <- x[,grep("^Ustar\\[", colnames(x))]

  # Extract numbers inside the brackets
  numbers <- gsub("Ustar\\[|\\]", "", colnames(samples_Ustar))        # Remove "Ustar[" and "]"
  index_matrix <- do.call(rbind, strsplit(numbers, ", ")) # Split by comma and convert to matrix

  # Convert to numeric and find the maximum index
  index_matrix_numeric <- apply(index_matrix, 2, as.numeric)
  n_traits <- max(index_matrix_numeric)

  # Get number of posterior draws
  n_samples <- dim(samples_Ustar)[1]

  # Create array to hold correlation matrices
  corr_array <- array(NA, dim = c(n_traits, n_traits, n_samples))

  for (s in 1:n_samples) {
    # Extract Cholesky sample
    L <- matrix(samples_Ustar[s, ], nrow = n_traits, byrow = FALSE)
    # Reconstruct correlation matrix
    corr_array[,,s] <- crossprod(L)
  }

  # Compute posterior mean correlation matrix
  corr_mean <- apply(corr_array, c(1,2), mean)

  return(corr_mean)
}
