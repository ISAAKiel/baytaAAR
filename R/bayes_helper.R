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
#' @param minimum_age numeric. Minimum age as specified in the
#' model. Default = NULL.
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
    minimum_age = NULL, # for M
    HDImass = 0.95,
    gelman_diag = TRUE,
    gelman_diag_multivariate = TRUE
    ) {
  mcmcMat = as.matrix(codaMCMClist,chains=TRUE)

  # add M from post-processing
  if(all(c("a","b") %in% colnames(mcmcMat)) & length(minimum_age) > 0) {
    M_draws <- 1/mcmcMat[,"b"] * log(mcmcMat[,"b"]/mcmcMat[,"a"]) +
      minimum_age
    mcmcMat <- cbind(mcmcMat, M = M_draws)
  }

  # add thresh_log and thresh_age from post-processing
  thresh_cols <- grep("^thresh\\[", colnames(mcmcMat), value=TRUE)
  thresh_len <- length(thresh_cols)
  beta0_cols <- grep("^beta0", colnames(mcmcMat), value=TRUE)

  for(col in thresh_cols) {
    # extract indices m,k
    idx <- gsub("thresh\\[|\\]", "", col)
    idx <- strsplit(idx, ",")[[1]]
    m <- idx[1]
    k <- idx[2]
    if (length(beta_cols) > 1) {
      beta_col  <- paste0("beta[", m, "]")
      beta0_col <- paste0("beta0[", m, "]")
    } else {
      beta0_col <- "beta0"
      beta_col <- "beta"
    }

    if(all(c(beta_col, beta0_col) %in% colnames(mcmcMat))) {
      thresh_log_draws <- (unlist(mcmcMat[,col]) - unlist(mcmcMat[,beta0_col])) /
        unlist(mcmcMat[,beta_col])
      new_name_log  <- paste0("thresh_log[", m, ",", k, "]")
      new_name_age  <- paste0("thresh_age[", m, ",", k, "]")
      mcmcMat <- cbind( mcmcMat, thresh_log_draws, exp(thresh_log_draws) )
      colnames(mcmcMat)[(ncol(mcmcMat)-1):ncol(mcmcMat)] <-
        c(new_name_log, new_name_age)
    }
  }
  parameterNames <- colnames(mcmcMat)

  summaryInfo = NULL
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName], credMass = HDImass ) )
    thisRowName = parName
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  summaryInfo_df <- as.data.frame(summaryInfo)
  if(gelman_diag == TRUE) {
    psrf_df <- as.data.frame((gelman.diag(chelsea_ps_list, multivariate = F))$psrf)
    thresh_NA <- data.frame(matrix(NA, 1 + 2 * thresh_len, 2))
    colnames(thresh_NA) <- colnames(psrf_df)
    rownames(thresh_NA) <- colnames(mcmcMat)[-c(1:(ncol(mcmcMat)- (1 + 2 * thresh_len)))]
    psrf_df <- rbind(psrf_df,thresh_NA)
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
  dres =  tryCatch({
    density( paramSampleVec )
    }, error = function(e) NA)
  modeParam = tryCatch({
    dres$x[which.max(dres$y)]
    }, error = function(e) NA)
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


#' @title Maximum and minimum diagnostic values
#'
#' @description
#' max. and min diagnostics values
#'
#' @rdname diagnostics.max.min
#' @export
#'
#' @examples
#' NULL
#'
diagnostics.max.min <- function(x) {
  result <- data.frame(PSRF_max = max(x$`PSRF Point est.`[which(x$MCSE > 0)]),
                       PSRF_upper_max = max(x$`PSRF Upper C.I.`[which(x$MCSE > 0)]),
                       ESS_min = min(x$ESS[which(x$MCSE > 0)]
                       ))
  return(result)
}


#' @title HDI median range of age estimates
#'
#'@description
#' Comparison with known age-at-death
#'
#' @param x a mcmc list.
#'
#' @param age_identifier a character string of either "age.s" or "age.s_c" to
#' select the uncalibrated or calibrated age estimates. Default: "age.s".
#'
#' @return a vector with mean, median and mode of the HDI range
#'
#' @export
#'
#' @examples
#'NULL
#'
hdi.agerange <- function(x, age_identifier = "age.s")
  {
  age_identifier_grep <- ifelse(age_identifier == "age.s",
                                "^age.s\\[", "^age.s_c")

  x_diag_red <- x[grep(age_identifier_grep,rownames(x)),]
  hdi_diff <- x_diag_red$HDIhigh - x_diag_red$HDIlow
  hdi_mean <- mean(hdi_diff)
  hdi_median <- median(hdi_diff)
  hdi_dens <- density( hdi_diff )
  hid_mode <- hdi_dens$x[which.max(hdi_dens$y)]
  return(c(hdi_mean, hdi_median, hid_mode))
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
    mat[mat < 0] <- NA
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
corr.mat.mean <- function(x) {
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


#' @title Summed or mean probability densities per category
#'
#'@description
#' Summing or averaging probability densities per category.
#'
#' @param mcmc_list a coda mcmc.list object.
#'
#' @param df_orig a data.frame of the original raw data
#'
#' @param group_col a string specifying the grouping category in `df_orig`.
#'
#' @param mode a string specifying the resulting data.frame of summed
#' probabilities or mean probabilities per category. Either `mean` or `summed`.
#'
#' @param age_identifier a character string of either "age.s" or "age.s_c" to
#' select the uncalibrated or calibrated age estimates. Default: "age.s".
#'
#' @return a data.frame with either probability summed by category or mean
#' per category.
#'
#' @export
#'
#' @examples
#'NULL
#'
prob.cat <- function(
    mcmc_list,
    age_identifier = "age.s",
    df_orig,
    group_col,
    mode = c("mean", "summed")
) {

  mode <- match.arg(mode)
  nChain <- length(mcmc_list)
  t_length <- nrow(df_orig)
  age_identifier_grep <- ifelse(age_identifier == "age.s",
                                "age.s[", "age.s_c[")

  df_orig$group_factor <- droplevels(as.factor(df_orig[, group_col]))
  df_orig$factor_cat <- as.numeric(df_orig$group_factor)
  levels_cat <- sort(unique(df_orig$factor_cat))

  if (mode == "mean") {

    # ---- Category-level posterior samples ----
    cat_post <- list()

    for (cLevels in levels_cat) {
      cat_samples <- NULL

      for (cIdx in 1:nChain) {
        chain_samples <- NULL

        for (t in 1:t_length) {
          if (cLevels == df_orig$factor_cat[t]) {
            param_name <- paste0(age_identifier_grep, t, "]")
            chain_samples <- cbind(chain_samples,
                                   mcmc_list[, param_name][[cIdx]])
          }
        }

        cat_samples <- c(cat_samples, rowMeans(chain_samples))
      }

      cat_post[[as.character(cLevels)]] <- cat_samples
    }

    dense_xy <- data.frame(
      category = factor(rep(levels(df_orig$group_factor),
                            each = length(cat_post[[1]])),
                        levels = levels(df_orig$group_factor)),
      value = unlist(cat_post)
    )

    return(dense_xy)

  } else {

    # ---- Mixture density per category ----
    xMat <- NULL
    yMat <- NULL

    for (cLevels in levels_cat) {
      coda_object_simplified <- NULL

      for (t in 1:t_length) {
        if (cLevels == df_orig$factor_cat[t]) {
          for (cIdx in 1:nChain) {
            param_name <- paste0(age_identifier_grep, t, "]")
            onecolumn <- mcmc_list[, param_name][[cIdx]]
            coda_object_simplified <- rbind(coda_object_simplified, onecolumn)
          }
        }
      }

      densInfo <- density(coda_object_simplified)
      xMat <- cbind(xMat, densInfo$x)
      yMat <- cbind(yMat, densInfo$y)
    }

    xMat_melt <- reshape::melt(as.data.frame(xMat), id.vars = NULL)
    yMat_melt <- reshape::melt(as.data.frame(yMat), id.vars = NULL)
    dense_xy <- cbind(xMat_melt, yMat_melt[, 2])
    colnames(dense_xy) <- c("category", "x", "y")

    levels(dense_xy$category) <- levels(df_orig$group_factor)

    df_orig_cat_n <- df_orig %>%
      dplyr::group_by(group_factor) %>%
      dplyr::summarize(n = n(), .groups = "drop")

    dense_xy$y_prop <- NULL
    for (i in 1:length(dense_xy$y)) {
      dense_xy$y_prop[i] <-
        as.numeric(df_orig_cat_n[as.numeric(dense_xy[i, 1]), 2] *
                     dense_xy$y[i]) / nrow(df_orig_cat_n)
    }

    return(dense_xy)
  }
}

