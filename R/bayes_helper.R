#' @title Diagnostic summary of MCMC samples
#'
#' @description
#' Summarising diagnostics from a \code{coda::mcmc.list}, partly derived from
#'  \href{https://jkkweb.sitehost.iu.edu/DoingBayesianDataAnalysis/}{Kruschke
#'  (2015)}.
#'
#' @details
#'  Because the first threshold is fixed, the Gelman-Rubin multivariate PSRF
#'  will always throw an error, so this is automatically set to \code{FALSE}. If
#'  the gelman diagnostics still produce an error, deactivate \code{gelman_diag}
#'  altogether by setting it to \code{FALSE}, too.
#'
#' @inheritParams age.comp.summary
#' @param HDImass numeric. Value within 0 and 1. Default = 0.95.
#' @param gelman_diag logical. If TRUE, the Gelman-Rubin diagnostics for computing
#' the PSRF is invoked. Default: TRUE.
#'
#' @return
#' A data.frame of class \code{diagnostic_summary} with the row names according
#' to the parameters to be monitored and the following numeric columns:
#' \itemize{
#' \item{ \code{PSRF Point est.} \emph{Potential scale reduction factor}
#' (= Gelman-Rubin statistic), a measure of the mixing of chains. }
#' \item{ \code{PSRF Upper C.I.} The upper limit of the 0.95-confidence interval
#' of the PSRF.}
#' \item{ \code{Mean} Arithmetic mean of the estimates. }
#' \item{ \code{Median} Median of the estimates. }
#' \item{ \code{Mode} Mode of the estimates. }
#' \item{ \code{ESS} \emph{Effective sample size}, a control of autocorrelation.}
#' \item{ \code{MCSE} \emph{Monte Carlo standard error}.}
#' \item{ \code{HDImass} Credibility level of the \emph{highest density
#' interval}.}
#' \item{ \code{HDIlow} sStart of the \emph{highest density interval}.}
#' \item{ \code{HDIhigh} End of the \emph{highest density interval}.}
#' }
#'
#' @rdname diagnostic.summary
#' @export
#'
#' @examplesIf interactive()
#' # select Sorsum data with auricular surface after Lovejoy et al. 1985 and
#' # convert to matrix
#' sorsum <- as.matrix(sorsum_as[,2])
#'
#' # example with default settings, please be patient
#' sorsum_res <- bay.ta(method = sorsum)
#'
#' # compute diagnostics of the MCMC samples
#' sorsum_diag <- diagnostic.summary(sorsum_res)
#'
#' # show first rows
#' head(sorsum_diag)
#'
diagnostic.summary <- function(mcmc_list, HDImass = 0.95, gelman_diag = TRUE)
  {
  checkmate::assertClass(mcmc_list, "mcmc.list")
  checkmate::assertNumeric(HDImass, lower = 0, upper = 1, len = 1)
  checkmate::assertLogical(gelman_diag)

  parameterNames = coda::varnames(mcmc_list)
  mcmcMat = as.matrix(mcmc_list,chains=TRUE)
  summaryInfo = NULL
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName],
                                                      credMass = HDImass ) )
    thisRowName = parName
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  summaryInfo_df <- as.data.frame(summaryInfo)
  if(gelman_diag == TRUE) {
    psrf_df <- as.data.frame((coda::gelman.diag(mcmc_list, multivariate =
                                            FALSE))$psrf)
    colnames(psrf_df) <- c("PSRF Point est.", "PSRF Upper C.I.")
    diagnostic_summary <- cbind(psrf_df, summaryInfo_df)
  }  else {
    diagnostic_summary <- summaryInfo_df
  }
  class(diagnostic_summary) <- c("diagnostic_summary", class(diagnostic_summary))
  return(diagnostic_summary)
}

# simplified version of a similar function in Kruschke 2015
#' @rdname summarizePost
#' @noRd
summarizePost = function(
    paramSampleVec,
    credMass=0.95
    ) {
  paramSampleVec <- stats::na.omit(paramSampleVec)
  meanParam = mean( paramSampleVec )
  medianParam = stats::median( paramSampleVec )
  dres =  tryCatch({
    stats::density( paramSampleVec )
    }, error = function(e) NA)
  modeParam = tryCatch({
    dres$x[which.max(dres$y)]
    }, error = function(e) NA)
  mcmcEffSz = tryCatch({
    es <- round(coda::effectiveSize(paramSampleVec), 1)
    unname(es)
    }, error = function(e) NA)

  MCSE = if (!is.na(mcmcEffSz)) stats::sd(paramSampleVec)/sqrt(mcmcEffSz) else NA

  hdiLim = tryCatch({
    HDIofMCMC(paramSampleVec, credMass = credMass)
  }, error = function(e) c(NA, NA))
  return( c( Mean=meanParam , Median=medianParam , Mode=modeParam ,
             ESS=mcmcEffSz , MCSE = MCSE,
             HDImass=credMass , HDIlow=hdiLim[1] , HDIhigh=hdiLim[2]) )
}

# simplified version of a similar function in Kruschke 2015
#' @rdname HDIofMCMC
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
#' Convenience function to quickly extract maximum and mininum diagnostics
#' values of the function \code{diagnostics.summary()} over all parameters. The
#' maximum values of the PSRF should be below 1.1 while the minumum ESS should
#' be above 10,000. If either of this is not the case, consider to increase the
#' length of the chains, i. e. the number of iterations.
#'
#' @param x output from function \code{diagnostics.summary()}
#'
#' @return
#' A data.frame with one row and the following numeric columns:
#' \itemize{
#' \item{ \code{PSRF_max} Maximum value of the \emph{potential scale reduction
#' factor}. }
#' \item{ \code{PSRF_upper_max} Mximum value of the upper limit of the
#' 0.95-confidence interval of the PSRF.}
#' \item{ \code{ESS_min} Minimum of the \emph{effective sample size}. }
#' }
#'
#' @rdname diagnostics.max.min
#' @export
#'
#' @examplesIf interactive()
#' # select Sorsum data with auricular surface after Lovejoy et al. 1985 and
#' # convert to matrix
#' sorsum <- as.matrix(sorsum_as[,2])
#'
#' # example with default settings, please be atient
#' sorsum_res <- bay.ta(method = sorsum)
#'
#' # compute diagnostics of the MCMC samples
#' sorsum_diag <- diagnostic.summary(sorsum_res)
#'
#' # show maximum and minimum values
#' diagnostics.max.min(sorsum_diag)
#'
diagnostics.max.min <- function(x) {
  checkmate::assertClass(x, "diagnostic_summary")

  result <-
    data.frame(PSRF_max = max(x$`PSRF Point est.`[which(x$MCSE > 0)]),
               PSRF_upper_max = max(x$`PSRF Upper C.I.`[which(x$MCSE > 0)]),
                       ESS_min = min(x$ESS[which(x$MCSE > 0)]
                       ))
  return(result)
}


#' @title Summary of age estimates
#'
#'@description
#' Convenience function to quickly extract the age-related estimates from the
#' result of the function \code{diagnostics.summary()}.
#'
#' @param x output from function \code{diagnostics.summary()}
#'
#' @inheritParams age.comp.summary
#'
#' @return A data.frame with the chosen mean measure and the HDI ranges as
#' specified in the output of \code{diagnostics.summary()} as columns and the
#' following rows:
#' \itemize{
#' \item{ \code{age_mean} Mean of the mean ages.}
#' \item{ \code{hdi_diff} Mean of the \emph{highest density intervals}.}
#' \item{ \code{b} Mean of the Gompertz parameter \eqn{\beta}. }
#' \item{ \code{a} Mean of the Gompertz parameter \eqn{\alpha}. }
#' \item{ \code{M} Modal age, derived from the Gompertz parameters
#' \eqn{\alpha} and \eqn{\beta} according to the equation
#' (1 / \eqn{\beta}) * log(\eqn{\beta} / \eqn{\alpha}) + minimum_age. }
#' }
#'
#' @export
#'
#' @examplesIf interactive()
#' # select Sorsum data with auricular surface after Lovejoy et al. 1985 and
#' # convert to matrix
#' sorsum <- as.matrix(sorsum_as[,2])
#'
#' # example with default settings, please be a little bit patient
#' sorsum_res <- bay.ta(method = sorsum)
#'
#' sorsum_diag <- diagnostic.summary(sorsum_res)
#'
#' # show summary of age-related estimates
#' age.estim.summary(sorsum_diag)
#'
age.estim.summary <- function(x,
                              mean_choice = "Mode",
                              age_identifier = "age.s")
  {
  checkmate::assertClass(x, "diagnostic_summary")
  checkmate::assertChoice(age_identifier, c("age.s", "age.s_c"))
  checkmate::assertChoice(mean_choice, c("Mode", "Median", "Mean"))

  age_identifier_grep <- ifelse(age_identifier == "age.s",
                                "^age.s\\[", "^age.s_c")
  hdi_mass <- (1 - x["a", "HDImass"]) / 2
  a <- x["a", mean_choice]
  a_low <- x["a", "HDIlow"]
  a_high <- x["a", "HDIhigh"]
  b <- x["b", mean_choice]
  b_low <- x["b", "HDIlow"]
  b_high <- x["b", "HDIhigh"]
  M <- x["M", mean_choice]
  M_low <- x["M", "HDIlow"]
  M_high <- x["M", "HDIhigh"]
  x_diag_red <- x[grep(age_identifier_grep, rownames(x)),]

  hdi_diff <- x_diag_red$HDIhigh - x_diag_red$HDIlow
  hdi_dens <- stats::density( hdi_diff )
  age_dens <- stats::density( x_diag_red[,mean_choice] )

  if (mean_choice == "Mode") {
    hdi <- hdi_dens$x[which.max(hdi_dens$y)]
    age_mean <- age_dens$x[which.max(age_dens$y)]
  } else if (mean_choice == "Median") {
    hdi <- stats::median(hdi_diff)
    age_mean <- stats::median(x_diag_red[,mean_choice])
  } else {
    hdi <- mean(hdi_diff)
    age_mean <- mean(x_diag_red[,mean_choice])
  }
  hdi_low <- stats::quantile(hdi_diff, probs = c( hdi_mass))
  hdi_high <- stats::quantile(hdi_diff, probs = c( 1 - hdi_mass))
  age_low <- stats::quantile(x_diag_red[,mean_choice], probs = c( hdi_mass))
  age_high <- stats::quantile(x_diag_red[,mean_choice], probs = c( 1 -hdi_mass))

  age_result <- data.frame(rbind(M = cbind(M, M_low, M_high),
                                 age_mean = cbind(age_mean, age_low, age_high),
                                 b = cbind(b, b_low, b_high),
                                 a = cbind(a, a_low, a_high),
                                 hdi = cbind(hdi, hdi_low, hdi_high)))
  colnames(age_result) <- c(mean_choice, hdi_mass, 1 - hdi_mass)
  rownames(age_result) <- c("age_mean", "hdi_diff", "b", "a", "M")
  return(age_result)
  }


#' @title gomp.a0
#'
#' @description
#' Internal function for generating starting values for the Gompertz model if
#' the starting age is not 15 years. Not run if the minimum age is actually 15.
#' The original formula derives from ##.
#'
#' @param sampling integer. Number of sampling steps. Default: 100000.
#' @param b_min numeric. Minimum of Gompertz \eqn{\beta} parameter. Default: 0.02.
#' @param b_max numeric. Maximum of Gompertz \eqn{\beta} parameter. Default: 0.1.
#' @param minimum_age numeric. Minimum age in years. Default: 15.
#'
#' @return vector with coefficients for generating \eqn{\alpha} and \eqn{\beta}
#' parameters for Gompertz function.
#'
#' @noRd
#' @examples
#'
#' gomp.a0()
#'
gomp.a0 <- function(
    sampling = 100000,
    b_min = 0.02,
    b_max = 0.1,
    minimum_age = 15) {

  # we do not want too much overhead so no computation if the default age == 15
  if (minimum_age == 15) {
    fit_coeff <- c(-66.76844784, -2.32502545, 0.0823)
  } else {
    null_age <- minimum_age - 15

    ind_df <- data.frame(b = stats::runif(n = sampling, min = b_min, max = b_max)) |>
      dplyr::mutate(a = exp(stats::rnorm(dplyr::n(),
                                  (-66.76844784 * (b - 0.0718) - 7.119),
                                  sqrt(0.0823) ))) |>
      dplyr::mutate(a0 = a * exp(b * null_age))

    fit <- stats::lm(log(a0) ~ b, data = ind_df)
    rse <- sum(fit$residuals**2)/fit$df.residual # without squaring
    fit_coeff <- c(fit$coefficients[2], fit$coefficients[1], rse )
    fit_coeff <- unname(fit_coeff)
  }
  return(fit_coeff)
}


#' @title Compute thresholds for chains
#'
#'@description
#' The computation of the thresholds on the log- and the age-scale is done
#' outside of the MCMC simulation to spare the computation cost and the memory.
#' The function returns a \code{coda::mcmc.list} which can be further processed.
#'
#' @inheritParams age.comp.summary
#'
#' @return A \code{coda::mcmc.list} for threshold values of traits on the
#' age scale.
#' @export
#'
#' @examplesIf interactive()
#' # select Sorsum data with auricular surface after Lovejoy et al. 1985 and
#' # convert to matrix
#' sorsum <- as.matrix(sorsum_as[,2])
#'
#' # example with default settings, please be patient
#' sorsum_res <- bay.ta(method = sorsum)
#'
#' # compute threshold chains
#' threshold_chains <- threshold.chains(sorsum_res)
#'
threshold.chains <- function(mcmc_list) {
  checkmate::assertClass(mcmc_list, "mcmc.list")

  out <- lapply(mcmc_list, function(chain) {
    s <- as.matrix(chain)

    # identify columns
    thresh_cols <- grep("^thresh\\[", colnames(s), value = TRUE)
    beta0_cols  <- grep("^beta0\\[", colnames(s), value = TRUE)
    beta_cols   <- grep("^beta\\[", colnames(s), value = TRUE)

    if (length(beta_cols) == 0) {
      beta_cols  <- "beta"
      beta0_cols <- "beta0"
    }

    n_methods <- length(beta_cols)

    # extract indices from thresh[i,j]
    thresh_idx <- do.call(rbind, regmatches(thresh_cols,
                                            gregexpr("\\d+", thresh_cols)))
    thresh_idx <- matrix(as.integer(thresh_idx), ncol = 2)

    thresholds_list <- vector("list", n_methods)

    for (m in seq_len(n_methods)) {

      beta0 <- s[, beta0_cols[m]]
      beta  <- s[, beta_cols[m]]

      # select ALL thresholds for method m
      cols_m <- thresh_cols[thresh_idx[,1] == m]

      if (length(cols_m) == 0) next

      thresh_m <- s[, cols_m, drop = FALSE]

      # drop NA-only thresholds (ragged structure cleanup)
      keep <- !is.na(thresh_m[1, ])
      thresh_m <- thresh_m[, keep, drop = FALSE]

      if (ncol(thresh_m) == 0) next

      # vectorized transform
      thresholds_m <- exp((thresh_m - beta0) / beta)

      # naming
      n_thresh <- ncol(thresholds_m)
      colnames(thresholds_m) <-
        paste0("thresh_age[", m, ",", seq_len(n_thresh), "]")

      thresholds_list[[m]] <- thresholds_m
    }

    thresholds_mat <- do.call(cbind, thresholds_list)

    coda::as.mcmc(
      thresholds_mat,
      start = stats::start(chain),
      end   = stats::end(chain),
      thin  = coda::thin(chain)
    )
  })

  return(coda::as.mcmc.list(out))
}


#' @title Extract thresholds
#'
#'@description
#' A convenience function to extract mean thresholds values from the output of
#' \code{diagnostic.summary()} which in turn was derived from a
#' \code{coda::mcmc.list} computed with \code{threshold.chains()}
#'
#' @param x output from function \code{diagnostic.summary()}
#'
#' @inheritParams age.comp.summary
#'
#' @return A matrix with threshold values of traits. The number of rows
#' corresponds to the number of traits, and the number of columns to the
#' maximum number of levels of one of the traits.
#'
#' @export
#'
#' @examplesIf interactive()
#' # select Sorsum data with auricular surface after Lovejoy et al. 1985 and
#' # convert to matrix
#' sorsum <- as.matrix(sorsum_as[,2])
#'
#' # example with default settings, please be patient
#' sorsum_res <- bay.ta(method = sorsum)
#'
#' # compute threshold chains
#' threshold_chains <- threshold.chains(sorsum_res)
#'
#' # compute summary diagnostics
#' threshold_diag <- diagnostic.summary(threshold_chains)
#'
#' # extract threshold matrix (for sorsum only 1 row)
#' threshold.matrix(threshold_diag)
#'
threshold.matrix <- function(
    x,
    mean_choice = "Mode")
  {
  checkmate::assertClass(x, "diagnostic_summary")
  checkmate::assertChoice(mean_choice, c("Mode", "Median", "Mean"))
  rownames_x <- substr(rownames(x), 1, 10)
  checkmate::assert_names(rownames_x, must.include = "thresh_age")

  rn <- rownames(x)
  sel <- grepl("^thresh_age\\[", rn)

  vals <- x[sel, mean_choice]

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
#' As the LKJ prior for the correlation matrix uses the Cholesky decomposition
#' of the correlation matrix, getting the correlation indices from the coda
#' chains is less straightforward than it seems. It involves taking the cross
#' product from the resulting coda estimates.
#'
#' @inheritParams age.comp.summary
#'
#' @return A symmetric matrix with correlations between traits. The number of
#' rows and columns corresponds to the number of traits.
#'
#' @export
#'
#' @examplesIf interactive()
#'
#'   # select Spitalfields data with multiple traits and convert to matrix
#'   spitalfields_traits <- as.matrix(spitalfields[,c(2:6)])
#'
#'   # example with multinormal likelihood, please be patient
#'   spitalfields_res <- bay.ta(algorithm = "mnorm",
#'   method = spitalfields_traits)
#'
#'   # compute correlation matrix
#'   corr.mat.mean(spitalfields_res)
#'
corr.mat.mean <- function(mcmc_list) {
  checkmate::assertClass(mcmc_list, "mcmc.list")

  x_matrix <- as.matrix(mcmc_list, chains=TRUE)

  samples_Ustar <- x_matrix[,grep("^Ustar\\[", colnames(x_matrix))]

  # Extract numbers inside the brackets
  numbers <- gsub("Ustar\\[|\\]", "", colnames(samples_Ustar))
  index_matrix <- do.call(rbind, strsplit(numbers, ", "))

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
#' Summing or averaging probability densities per category. The resulting
#' data.frames can be used, for example, to produce illustrative diagrams. See
#' the vignettes for some examples.
#'
#' @param group_vec a vector specifying the grouping category.
#'
#' @param mode a string specifying the resulting data.frame of summed
#' probabilities or mean probabilities per category. Either \code{mean} or
#' \code{summed}.
#'
#' @inheritParams age.comp.summary
#'
#' @return A data.frame with either probability summed by category or mean
#' per category.
#'
#' @export
#'
#' @examplesIf interactive()
#'
#'   # select Spitalfields data with multiple traits and convert to matrix
#'   spitalfields_traits <- as.matrix(spitalfields[,c(2:6)])
#'
#'   # example with multinormal likelihood, please be patient
#'   spitalfields_res <- bay.ta(framework = "NIMBLE", algorithm = "mnorm",
#'   method = spitalfields_traits)
#'
#'   # compute averaging probabilities per category Sex
#'   prob_cat_mean <- prob.cat(spitalfields_res, group_vec = spitalfields$Sex,
#'   mode = "mean")
#'
#'   # compute summed probabilities per category Sex
#'   prob_cat_summed <- prob.cat(spitalfields_res, group_vec = spitalfields$Sex,
#'   mode = "summed")
#'
prob.cat <- function(
    mcmc_list,
    age_identifier = "age.s",
    group_vec,
    mode = c("mean", "summed")
) {
  checkmate::assertClass(mcmc_list, "mcmc.list")
  checkmate::assertChoice(age_identifier, c("age.s", "age.s_c"))
  checkmate::assertVector(group_vec)
  checkmate::assertChoice(mode, c("mean", "summed"))

  mode <- match.arg(mode)
  nChain <- length(mcmc_list)
  t_length <- length(group_vec)
  age_identifier_grep <- ifelse(age_identifier == "age.s",
                                "age.s[", "age.s_c[")

  group_factor <- droplevels(factor(group_vec))
  factor_cat <- as.numeric(group_factor)
  levels_cat <- sort(unique(factor_cat))

  if (mode == "mean") {

    # ---- Category-level posterior samples ----
    cat_post <- list()

    for (cLevels in levels_cat) {
      cat_samples <- NULL

      for (cIdx in 1:nChain) {
        chain_samples <- NULL

        for (t in 1:t_length) {
          if (cLevels == factor_cat[t]) {
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
      category = factor(rep(levels(group_factor),
                            each = length(cat_post[[1]])),
                        levels = levels(group_factor)),
      value = unlist(cat_post)
    )
  } else {

    # ---- Mixture density per category ----
    xMat <- NULL
    yMat <- NULL

    for (cLevels in levels_cat) {
      coda_object_simplified <- NULL

      for (t in 1:t_length) {
        if (cLevels == factor_cat[t]) {
          for (cIdx in 1:nChain) {
            param_name <- paste0(age_identifier_grep, t, "]")
            onecolumn <- mcmc_list[, param_name][[cIdx]]
            coda_object_simplified <- rbind(coda_object_simplified, onecolumn)
          }
        }
      }

      densInfo <- stats::density(coda_object_simplified)
      xMat <- cbind(xMat, densInfo$x)
      yMat <- cbind(yMat, densInfo$y)
    }

    xMat_melt <- tidyr::pivot_longer(as.data.frame(xMat),
                                     cols = tidyr::everything(),
                                     names_to = "variable",
                                     values_to = "value")
    yMat_melt <- tidyr::pivot_longer(as.data.frame(yMat),
                                     cols = tidyr::everything(),
                                     names_to = "variable",
                                     values_to = "value")
    dense_xy <- cbind(xMat_melt, yMat_melt[, 2])
    colnames(dense_xy) <- c("category", "x", "y")

    dense_xy$category <- as.factor(dense_xy$category)
    levels(dense_xy$category) <- levels(group_factor)
    cat_counts <- as.matrix(table(group_factor))

    df_orig_cat_n <- as.data.frame(cat_counts)
    df_orig_cat_n$group_factor <- rownames(df_orig_cat_n)
    df_orig_cat_n <- df_orig_cat_n[,c(2,1)]
    colnames(df_orig_cat_n) <- c("group_factor", "n")

    dense_xy$y_prop <- NULL
    for (i in 1:length(dense_xy$y)) {
      dense_xy$y_prop[i] <-
        as.numeric(df_orig_cat_n[as.numeric(dense_xy[i, 1]), 2] *
                     dense_xy$y[i]) / nrow(df_orig_cat_n)
    }
  }
  return(dense_xy)
}

