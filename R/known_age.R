#' @title Quality measures of age estimation
#'
#'@description
#' Comparison of estimated age with known age-at-death.
#'
#' @param mcmc_list MCMC output from coda chains.
#'
#' @param known_age a vector of known age-at-death. NAs are allowed and those
#' individuals will subsequently be ignored.
#'
#' @param mean_choice a character string of either "Mode", "Median" or "Mode".
#' Default: "Mode".
#'
#' @param age_identifier a character string of either "age.s" or "age.s_c" to
#' select the uncalibrated or calibrated age estimates. Default: "age.s".
#'
#' @return a dataframe with age estimation quality parameters
#'
#' @export
#'
#' @examples
#'NULL
#'
age.comp.summarize <- function(mcmc_list,
                               known_age,
                               mean_choice = "Mode",
                               age_identifier = "age.s",
                              ...) {
  checkmate::assertClass(mcmc_list, "mcmc.list")
  checkmate::assertChoice(age_identifier, c("age.s", "age.s_c"))
  checkmate::assertAtomic(known_age, all.missing = FALSE)
  checkmate::assertChoice(mean_choice, c("Mean", "Median", "Mode"))

  age_identifier_grep <- ifelse(age_identifier == "age.s",
                                "^age.s\\[", "^age.s_c")
  idx <- which(!is.na(known_age))
  x_mcmcMat = as.matrix(mcmc_list, chains=TRUE)
  x_diag <- diagnostic.summary(mcmc_list, ...)
  x_diag_red <- x_diag[grep(age_identifier_grep,rownames(x_diag)),][idx,]
  ages <- x_mcmcMat[,grep(age_identifier_grep,colnames(x_mcmcMat))][,idx]
  known_age <- known_age[idx]

  estimated_age <- x_diag_red[,mean_choice]
  Residual_model  <-  lm((known_age - estimated_age) ~ known_age)

  corrPearson <- cor.test(estimated_age, known_age, method="pearson")

  tmnlp_res <- tmnlp(known_age, ages)
  crps_res <- mean(scoringRules::crps_sample(known_age, t(ages)) )

  age_estimation <- data.frame(Mean_estimated = mean(estimated_age),
                               Bias = mean(estimated_age - known_age),
                               corrPearson = corrPearson$estimate,
                               corr_p = corrPearson$p.value,
                               Residual_age_slope =
                                 Residual_model$coefficients[2],
                               Inaccuracy =
                                 mean(abs(estimated_age - known_age)),
                               RMSE = Metrics::rmse(known_age, estimated_age),
                               TMNLP = tmnlp_res,
                               CRPS = crps_res)
    HDIhigh <- x_diag_red[,"HDIhigh"]
    HDIlow <- x_diag_red[,"HDIlow"]
    age_estimation$Coverage <-
      sum(ifelse(known_age >= HDIlow - 0.05 & known_age <=
                   HDIhigh + 0.05, 1, 0)) / length(known_age) * 100
    age_estimation$HDI_Diff_median <- stats::quantile(HDIhigh - HDIlow,
                                                      probs = c(0.5))
    age_estimation$HDI_Diff_quant_025 <- stats::quantile(HDIhigh - HDIlow,
                                                         probs = c( 0.025))
    age_estimation$HDI_Diff_quant_975 <- stats::quantile(HDIhigh - HDIlow,
                                                         probs = c(0.975))
  return(age_estimation)
}


#' @title TMNLP
#'
#' @description
#' Internal function for computing the TMNLP.
#'
#' @inheritParams age.comp.summarize
#' @param mcmcMat MCMC matrix.
#'
#' @return numeric. Value of TMNLP, smaller is better.
#'
#'
#' @export
#' @noRd
#' @examples
#' NULL
#'
tmnlp <- function(known_age, mcmcMat) {
  x_length <- length(known_age)
  log_vals <- numeric(x_length)

  for (i in seq_len(x_length)) {
    age_dens <- density(mcmcMat[i, ], n = 512 * 8)
    vec_i <- approx(age_dens$x, age_dens$y, xout = known_age[i], rule = 2)$y
    vec_i <- pmax(vec_i, .Machine$double.eps)
    log_vals[i] <- log(vec_i)
  }
  tmnlp <- round(-mean(log_vals), 4)
  return(tmnlp)
}


#' @title Plots of quality measures of age estimation
#'
#'@description
#' Comparison of estimated age with known age-at-death.
#'
#' @param x output from the function `diagnostic.summary`
#'
#' @inheritParams age.comp.summarize
#'
#' @return a ggplot object with 2 x 2 single plots, showing
#'
#' @export
#'
#' @examples
#'NULL
#'
age.comp.plot <- function(x,
                        age_identifier = "age.s",
                        known_age,
                        mean_choice = "Mode" ) {
  checkmate::assertClass(x, "diagnostic_summary")
  checkmate::assertChoice(age_identifier, c("age.s", "age.s_c"))
  checkmate::assertAtomic(known_age, all.missing = FALSE)
  checkmate::assertChoice(mean_choice, c("Mean", "Median", "Mode"))

  x$chosen_mean <- x[,mean_choice]
  age_identifier_grep <- ifelse(age_identifier == "age.s",
                                "^age.s\\[", "^age.s_c")
  idx <- which(!is.na(known_age))
  x_red <- x[grep(age_identifier_grep,rownames(x)),][idx,]
  age_min <- round(min(x_red$HDIlow))

  x_red$known_age <- known_age[idx]
  known_age_density <- density(x_red$known_age, bw = 5)
  known_age_density_df <- data.frame(x = known_age_density$x,
                                     y = known_age_density$y)
  x_length <- length(x_red$known_age)
  x_ordered <- x_red[order(x_red$known_age),]
  x_ordered$id <- c(1:x_length)
  alpha_mean <- x['a',"Mean"]
  beta_mean <- x['b',"Mean"]

  library(ggplot2)

  plot1 <- ggplot(x_ordered, aes(x = id)) +
    geom_errorbar(aes(ymin=HDIlow, ymax=HDIhigh,
                      color= known_age > HDIlow - 0.05 &
                        known_age < HDIhigh + 0.05), lwd = 0.5 ,width = 0) +
    geom_point(aes(y=known_age), shape = 3, colour = "black") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_colour_manual(name = 'Age in range',
                        values = setNames(c('chartreuse4','coral2'),c(T, F))) +
    xlab("\nIndividuals ordered by known age-at-death") +
    ylab("HDIlow to HDIhigh, known age-at-death\n") +
    theme_light()

  plot2 <- ggplot() +
    geom_line(data = known_age_density_df,
              aes(x = x, y = y, col = "density of actual ages\n(bw = 5)\n")) +
    xlim(age_min, 100) +
    geom_function(fun = function(x) flexsurv::dgompertz(x - age_min,
                                                        beta_mean, alpha_mean),
                  aes(col = "Gompertz parameters\nfrom estimates")) +
    theme_light() +
    xlab("age-at-death") +
    ylab("Density\n") +
    scale_colour_manual(values = c("red","black")) +
    theme( legend.title = element_blank(),
           legend.spacing.y = unit(1.0, 'cm')) +
    guides(fill = guide_legend(byrow = F))

  plot3 <- ggplot (x_ordered, aes(x = known_age, y = chosen_mean)) +
    geom_point(shape = 21) + xlim(15,95) + ylim(15,95) +
    geom_smooth(method='lm', formula= y~x) + theme_light() +
    xlab("\nknown age-at-death") +
    ylab(paste0(mean_choice, " of estimated age-at-death\n")) +
    geom_abline(slope = 1, intercept = 0, linetype = 3)

  plot4 <- ggplot(x_ordered, aes(x=known_age, y = ( known_age - chosen_mean))) +
    geom_point(shape = 21) + geom_smooth(method='lm', formula= y~x) +
    geom_hline(yintercept = 0, linetype = 3) + theme_light() +
    xlab("\nknown age-at-death") + ylab("Residual\n")

  plot_result <- ggpubr::ggarrange( plot1, plot2, plot3, plot4, nrow = 2,
                                    ncol = 2)
   return(plot_result)
}

#' @title Sequential output of cumulative binomial test
#'
#'@description
#' Sequential output of cumulative binomial test
#'
#' @param HDImass a numeric or a vector with the probability range.
#'
#' @inheritParams age.comp.summarize
#'
#' @return a data.frame
#'
#' @export
#'
#' @examples
#'NULL
#'
sequential.binom.test <- function(mcmc_list,
                                  known_age,
                                  HDImass = 0.95,
                                  age_identifier = "age.s") {
  checkmate::assertClass(mcmc_list, "mcmc.list")
  checkmate::assertAtomic(known_age, all.missing = FALSE)
  checkmate::assertAtomic(HDImass, all.missing = FALSE, any.missing = FALSE,
                          unique = TRUE)
  checkmate::assertChoice(age_identifier, c("age.s", "age.s_c"))

  result_df <- data.frame()
  age_identifier_grep <- ifelse(age_identifier == "age.s",
                                "^age.s\\[", "^age.s_c")

  idx <- which(!is.na(known_age))
  known_age <- known_age[idx]

  for (i in HDImass) {
    MCMC_diag  <-  diagnostic.summary(mcmc_list, HDImass = i, gelman_diag = F)
    MCMC_diag_age <-
      MCMC_diag[grep(age_identifier_grep,rownames(MCMC_diag)),][idx,]
    MCMC_diag_age$known_age <- known_age

    MCMC_diag_age$in_HDI <-
      MCMC_diag_age$known_age >= MCMC_diag_age$HDIlow - 0.05 &
      MCMC_diag_age$known_age <= MCMC_diag_age$HDIhigh + 0.05
    n_total <- nrow(MCMC_diag_age)
    n_in <- sum(MCMC_diag_age$in_HDI)

    # Binomial test: is observed coverage significantly different from i?
    binom_coverage <- binom.test(n_in, n_total, p = i)
    binom_result <- data.frame(coverage = i, n_in = n_in,
                               perc = unname(binom_coverage$estimate),
                               CI_low = binom_coverage$conf.int[1],
                               CI_up = binom_coverage$conf.int[2],
                               p_value = binom_coverage$p.value)
    result_df <- rbind(result_df, binom_result)
  }
  return(result_df)
}
