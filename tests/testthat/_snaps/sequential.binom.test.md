# sequential.binom.test() produces correct output with default HDImass value

    Code
      sequential.binom.test(bay.ta_compare, known_age = spitalfields$Age)
    Output
        coverage n_in      perc    CI_low     CI_up      p_value
      1     0.95  152 0.8444444 0.7830858 0.8940772 1.057047e-07

# sequential.binom.test() produces correct output with several HDImass values

    Code
      sequential.binom.test(bay.ta_compare, known_age = spitalfields$Age, HDImass = c(
        0.75, 0.95))
    Output
        coverage n_in      perc    CI_low     CI_up      p_value
      1     0.75  114 0.6333333 0.5584181 0.7037779 5.328385e-04
      2     0.95  152 0.8444444 0.7830858 0.8940772 1.057047e-07

