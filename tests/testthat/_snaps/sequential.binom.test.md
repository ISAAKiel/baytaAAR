# sequential.binom.test() produces correct output with default HDImass value

    Code
      sequential.binom.test(bay.ta_compare, known_age = spitalfields$Age)
    Output
        coverage n_in      perc    CI_low     CI_up   p_value
      1     0.95  164 0.9111111 0.8596638 0.9483365 0.0245341

# sequential.binom.test() produces correct output with several HDImass values

    Code
      sequential.binom.test(bay.ta_compare, known_age = spitalfields$Age, HDImass = c(
        0.75, 0.95))
    Output
        coverage n_in      perc    CI_low     CI_up   p_value
      1     0.75  128 0.7111111 0.6390067 0.7761161 0.2287218
      2     0.95  164 0.9111111 0.8596638 0.9483365 0.0245341

