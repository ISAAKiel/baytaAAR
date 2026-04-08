# sequential.binom.test() produces correct output with default HDImass value

    Code
      sequential.binom.test(bay.ta_compare, known_age = spitalfields$Age)
    Output
        coverage n_in      perc    CI_low     CI_up      p_value
      1     0.95  151 0.8388889 0.7768829 0.8893679 2.874477e-08

# sequential.binom.test() produces correct output with several HDImass values

    Code
      sequential.binom.test(bay.ta_compare, known_age = spitalfields$Age, HDImass = c(
        0.75, 0.95))
    Output
        coverage n_in      perc    CI_low     CI_up      p_value
      1     0.75  100 0.5555556 0.4797835 0.6294675 1.980641e-08
      2     0.95  151 0.8388889 0.7768829 0.8893679 2.874477e-08

