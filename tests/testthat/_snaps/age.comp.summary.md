# age.comp.summarize() produces correct output

    Code
      age.comp.summary(bay.ta_compare, known_age = spitalfields$Age)
    Output
            Bias corrPearson       corr_p Residual_slope Inaccuracy     RMSE    TMNLP
      1 2.014126    0.541778 4.039942e-15      0.3994508   14.04009 17.64114 4.461029
            CRPS
      1 9.072315

# age.comp.summarize() produces correct output wit reduced dataset

    Code
      age.comp.summary(bay.ta_compare, known_age = spitalfields_Age)
    Output
            Bias corrPearson       corr_p Residual_slope Inaccuracy     RMSE    TMNLP
      1 2.215934     0.54988 3.237091e-15      0.3951506   13.89124 17.51473 4.446444
            CRPS
      1 9.015975

