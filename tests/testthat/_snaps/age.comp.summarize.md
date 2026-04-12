# age.comp.summarize() produces correct output

    Code
      age.comp.summarize(bay.ta_compare, known_age = spitalfields$Age)
    Output
          Mean_estimated     Bias corrPearson       corr_p Residual_age_slope
      cor       58.36968 2.014126    0.541778 4.039942e-15          0.3994508
          Inaccuracy     RMSE TMNLP     CRPS Coverage HDI_Diff_median
      cor   14.04009 17.64114 4.461 9.072315 84.44444        45.23367
          HDI_Diff_quant_025 HDI_Diff_quant_975
      cor           32.31223           57.74559

# age.comp.summarize() produces correct output wit reduced dataset

    Code
      age.comp.summarize(bay.ta_compare, known_age = spitalfields_Age)
    Output
          Mean_estimated     Bias corrPearson       corr_p Residual_age_slope
      cor       58.79308 2.215934     0.54988 3.237091e-15          0.3951506
          Inaccuracy     RMSE  TMNLP     CRPS Coverage HDI_Diff_median
      cor   13.89124 17.51473 4.4464 9.015975 84.57143        45.14256
          HDI_Diff_quant_025 HDI_Diff_quant_975
      cor           32.19013           57.87153

