# age.comp.summarize() produces correct output

    Code
      age.comp.summarize(bay.ta_compare, known_age = spitalfields$Age)
    Output
          Mean_estimated     Bias corrPearson       corr_p Residual_age_slope
      cor       63.80984 7.454284   0.4522719 1.846485e-10          0.4610876
          Inaccuracy     RMSE  TMNLP     CRPS Coverage HDI_Diff_median
      cor   17.13358 21.37023 4.7709 9.469478 91.11111        58.39994
          HDI_Diff_quant_025 HDI_Diff_quant_975
      cor           39.68785           70.44868

# age.comp.summarize() produces correct output wit reduced dataset

    Code
      age.comp.summarize(bay.ta_compare, known_age = spitalfields_Age)
    Output
          Mean_estimated    Bias corrPearson       corr_p Residual_age_slope
      cor       64.27506 7.69792   0.4565118 2.160107e-10          0.4629466
          Inaccuracy     RMSE  TMNLP     CRPS Coverage HDI_Diff_median
      cor   17.08325 21.32045 4.7648 9.434711 91.42857        58.28629
          HDI_Diff_quant_025 HDI_Diff_quant_975
      cor           39.63137           70.47632

