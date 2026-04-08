# age.comp.summarize() produces correct output

    Code
      age.comp.summarize(bay.ta_compare, known_age = spitalfields$Age)
    Output
          Mean_estimated     Bias corrPearson       corr_p Residual_age_slope
      cor       71.56723 15.21167   0.5294805 2.138456e-14          0.5492014
          Inaccuracy     RMSE  TMNLP     CRPS Coverage HDI_Diff_median
      cor   17.70909 21.85109 4.7709 11.02673 83.88889        48.33585
          HDI_Diff_quant_025 HDI_Diff_quant_975
      cor           34.28944           65.00493

# age.comp.summarize() produces correct output wit reduced dataset

    Code
      age.comp.summarize(bay.ta_compare, known_age = spitalfields_Age)
    Output
          Mean_estimated     Bias corrPearson       corr_p Residual_age_slope
      cor       71.99159 15.41445    0.538056 1.605603e-14          0.5560766
          Inaccuracy    RMSE  TMNLP     CRPS Coverage HDI_Diff_median
      cor   17.66878 21.8376 4.7648 11.01945       84        47.85713
          HDI_Diff_quant_025 HDI_Diff_quant_975
      cor            34.2472           65.13761

