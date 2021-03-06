# ICBM-OCEAN 
Institute for Chemistry and Biology of the Marine Environment - **O**ldenburg **C**omplex molecular mixture **E**valuation & **A**nalysis

# Supplementary Code Modules

These functions provide details about several calculations within ICBM-OCEAN, which are not covered by already published "R" packages mentioned inside the ICBM-OCEAN publication.

These functions include:

*fastjoin.R:* 

Merge masses considered to be equal from different samples (spectra) into a large table. Calculate the (intensity weighted) mean mass as an improved estimator of the true mass. The algorithm used here is fast, but can lead to chaining of masses, when the precision of the mass spectrometer is not high enough or when masses of equal origin can have a higher error than the distance to masses of different origin from another sample. 

*precisejoin.R:* 

Merge masses considered to be equal from different samples (spectra) into a large table. Calculate the (intensity weighted) mean mass as an improved estimator of the true mass. The algorithm used here is magnitudes slower than fast join, but can compensate chaining of masses better than *fastjoin.R*, when the precision of the mass spectrometer is not high enough or when masses of equal origin can have a higher error than the distance to masses of different origin from another sample. 
Always use fastjoin, if your instrument has a high enough precision.

*ResPow_outlier.R:*

Estimates and eliminates side peaks that have a suspiciously high resolution power. Calculations are based on a median regression of mass as explaratory variable and resolution power as response. The slope of the kernel density estimate of the resulting residuals indicates the point where peaks are considered to be outliers. 

*isotope_ratio_deviance_function.R:*

Calculates isotope ratio deviances e.g. Delta-C-13. using the ratio of observed intensities and expected intensities, where expected intensities are calculated based on binomial probabilities. A violin plot function is integrated to visualize results. 

*create_homologous_series_network.R:*

This function creates homologous series networks from molecular formulae. For this it uses the formulae instead of Kendrick mass defect calculations. Different links (e.g. CH2, O) can be selected. The network can be visualized interactively. 



