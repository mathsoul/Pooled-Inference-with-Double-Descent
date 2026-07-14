## Scaled Dataset
The raw data are scaled using the same method as the mean absolute scaled error (MASE), following the approach documented in [Makridakis et al. (2020)](https://doi.org/10.1016/j.ijforecast.2019.04.014). For each frequency (e.g., hourly, daily), three files are provided: forecasts, true values (truth), and errors, where errors are defined as forecasts minus true values.

In additional to the scaled datasets, we also include the 'M4-info.csv', downloaded from [Here](https://github.com/Mcompetitions/M4-methods/tree/master/Dataset). This file contains the starting date of each time series, which is used to create pooled variables.

