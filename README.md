## R0-methods-comparison

Code and data to reproduce the analysis in:

### Code

The github directory is structured as follows:

- The `scripts` directory contains all the file necessary to reproduce the analyses described in the paper, as well as the main and supplementary figures.
- The `data` directory contains the SEIR simulations used in our analysis as well as the COVID-19 data obtained from the [Hopkins University TODO]( https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series)
- The `results` data contains html reports with code and plots.

#### `scripts`	

- `report_simulations.Rmd` contains the code to run the simulation experiments.
- `real_data.R` containst the code to run the analysis on COVID data.
- `methods.R` contains helper functions that fit different methods for the estimation of $R_0$.
- `produce_simulations.R` makes use the `simulation.R` script to generate simulated epidemic curves, save them in the data folder and visualise them

#### `data`

#### `results`



### License

- The code in this repository is licensed under the Creative Commons Attribution 4.0  International ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)) by TODO. 

- The data in the `data/time_series_covid19_confirmed_global.csv` file obtained from the  [Johns Hopkins University](https://github.com/CSSEGISandData/COVID-19.) are redistributed under the Creative Commons Attribution 4.0  International ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).  Copyright Johns Hopkins University 2020.

### Acknowledgement

### Funding

