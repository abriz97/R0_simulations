## R0-methods-comparison

Code and data to reproduce the analysis in: [Refining reproduction number estimates to account for unobserved generations of infection in emerging epidemics] (https://www.medrxiv.org/content/10.1101/2021.11.08.21266033v1)

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

-  `sims_BP` contains epidemics simulated via a Branching Process
-  `sims_SEIR` contains epidemics simulated using a stochastic SEIR model
-  `time_series_covid19_confirmed_global.csv` contains official COVID-19 incidence data from the Johns Hopkins University


#### `results`

-   `produce_simulations.html` contains information on how the simulated data were obtained.
-   `EpiEstim_comparison.html` compares the estimates otained from EpiEstim and our proposed adjusted method. 
-   `report_simulations.html` compares the performance of the Linear Exponential Growth method and EpiEstim with and without adjustment on simulated data.


### License

- The code in this repository is licensed under the Creative Commons Attribution 4.0  International ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)). 

- The data in the `data/time_series_covid19_confirmed_global.csv` file obtained from the  [Johns Hopkins University](https://github.com/CSSEGISandData/COVID-19.) are redistributed under the Creative Commons Attribution 4.0  International ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)).  Copyright Johns Hopkins University 2020.

### Funding

This work was supported by the Engineering and Physical Science Research Council through the EPSRC Centre for Doctoral Training in Modern Statistics and Machine Learning at Imperial and Oxford. We acknowledge the Abdul Latif Jameel Institute for Disease and Emergency Analytics, funded by Community Jameel and the MRC Centre for Global Infectious Disease Analysis (reference MR/R015600/1) jointly funded by the UK Medical Research Council (MRC) and the UK Department for International Development (DFID), under the MRC/DFID Concordat agreement and is also part of the EDCTP2 programme supported by the European Union. ID acknowledges research funding from a Sir Henry Dale Fellowship funded by the Royal Society and Wellcome Trust [gran 213494/Z/18/Z].

