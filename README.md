# BM-Spatial-SIR

Replication package for Alberto Bisin and Andrea Moro, "Learning Epidemiology by Doing: The Empirical Implications of a Spatial SIR Model with Behavioral Responses," NBER Working Paper 27590, February 2021

![alt='size comparison'](output/images/nc5-app-animation-sizecomp.gif)

Overview
--------

The code in this replication package generates the data and figures for the paper using Python and Stata. The replicator should expect the code to run for about 5-10 hours (including the code necessary to reproduce appendix data and figures). This README file does not contain documentation describing how to generate new simulations. Please refer to the code in class_spatialModels.py and its comments for details.

We kindly ask authors using or adapting this code in scientific publications to cite our paper:
  Alberto Bisin and Andrea Moro. "Learning Epidemiology by Doing: The Empirical Implications of a Spatial SIR Model with Behavioral Responses," NBER Working Paper 27590, February 2021

  BibTeX entry:

  ```
  @misc{bisinmoro2020geo,
    title = "Learning Epidemiology by Doing: The Empirical Implications of a Spatial-SIR Model with Behavioral Responses",
    author = "Bisin, Alberto and Moro, Andrea", institution = "National Bureau of Economic Research",
    type = "Working Paper",
    series = "Working Paper Series",
    number = "27590", year = "2021",
    month = "February",
    doi = {10.3386/w27590},
    URL = "http://www.nber.org/papers/w27590",
    howpublished = "NBER Working Paper Series \# 27590"
    }
```

Data Availability and Provenance Statements
-------------------------------------------

- Data used for the model calibration is downloaded from the Italian government github web site https://github.com/pcm-dpc/COVID-19 and analyzed using Stata. The code in stata/importdata.do downloads the raw data and creates a moving average of growth rates of infections in the Lombardy region, saved in input/drlombardia5.txt. A copy of the original source data file downloaded on 12/16/2020 is saved in csv format in stata/rawdata.csv

### Statement about Rights

- We certify that the authors of the manuscript have legitimate access to and permission to use the data used in this manuscript. The data is licensed under the CreativeCommons CC-BY-4.0 (Attribution 4.0) International license, https://creativecommons.org/licenses/by/4.0/deed.it


Dataset list
------------

All files are provided

| Data file               | Source    | Notes                        |
|-------------------------|-----------|------------------------------|
| `stata/rawdata.csv`     | generated | see stata/importdata.do      |
| `stata/rawdata.dta`     | generated | stata version of rawdata.csv |
| `input/drlombardia.txt`| generated | used for calibration         |

Computational requirements
---------------------------

This code has been run on a MacBook Pro, 1,7 GHz Quad-Core Intel Core i7 running on MacOS Catalina 10.15.6.

### Software Requirements

- Python 3.8.5
  - `numpy` 1.19.2
  - `matplotlib` 3.3.2
  - `pandas` 1.1.3

- (optional) The python code was run using the IDE Spyder version 4.2.0

- (optional) The code used to download and generate calibration data (provided in the distribution) was run with Stata 14.

### Memory and Runtime Requirements

Most simulations take minutes to run on a standard (2020) desktop. The simulation of the model with city twice the size of the benchmark model may take a 1-2 hours depending on model parameters. We run 20 replications of each simulation using the multiprocessing package running 4 to 6 simulations in parallel. Depending on hardware, runtime can be reduced by increasing the value of the nprocs variable.

Simulation data is saved in pickle format under output/

Figures from such data is generated in seconds

#### Summary

Approximate time needed to reproduce the analyses on a standard 2020 desktop machine: 5-10 hours
The execution time can be shortened by reducing the number of replications of each simulation (variable nboots in each simulation). Results may vary slightly

#### Details

The code was last run on a **MacBook Pro, 1,7 GHz Quad-Core Intel Core i7** running on **MacOS Catalina 10.15.6**

The code assumes the root directory is the working directory, loads data from input/ and saves output in output/

Description of programs/code
----------------------------

Files marked with \* contain modules imported by other files and are not meant to be executed

- \* class_spatialModels.py contains the following classes

    1) spatialSAYDR. This class contains all methods necessary to simulate a 5-state (behavioral) spatial-SIR model, and generates object attributes containing contagion statistics. Simulation data is saved in object attributes.
    2) spSAYDR_randLoc (child of spatialSAYDR). Simulates spatial-SIR with agents moved in randomly drawn locations every day
    3) spSAYDR_hetDensity (child of spatialSAYDR). Simulates spatial-SIR with initial density of agents decreasing from center of city
    4) spSAYDR_behav (child of spatialSAYDR). Simulates spatial-SIR with behavioral responses

- \* class_SIRmodel.py

    contains the SIRmodel class, with all methods necessary to simulate a standard (behavioral) SIR model (without the spatial elements)

- \* class_averageStats.py contains the averageStats class used to compute average statistics for several replications of each simulation. The class accepts as input a list of objects of the spatialSAYDR class

- \* class_simul_policy.py contains class and code to generate policy simulations reported in the external appendix

- sim-paper.py contains all code used to run simulations of the spatial_SIR. All simulation data generated in sim_base is saved in gzipped pickle format in output/ The file contains lists of objects corresponding to the estimated model, each object is one replication of a model simulation

- fig-paper.py contains all code used to generate the figures in the paper

- sim-appendix.py code to run simulations of spatial_SIR models necessary to generate appendix figures

- fig-appendix.py code to generate appendix figures

- input/ directory containing input data

- stata/ directory containing stata code to generate input data

- output/ empty directory, storing generated simulation data

- output/images empty directory, storing generated figures

The random seed is set as spatialSAYDR attribute q_seed (row 33 of class_spatialModels.py) and passed as a parameter when the attribute is called in sim_base.py. When running multiple replications of the same model, the seed is changed accordingly
(see e.g. file sim_base.py, row 53)

### License for Code

This code is licensed under a Creative Commons/CC-BY-SA 4.0 (Attribution-ShareAlike 4.0 International) license. See LICENSE.txt for details.

We kindly ask academics using or adapting this code in scientific publications to cite our paper:
  Alberto Bisin and Andrea Moro. "Learning Epidemiology by Doing: The Empirical Implications of a Spatial SIR Model with Behavioral Responses," NBER Working Paper 27590, February 2021

  BibTeX entry:

  ```
  @misc{bisinmoro2020geo,
    title = "Learning Epidemiology by Doing: The Empirical Implications of a Spatial-SIR Model with Behavioral Responses",
    author = "Bisin, Alberto and Moro, Andrea", institution = "National Bureau of Economic Research",
    type = "Working Paper",
    series = "Working Paper Series",
    number = "27590", year = "2021",
    month = "February",
    doi = {10.3386/w27590},
    URL = "http://www.nber.org/papers/w27590",
    howpublished = "NBER Working Paper Series \# 27590"
    }
```

Instructions to Replicators
---------------------------

(Optional step) To download and regenerate calibration data (optional)
- Run stata/importdata.do

To simulate the model and generate figures, except Section 7:

1) Run sim_base.py
2) Run fig-shortpaper.py

To generate figures for Section 7
(apologies... this has not been cleaned up)
3) Run sim_estimation.py up to row 254
4) Run sim_estimate.do in Stata
5) Run sim_estimation.py up to row 847
6) Run sim_estimate_policies_multitimes-pdate20.do in Stata
7) Run sim_estimation.py
8) Run sim_estimate_policies_multitimes.do (this generates figures in the last table)

(To generate appendix figures)

9) Run fig-appendix.py
10) Run class_simul_policy.py



All figures are saved in output/images

List of figures
---------------------------

| Figure/Table #    | Program                  | Line Number | Output file                      | Note                            |
|-------------------|--------------------------|-------------|----------------------------------|---------------------------------|

please refer to latex files under \latex for figure names and search code for that name
(note: bib and style files are missing therefore these files are reproduced here only fore reference and not for compilation)

---

## Acknowledgements

In generating this README we followed the template README for social science replication packages, available at https://social-science-data-editors.github.io/template_README/
