# BM-Spatial-SIR

Replication package for Bisin, Alberto, and Andrea Moro. "JUE Insight: Learning Epidemiology by Doing: The Empirical Implications of a Spatial SIR Model with Behavioral Responses," Journal of Urban Economics 127, January 2022

![alt='size comparison'](output/images/nc5-app-animation-sizecomp.gif)

Overview
--------

The code in this replication package generates the data and figures for the paper using Python and Stata. The replicator should expect the code to run for about 5-10 hours (including the code necessary to reproduce appendix data and figures). This README file does not contain documentation describing how to generate new simulations. Please refer to the code in class_spatialModels.py and its comments for details.

We kindly ask authors using or adapting this code in scientific publications to cite our paper:
Bisin, Alberto, and Andrea Moro. "JUE Insight: Learning Epidemiology by Doing: The Empirical Implications of a Spatial SIR Model with Behavioral Responses," Journal of Urban Economics 127, January 2022

  BibTeX entry:

  ```
  @article{bisin-moro-JUE-21,
    title = {JUE Insight: Learning Epidemiology by Doing: The Empirical Implications of a Spatial-SIR Model with Behavioral Responses},
    journal = {Journal of Urban Economics},
    pages = {103368},
    year = {2021},
    issn = {0094-1190},
    doi = {https://doi.org/10.1016/j.jue.2021.103368},
    url = {https://www.sciencedirect.com/science/article/pii/S0094119021000504},
    author = {Alberto Bisin and Andrea Moro}
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
(use Environment/requirements.txt if using a virtual environment)

- (optional) The python code was run using the IDE Spyder version 4.2.0

- Stata 14 (used only to generate calibration data -provided in the distribution- and to generate data for figures and table for Section 7)

### Memory and Runtime Requirements

Most simulations take minutes to run on a standard (2020) desktop. The simulation of the model with city twice the size of the benchmark model may take a 1-2 hours depending on model parameters. We run 20 replications of each simulation using the multiprocessing package running 7 simulations in parallel. Depending on hardware, runtime can be reduced by increasing the value of the processors variable.

Figures from saved simulation data are generated in minutes

#### Summary

Approximate time needed to reproduce the analyses on a standard 2020 quad-core desktop machine with 4 cores / 8 threads: 10-12 hours
The execution time can be shortened by reducing the number of replications of each simulation (variable nboots in each simulation) or increasing the number of processors (variable processors).

#### Details

The code was last run on a **MacBook Pro, 1,7 GHz Quad-Core Intel Core i7** running on **MacOS Big Sur 11.4**

The code assumes the root directory is the working directory, loads data from input/ and saves output in output/

Description of programs/code
----------------------------

Files marked with \* contain modules imported by other files and are not meant to be executed

- \* class_spatialModels.py contains the following classes

    1) spatialSAYDR. This class contains all methods necessary to simulate a 5-state (behavioral) spatial-SIR model, and generates object attributes containing contagion statistics. Simulation data is saved in object attributes.
    2) spSAYDR_randLoc (child of spatialSAYDR). Simulates spatial-SIR with agents moved in randomly drawn locations every day
    3) spSAYDR_hetDensity (child of spatialSAYDR). Simulates spatial-SIR with initial density of agents decreasing from center of city
    4) spSAYDR_behav (child of spatialSAYDR). Simulates spatial-SIR with behavioral responses
    5) spSAYDR_behav_local (child of spSAYDR_behav). Simulates spatial-SIR with behavioral responses only around the contagion circle

- \* class_SIRmodel.py

    contains the SIRmodel class, with all methods necessary to simulate a (behavioral) SIR model (without the spatial elements)

- \* class_averageStats.py contains the averageStats class used to compute average statistics for several replications of each simulation. The class accepts as input a list of objects of the spatialSAYDR class

- \* class_simul_policy.py contains class and code to generate policy simulations reported in the external appendix

- sim-paper.py contains all code used to run simulations of the spatial_SIR. All simulation data generated in sim_base is saved in gzipped pickle format in output/ The file contains lists of objects corresponding to the estimated model, each object is one replication of a model simulation

- fig-paper.py    contains all code used to generate the figures in the paper

- sim-appendix.py code to run simulations of spatial_SIR models necessary to generate appendix figures

- fig-appendix.py code to generate appendix figures

- input/ directory containing input data

- stata/ directory containing stata code to generate input data and estimation of simulated data

- output/ empty directory, storing generated simulation data

- output/images empty directory, storing generated figures

The random seed is set as spatialSAYDR attribute q_seed (row 33 of class_spatialModels.py) and passed as a parameter when the attribute is called in sim_base.py. When running multiple replications of the same model, the seed is changed accordingly
(see e.g. file sim_base.py, row 53)

### License for Code

This code is licensed under a Creative Commons/CC-BY-SA 4.0 (Attribution-ShareAlike 4.0 International) license. See https://creativecommons.org/licenses/by-sa/4.0/ for details.

We kindly ask academics using or adapting this code in scientific publications to cite our paper:
  Alberto Bisin and Andrea Moro. "JUE Insight: Learning Epidemiology by Doing: The Empirical Implications of a Spatial SIR Model with Behavioral Responses," Journal of Urban Economics, Forthcoming 2021

  BibTeX entry:

```
  @article{BISIN2021103368,
  title = {JUE Insight: Learning Epidemiology by Doing: The Empirical Implications of a Spatial-SIR Model with Behavioral Responses},
  journal = {Journal of Urban Economics},
  pages = {103368},
  year = {2021},
  issn = {0094-1190},
  doi = {https://doi.org/10.1016/j.jue.2021.103368},
  url = {https://www.sciencedirect.com/science/article/pii/S0094119021000504},
  author = {Alberto Bisin and Andrea Moro}
  }
```

Instructions to Replicators
---------------------------

If you are running the files interactively, the working directory should be the root 
of the project. 

To generate an appropriate Python virtual environment, a requirements.txt 
file is provided under Environment/ 

All figures are saved in output/images. Intermediated data is saved in output/

0. (Optional step) To download and regenerate calibration data
- Run stata/importdata.do

1. Edit the first line in doall.sh to indicate the location of your
stata executable 

2. (all at once) To simulate the model and generate figures, all code can be run
at once by executing from a shell (after appropriately modifying
the path to the stata executable, see step 1) from the root directory 
of the project.

```
bash doall.sh
```

2. (manual runs) Alternatively, run the following files can be run from the root
of the project in this order:

```
1) run python sim_paper.py
2) run stata/sim_estimate_dens.do
3) run stata/sim_estimate_policies_multitimes-pdate20.do
4) run stata/sim_estimate_policies_multitimes.do (this generates latex code for the last table)
5) run python fig-paper.py (generates all figures)
```

To generate appendix figures:

```
6) run python class_simul_policy.py
7) run python fig-appendix.py
```


List of figures
---------------------------

| Figure # | Program      | Line Number | Output file
|----------|--------------|-------------|-----------------------------------------------
| 1        | fig-paper.py | 76          | output/images/nc5-baseline_pos-day10.png
| 1        | fig-paper.py | 76          | output/images/nc5-baseline_pos-day20.png
| 1        | fig-paper.py | 76          | output/images/nc5-baseline_pos-day30.png
| 2        | fig-paper.py | 234         | output/images/nc5-density_contagion2.pdf
| 3        | fig-paper.py | 76          | output/images/nc5-baseline_pos-day10.png
| 3        | fig-paper.py | 76          | output/images/nc5-baseline_pos-day20.png
| 3        | fig-paper.py | 76          | output/images/nc5-baseline_pos-day30.png
| 3        | fig-paper.py | 96          | output/images/nc5-randomcluster_pos-day10.png
| 3        | fig-paper.py | 96          | output/images/nc5-randomcluster_pos-day20.png
| 3        | fig-paper.py | 96          | output/images/nc5-randomcluster_pos-day30.png
| 4        | fig-paper.py | 307         | output/images/nc5-short-random-rates.pdf
| 5        | fig-paper.py | 405         | output/images/nc5-SIR-citysize-rates.pdf
| 6        | fig-paper.py | 496         | output/images/nc5-short-density_contagion1.pdf
| 7        | fig-paper.py | 619         | output/images/nc5-short-3densities.pdf
| 8        | fig-paper.py | 657         | output/images/nc5-hetdens1.pdf
| 8        | fig-paper.py | 129         | output/images/nc5-hetdens_pos-day0.png
| 9        | fig-paper.py | 76          | output/images/nc5-baseline_pos-day10.png
| 9        | fig-paper.py | 76          | output/images/nc5-baseline_pos-day20.png
| 9        | fig-paper.py | 76          | output/images/nc5-baseline_pos-day30.png
| 9        | fig-paper.py | 117         | output/images/nc5-nomove_pos-day10.png
| 9        | fig-paper.py | 117         | output/images/nc5-nomove_pos-day20.png
| 9        | fig-paper.py | 117         | output/images/nc5-nomove_pos-day30.png
| 9        | fig-paper.py | 117         | output/images/nc5-nomove_pos-day50.png
| 9        | fig-paper.py | 117         | output/images/nc5-nomove_pos-day150.png
| 9        | fig-paper.py | 117         | output/images/nc5-nomove_pos-day250.png
| 10       | fig-paper.py | 774         | output/images/nc5-short-nomovement-rateslarge.pdf
| 11       | fig-paper.py | 807         | output/images/nc5-SIR_beh_responses.pdf
| 12       | fig-paper.py | 891         | output/images/nc5-SIR_beh.pdf
| 13       | fig-paper.py | 960         | output/images/nc5-SIR_beh_local.pdf
| 14       | fig-paper.py | 989         | output/images/nc5-est_densitybetas.pdf
| 15       | fig-paper.py | 1138        | output/images/nc5-estimatedbeta_policies.pdf
| 15       | fig-paper.py | 1287        | output/images/nc5-estimatedbeta_beh-policies.pdf
| 16       | fig-paper.py | 1339        | output/images/nc5-prediction_nobias.pdf
| 17       | fig-paper.py | 1384        | output/images/nc5-prediction_withbias.pdf

List of Tables
---------------------------
| Table # | Program                           | Line Number | Output file
|---------|-----------------------------------|-------------|-----------------------------------------------
| 1       |                                   |             | manually generated
| 2       |sim_estimate_policies_multitimes.do| 160         | output/nc5-estimates-multitimes

Statistics reported in Section 4
---------------------------
| Program                           | Line Number | Output file
|-----------------------------------|-------------|-----------------------------------------------
| fig-paper.py                      | 142-149     | output/nc5-stats_on_section_4.txt


## Acknowledgements

* To generate this README we followed the template README for social science replication packages, available at https://social-science-data-editors.github.io/template_README/

* To generate a list of figure numbers with corresponding figure filenames from the LaTeX source:
https://stackoverflow.com/questions/66551105/lazy-quantifier-for-an-exact-word-in-large-text-with-newlines/66555601#66555601
Next, to find the code row numbers where each filename is generated (needs some formatting of the file generates by the previous step) run the following on a shell
while read line; do grep -n "$line" *.py; done < figs.txt
