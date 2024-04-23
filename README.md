# SpeciesContributionToStability
This repository contains all R code and data used for the framework on measuring species net contribution to functional stability after perturbation

Zenodo: https://doi.org/10.5281/zenodo.10027062

R version 4.1.3 (2022-03-10)

Platform: x86_64-apple-darwin17.0 (64-bit)

Running under: macOS 14.4.1


How to work with this repository:

There are two folders containing analyses and data creation for model simulations “BEFD_Analysis” and empirical data “SITES_Analysis”. The empirical data are stored in folders “SITES_Data”. Output files are stored in the folder “OutputSubmission”, which will be created in the Main.R file

To run the code, please download the simulation data available in Zenodo and store them in BEFD_createdData. 

Download/Clone the project, set the working directory to the project location. The R environment (package versions + R version) used to develop the code is documented in the renv.lock-file contained in the repository and in the README on Zenodo. Package versions can either be looked up in the renv.lock-file or the entire session can be restored by installing the renv-package and running

renv::restore()
after opening the project.

All scripts should then work as they are when the required input-files are in the correct locations (look at the file paths specified in the scripts)
