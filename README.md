# SpeciesContributionToStability
This repository contains all R code and data used for the framework on measuring species net contribution to functional stability after perturbation


R version 4.1.3 (2022-03-10)

Platform: x86_64-apple-darwin17.0 (64-bit)

Running under: macOS 13.4


How to work with this repository:

There are two folders containing analyses and data creation for model simulations (BEFD_Analysis) and empirical data (SITES_Analysis).
To run the code, please download the data available in Zenodo and store them in BEFD_createdData. For empirical data create a folder 'SITES_Data' and download the raw data via the SITES data portal (Langenheder et al. 2020). 

Download/Clone the project, set the working directory to the project location. The R environment (package versions + R version) used to develop the code is documented in the renv.lock-file contained in the repository and in the README.docx. Package versions can either be looked up in the renv.lock-file or the entire session can be restored by installing the renv-package and running

renv::restore()
after opening the project.

All scripts should then work as they are when the required input-files are in the correct locations (look at the file paths specified in the scripts)
