# SpeciesContributionToStability
This repository contains all R code and data used for the manuscript entitled **"Partitioning species contributions to ecological stability in disturbed communities"**

**Author list:** Charlotte Kunze, Dominik Bahlburg, Pablo Urrutia-Cordero, Maren Striebel, Egle Kelpsiene, Silke Langenheder, Ian Donohue & Helmut Hillebrand
      - Inquiries about the MS or the code to: charlotte.kunze@uol.de

**R version:** 4.3.2 (2022-03-10)

Platform: x86_64-apple-darwin17.0 (64-bit)

Running under: macOS 14.4.1

Zenodo: https://doi.org/10.5281/zenodo.11130744 


## How to work with this repository:

Download/Clone the project, set the working directory to the project location. 

The R environment (package versions + R version) used to develop the code is documented in the renv.lock-file contained in the repository. Package versions can either be looked up in the renv.lock-file or the entire session can be restored by installing the renv-package and running renv::restore() after opening the project. Depending on the R version on your local device, restoring the renv will return errors, then please install packages one by one as indicated in the **Main.R** script.

All scripts should then work as they are, if the required input-files are in the correct locations (look at the file paths specified in the scripts)



## To run the code: 

Please open the Rproject **SpeciesContributionToStability.Rproj** and the file **Main.R** and follow the instructions: Download data from Zenodo and store in folders *BEFD_createdData* and *SITES_Data* that are created in the **Main.R** file.

Output files are stored in the folder *“OutputSubmission”*, which will also be created in the **Main.R** file



