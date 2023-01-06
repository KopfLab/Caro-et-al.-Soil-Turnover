# Supporting Information

This repository contains all source code needed to reproduce the calculations and plots of the following manuscript:

Caro et al. (In Review)

# What can I do with this code?

In publishing this repository, our hope is that this code is useful to other members of the scientific community. This repository is released under a Creative Commons BY (CC-BY) license, which means that all code published here can be shared and adapted for any purposes so long as appropriate credit and citation of the original paper is given. See attribution section for details.

# How do I run this code?

1. Download and install R for your operating system.
2. Download and install RStudio for your operating system.
3. Download a zip file of this repository and decompress it in a directory of your choosing on your computer.
4. Navigate to the directory and open `d2o.Rproj` file to start Rstudio and load this project's files.
5. Open the script(s) you would like to run. Scripts are numbered in the order they should be executed e.g, 01, 02, 03. Duplicate numbers mean those scripts can be run in any order relative to each other. We recommend beginning on script `02_LH_-SIP_Quantification.Rmd`, where growth rate calculations relevent to the manuscript take place.
6. Ensure that you have all of the required libraries installed by inspecting the `Setup` chunks. In these scripts, we note the CRAN/GitHub version/release that was used. If any libraries fail to install, note the name of the library and attempt to manually install its most recent version via CRAN or GitHub.
7. To generate an HTML report, select File --> Knit from the menu.


# Scripts

- `00_dD_label_measurements.Rmd`: Analysis of D2O label strength measurements.
- `00_soil_report.Rmd`: Plotting of soil geochemistry data.
- `01_LH-SIP_Data_Reduction.Rmd` initial data wrangling and reduction of GC-FID and GC-P-IRMS data.
- `02_LH-SIP_Quantification.Rmd`: Stable Isotope Probing Quantification
- `02_bacdive_taxonomy.Rmd`:
- `03_data_viz.Rmd`: figure generation and plot output
- `99_labeling_calculations.Rmd`: explanation of labeling calculations
- `99_labeling_model.Rmd`: explanation of the labeling model

**NOTE:** `01_LH-SIP_Data_Reduction.Rmd` will not successfully run on your local computer because this script reduces a large amount of raw GC-IRMS and GC-FID data. Because GitHub does not support uploading this amount of data, we are only providing this script for reference and inspection. Please contact the authors for raw data requests. Scripts `02` and onwards all run off of the cache that is supplied in this repository.

## Folders
- `cache` contains cached IRMS data files and script outputs
- `archive` contains depracated scripts
- `data` raw data (empty, as GitHub does not support publishing this volume of data).
- `data_output` exported data for sharing
- `fig_output` figures for sharing
- `html_output` html files for sharing
- `libs` libraries
