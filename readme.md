## Project title
Sodium salicylate improves detection of amplitude-modulated sound in mice

------------
## Description

This is the code to generate Figure panels for the manuscript "Sodium salicylate improves detection of amplitude-modulated sound in mice" by Van den Berg et al, 2024. 

--------

## Authors

- Maurits M. van den Berg*, m.m.vandenberg@erasmusmc.nl

- Aaron B. Wong*, a.wong@erasmusmc.nl (maintainer)

- Ghais Houtak, ghais.houtak@adelaide.edu.au

- Ross S. Williamson, ross.s.williamson@pitt.edu

- J. Gerard G. Borst, g.borst@erasmusmc.nl


\* equal contributions

---
## Last updated

12 Mar 2024

---
## BUILD INSTRUCTIONS

1) Dependencies

 - Obtain raw data files and extract the folder /Data/ in the root directory of this repository
 - Matlab (version 2019a or later)
 - R and RStudio (3.6.x)
   - Install the following R packages ("here","lme4","nlme","plyr","multcomp") or run the following code in the R console:
	```
	packages <- c("here","lme4","nlme","plyr","multcomp")
	install.packages(packages)
	```

2) Directory structure

  - `/Data/`  Stores raw and analyzed data.
  - `/Figures/`   Stores illustrator files to assemble final figures.
  - `/src/paper` Stores codes to generate figure panels and illustrations.
  - `/src/statistics` Stores codes for statistics related to linear mixed models (Fig .
  - `/functions/` Stores miscellaneous functions needed to run the code

3) How to run the project

 - Extract the folder `/Data/` from source data into the root directory
 - Run the script `run_me_AllFigures.m` in Matlab, with the root directory as the working directory.
 - Linear mixed-effect model statistics: Run the R markdown file `src\statistics\StatisticsSalicylate_2022-09-15.rmd` in R-studio. 
