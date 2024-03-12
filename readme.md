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

 - Matlab (version 2019a or later)
 - R and RStudio (3.6.x)
  Install the following R packages:

	packages <- c("data.table", "ggplot2")

	install.packages(packages)

- Obtain raw data files and put them into /Data/

2) Directory structure

/Data/ 			Stores raw and analyzed data
/Figures/ 			Stores illustrator files to assemble final figures.
/src/paper      Stores codes to generate figure panels and illustrations
/functions/     

3) How to run the project

Extract source data into the folder `/Data/`
Run the script `run_me_AllFigures.m` in Matlab, with the root directory as the working directory.
