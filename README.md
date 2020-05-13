# art-muse-nicotine

Welcome to the code and raw data repo for our project (link to preprint or published paper will go here)!

## Setup

You can use the green "Clone or download" button to download a zipped file of this repository, or fork a copy to your own GitHub account and git clone it to your local machine. Once you get the unzipped folder on your machine, you can use the included .Rproj file to open an RStudio instance for the project.

The project was run in R 3.6.1, and we've used the [renv](https://rstudio.github.io/renv/articles/renv.html) package to catalog package dependencies. If you don't currently have renv installed on your machine, you can install the newest version from CRAN. Upon opening the R project, running `renv::restore()` will download all the necessary packages, at the versions we used, without conflicting with existing package versions you might have in your primary R setup. Once you've done that, you're ready to reproduce our analyses.

## Reproducing our analyses

To reproduce the analyses reported in our text, you'll need to:

1. Download our intermediate R files or reproduce them from the provided raw data
2. Knit R/report.Rmd to generate the results, embedded in the body text of our results section

To reproduce the results in the main text, open R/report.Rmd in RStudio and knit it to HTML. 

(We haven't provided a make recipe for the report file, as we've encountered OS-dependent issues with knitting .Rmds through make. The issue seems to stem from a failure to locate R's pandoc installation when using the Rscript bash command to run R code through the terminal in some Mac operating systems.)

### Downloading our R files

The intermediate R data are uploaded to Dryad (link will go here) for downloading. You will need to download these separately if you wish to use them, as they're a bit too large to provide via GitHub. Once you download the files, put them into the sub-folder [ignore/data_R](#data_r), and report.Rmd should behave just fine.

### Reproducing our R files from scratch

If you wish to recreate these data by running the code, you can use the recipes in the included Makefile to make the key ones, by running the following commands in a terminal with the working directory set to this folder:

```
make ignore/data_R/preplots_rstanarm.rda
make ignore/data_R/aprime_by_ppm.rda
```

**Please be aware that the rstanarm model objects take several hours in total to sample, and require your machine to be running the whole time.** Proceed with caution if you do choose this path!

## What's in the repo

The below descriptions are organized according to the folder structure.

### materials

#### task

In this folder, find experiment-running Matlab scripts and stimulus images necessary to run the "art gallery" attention task. This is the [same task](https://github.com/alylab/artmusePatient) that was run in Ruiz, Meager, Agarwal, & Aly (2020). The task files are identical to those in the other study (and identically named, so they have "patient" in the names per the first study). These files _should_ run out-of-the box _if you have the following Matlab setup already on your machine._

- Matlab 2017a or 2017b, should run on newer editions but unconfirmed
- [PsychToolbox 3.0.15](http://psychtoolbox.org/download.html)

We're currently unable to provide an environment file that will auto-download all the dependencies for you. However, the code info text file in the folder provides some additional instructions once you've got everything downloaded.

#### data

In this folder, find:

- raw .txt and .mat task data files, one each per participant per session
- a data dictionary explicating each of the column contents of those raw task data files
- a .csv of demographics for each participant, including session-specific expired breath carbon monoxide parts-per-million measurements

If you choose to reproduce our intermediate data files and results from scratch, running the code will draw from the data in this folder, so you don't need to download additional files.

### R

In this folder, find analysis code for our main results.

- report.Rmd: This file generates main text figures and numerical results, embedded in a local copy of the results text.
- explore.Rmd: This file contains our various explorations into the data, provided for transparency. However, it's a bit messy and has not been cleaned up as much as the rest of the repository. Beware!
- *.R: These files contain code to preprocess raw data and estimate various model objects. They output various .rda files that are expected in the .Rmd report files.

### ignore

If you've just downloaded the repo, this folder won't exist! We've added this folder to the .gitignore, as it stores big data that is too unwieldy for GitHub. You'll need to create this subfolder, and the two child folders indicated below, inside the main art-muse-nicotine directory, if you want to re-run the code.

#### data_R

If you download our intermediate files from Dryad, those should be placed here. Alternatively, R files you reproduce by re-running our code will automatically be saved here.

#### figures

This folder will contain any figures rendered out from report.Rmd. These are results figures from the main text, as well as alternative versions with different levels of grid lines in the back of the plots.
