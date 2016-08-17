# Purpose

Reproducible workflow used to process and analyze the data associated with the working paper entitled: Staphylococcal Smooth Biofilm Morphology on Organotypic Human Keratinocyte Culture is Dependent on agrA.

The data for this project will be made available upon publication. This repo is currently linked to an OSF project.

# File structure

This is my first time making a paper like this. And I think it helps to have a written explanation of how I have organized things in my mind.
If not for others to follow than for me to internalize the lessons learned in the process of compling this paper.

## / (root)

The makefile is placed in the root directory.

Analysis scripts written in R are also placed here. Scripts can be run with `rmarkdown::render` and `knitr::spin` to produce html documents to human reading OR they can be run as is to get intermediate files.
The scripts can do the following:

1. %_datawrangle can be called to produce a "wrangled" data objects for later visualization
2. %_stats_vis can be called to produce the final "atomic" figures saved as a part of html files which detail the process

## /Data

Wrangled data saved as .RData forms are saved here.
Raw data is stored here under subfolders which roughly correspond to the experimental logic.

## /Figures

Final figures as .pdf and working figures as .svg are saved here.
Intermediate files are high quality raster images which are then saved into pdf and removed.
This allows me to set up the layout of multi panel figures with embedded images in Inkscape which will automatically update them if changes are made in the atomic figure.

## /analysis_html

HTML reports from calling a full knitting of the R analysis scripts are saved here. These can be browsed instead of the source code for (perhaps?) easier understanding.