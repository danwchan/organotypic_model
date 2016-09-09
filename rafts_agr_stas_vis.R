#' ---
#' title: "visualizationa and statistics for only the agr mutants"
#'output:
#'  html_document:
#'    toc: true
#'    theme: united
#' ---

#/*    Analysis script for statistics and visualization of organotypic raft cfu assay --agr subset*/
#/*    Copyright (C) 2016  Daniel Chan*/

#/*    This program is free software: you can redistribute it and/or modify*/
#/*    it under the terms of the GNU General Public License as published by*/
#/*    the Free Software Foundation, either version 3 of the License, or*/
#/*    (at your option) any later version.*/

#/*    This program is distributed in the hope that it will be useful,*/
#/*    but WITHOUT ANY WARRANTY; without even the implied warranty of*/
#/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the*/
#/*    GNU General Public License for more details.*/

#/*    You should have received a copy of the GNU General Public License*/
#/*    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/


#/*making a note book with knitr.spin is preferable because of the dynamic possible code debugging*/
#/* YAML header, #' for drop into Rmarkdown (#, ## headers), #+ for chunks (try not to disrupt code chunks with comments, place before)*/

#/*set global knitr options*/
#+ knitr-options, message=FALSE, echo=FALSE
knitr::opts_chunk$set(warning = FALSE, tidy = FALSE)

#' #Initial data and analysis set up
#'  
#+ import-libraries, message=FALSE
#import libraries
#statistical tests
require(stats)
require(dunn.test)
require(orddom)
#plotting
require(gplots)
require(ggplot2)
require(ggthemes)
require(scales)
require(grid)
#data wrangl'n
require(gdata)
require(reshape2)
require(data.table)
require(dplyr)
require(tidyr)
require(purrr)

#+ functions
source("sig_bars_generator.R")
source("numeric_to_label.R")

#+ style
source("visual_formatting.R")

#+ Session-info
sessionInfo() #for reproducibility

#+ some-other-options
original_par <- par() #for resetting to original par after generating the plot in the null device

#+ load-data
load("Data/merged_raft_cfu.RData")
