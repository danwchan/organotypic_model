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

date_index <- norm_data[sample_id %in% c("agrA_KO", "agrA_C123F_comp-20", "agrA_C123F_empty"),unique(date)] #find all experiments with the agr sample
norm_data <- norm_data[date %in% date_index][sample_id %in% c("wt", "agrA_KO", "agrA_C123F_comp-20", "agrA_C123F_empty")] #extract them from the data with thier matched wt

#' 
#' #Initial Visualize
#' 

#+ dirty-visualze
#+ data-summary
data_summary <- norm_data %>%
  group_by(sample_id) %>%
  summarise(
    mean = mean(cfu_log, na.rm = TRUE), # means comparison
    sdev = sd(cfu_log, na.rm = TRUE),
    ci_lower = t.test(cfu_log)$conf.int[1], #95% confidence intervals CANT DO IT CAUSE THE DATA?
    ci_upper = t.test(cfu_log)$conf.int[2])
data_summary

#+ summary-plot, message=FALSE, fig.width=12, fig.height=10
pos = position_dodge(width = 0.9)#for error bars to dodge dodging columns
data_summary_plot <- ggplot(data_summary, aes(sample_id, mean, ymin = ci_lower, ymax = ci_upper)) +
  geom_bar(aes(fill = sample_id), stat="identity", position = pos, width = 0.9) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  geom_errorbar(aes(fill = sample_id), width = 0.2, position = pos)
data_summary_plot

#' 
#' #Set values
#'
#' These variables affect the coming scripts.
#' blocking_factor is the column (data frame variable) in which all experimental conditions are tested in
#' test_factor1 is the column (data frame variable) which contains the experimental conditions
#' also reorder the factors to represent the presentation preferred order
#'

#+ set
transformed_data <- c("cfu", "cfu_log") #different transformations/normalizations
alpha_level <- 0.05
blocking_factor <- "timepoint"
test_factor1 <- "sample_id"
set_data <- transformed_data[[2]]
set_test <- "metric" #ordinal <- non-paramtric, cliff's D, metric <- parametric, cohen's D

norm_data$sample_id <- factor(norm_data$sample_id, levels(norm_data$sample_id)[c(2,1,4,3)])
norm_data$timepoint <- plyr::mapvalues(norm_data$timepoint, c("72", "120"), c("72 hours", "120 hours"))

