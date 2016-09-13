#' ---
#' title: "hla tidy data"
#'output:
#'  html_document:
#'    toc: true
#'    theme: united
#' ---
#' 
#/*    Data cleaning script for raft cfu assay, hla append section*/
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
#'going way back to replot and reanalyze data for the paper
#'
#'
#/*making a note book with knitr.spin is preferable because of the dynamic possible code debugging*/
#/* YAML header, #' for drop into Rmarkdown (#, ## headers), #+ for chunks (try not to disrupt code chunks with comments, place before)*/

#/*set global knitr options*/
#+ knitr-options, message=FALSE, echo=FALSE
knitr::opts_chunk$set(warning = FALSE, tidy = FALSE)

#' #Initial data and analysis set up
#'  
#+ import-libraries, message=FALSE
#import libraries
#plottting
require(gplots)
require(ggplot2)
require(ggthemes)
require(scales)
require(grid)
#data wrangl'n
require(purrr)
require(gdata)
require(reshape2)
require(data.table)
require(dplyr)
require(tidyr)

#+ Session-info
sessionInfo() #for reproducibility

#+ some-other-options
original_par <- par() #for resetting to original par after generating the plot in the null device

#' ##Data import and processing: 
#' 
#' copy and pasted with modification from my very first R scripts :P
#'

#+ importbind1
tablename <- "Data/hla/140516_WPH_test.csv"
source("old_import_chunk.R")
#' data from 
{{tablename}}
#' imported into melted_data, bind to working data...
str(melted_data)
working_data <- data.frame(melted_data, date = "140516")

#+ importbind2
tablename <- "Data/hla/140711_WPH_2nd.csv" #second table to import
source("old_import_chunk.R")
#' import2 chunk imported 
{{tablename}}
str(melted_data)
working_data <- rbind(working_data, data.frame(melted_data, date = "140711"))

#+ importbind3
tablename <- "Data/hla/141112_WHPacomp.csv" #3rd table to import
source("old_import_chunk.R")
#' data from 
{{tablename}}
#' imported into melted_data, bind to working data...
str(melted_data)
working_data <- rbind(working_data, data.frame(melted_data, date = "141112"))

#+ importbind4
tablename <- "Data/hla/150123_WHPcomp.csv" #last table to import
source("old_import_chunk.R")
#' data from 
{{tablename}}
#' imported into melted_data, bind to working data...
working_data <- rbind(working_data, data.frame(melted_data, date = "150123"))

#+ importbind5
tablename <- "Data/hla/150122_WPcomp.csv" #last table to import
source("old_import_chunk.R")
#' data from 
{{tablename}}
#' imported into melted_data, bind to working data...
working_data <- rbind(working_data, data.frame(melted_data, date = "150122"))

#+ process-data
working_data <- working_data %>%
  filter(sample_id == "WT" | sample_id == "hla" | sample_id == "PSM") %>% #vector used in these experiment is pTX isolated from strains from M.Otto
  filter(Time == 72 | Time == 120) %>%
  transform(timepoint = as.factor(Time)) %>% #make timepoint and CFU_delivered a factor
  mutate(id_merge = as.factor(paste(timepoint, sample_id, sep = '_'))) %>% #add an id column with all experimental variables
  na.omit() #remove NA's

#+ transform-data
max_cfu <- max(working_data$cfu)
min_cfu <- min(working_data$cfu)
norm_data <- working_data %>%
  mutate(cfu_log = log10(cfu), 
         cfu_scale = (cfu - min_cfu)/(max_cfu - min_cfu)) #use the log base 10 and scale to 0-1
norm_data <- norm_data[norm_data$cfu_log != -Inf,]

#+ expermental-variables
#getting more generalized
blocking_factor <- "conditions_block"  #the blocking factor is timepoint+mutants (operations are performed within blocks), required for ANOVA + orddom
#blocking_factor2 <- "timepoint" #back to just this in the reanalysis
test_factor1 <- "sample_id"
transformed_data <- c("cfu", "cfu_log", "cfu_scale")
alpha_level <- 0.05

#' #Quick Overview of Data
#' 

#+ overview-statistics
data_summary <- working_data %>%
  group_by(sample_id, timepoint) %>%
  summarize(mean = mean(cfu),
            sd = sd(cfu),
            ci_lower = t.test(cfu)$conf.int[1],
            ci_upper = t.test(cfu)$conf.int[2],
            n = length(cfu))
data_summary

#+ overview-plot
summary_plot <- ggplot(data_summary, aes(sample_id, mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = abs(ci_lower), ymax = abs(ci_upper))) + #negative CI are mostly because of n < 3
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  facet_grid(.~timepoint)
summary_plot

#+ first-plot,  message=FALSE, fig.width=12, fig.height=10
scatterplot <- ggplot(working_data, aes(sample_id, cfu)) +
  scale_color_brewer(palette = "YlOrBr") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = date, shape = sample_id), position = 'jitter', size = 3) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  facet_grid(.~timepoint)
scatterplot

#+ for latr merging, eval=FALSE
hla_tidy <- data.table(working_data)
#column normalization
hla_tidy[, Time := NULL][, ':=' (count = NA, dilution = NA, sample_strain = as.factor("LAC"), notes = "", exp_status = "", double_factor = "")]
#label normalization
hla_tidy[sample_id == "hla", sample_id := "hla_KO"][sample_id == "WT", sample_id := "wt"][sample_id == "PSM", sample_id := "3psm_KO"]
save(hla_tidy, file = "Data/hla_psm_tidy.RData")
