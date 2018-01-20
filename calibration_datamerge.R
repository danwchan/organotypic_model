#' ---
#' title: "two replicates of cfu calibration"
#'output:
#'  html_document:
#'    toc: true
#'    theme: united
#' ---

#/*    Script to combine and process two replicates of the cfu calibration assay before merging with the old assay replicate*/
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
#/* YAML header, #' for drop into roxygen (#, ## headers), #+ for chunks (broken by spaces)*/

#/*set global knitr options*/
#+ knitr-options, include=FALSE
knitr::opts_chunk$set(warning = FALSE)

#' #Initial set up

#+ import-libraries, message=FALSE
#import libraries
require(stats)
require(gplots)
require(ggplot2)
require(scales)
require(grid)
require(reshape2)
require(dplyr)

#+ style, include=FALSE
#apply a graphpad prism-like theme
#The palette with grey:
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#The palette with black:
#cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
themeColourLargetext <- theme_classic() +
  theme(line = element_line(size = rel(1.5), lineend = 'round'),
        legend.position = c(0.65, 0.3),
        legend.text  = element_text(size = rel(2.1)),
        legend.key = element_blank(),
        legend.key.size = unit(5, 'char'),
        legend.text.align = 0,
        legend.title = element_blank(),
        axis.title = element_text(size = rel(2.1), face = 'bold'),
        axis.title.x = element_text(vjust = -0.5),
        axis.title.y = element_text(hjust = 0.5),
        axis.text = element_text(size = rel(1.9)),
        plot.title = element_blank())
themeTextClean <- theme_classic() +
  theme(line = element_line(size = rel(1.5), lineend = 'round'),
        legend.text = element_blank(),
        legend.key = element_blank(),
        axis.text = element_blank())# The palette with grey:


original_par <- par() #for resetting to original par after generating the plot in the null device
#setwd("/home/danwchan/Documents/Scientific_Publishing/Organotypic_model/Figures/160601_remake_old_graphs") #why is this bad again? set the working directory

#' ##Data processing:

#+ import-data
#' import the data from csv
raw_data1 <- read.csv('Data/calibration/160620_rep2_cfu_calibration.csv',comment.char = '#')
raw_data2 <- read.csv('Data/calibration/160705_rep3_cfu_calibration.csv',comment.char = '#')
raw_data <- rbind(raw_data1, raw_data2)
str(raw_data)

#+ process-data
working_data <- raw_data %>%
  transform(timepoint = as.factor(timepoint),
            date = as.factor(date)) %>% #make timepoint a factor
  mutate(id_merge = as.factor(paste(timepoint, sample_id, sep = '_')),
         exp_status = as.factor(ifelse(grepl("agar.*", as.character(notes)) == TRUE,
                                       "calibration",
                                       "")), #flag the experimental data from those used to QA the initial innoculum
         double_factor = as.factor(ifelse(grepl(".*double", as.character(notes)) == TRUE,
                                          "double",
                                          "")), #flag double drops, where the 3uL merges into one
         is_infected = as.factor(ifelse(grepl("infected.*", as.character(notes)) == TRUE,
                                   "infected",
                                   "")),
         cfu = ifelse(double_factor == 'double',
                      10^dilution * 50 * count / 2, 
                      10^dilution * 50 * count)) #calculate cfu, where double is normalized to the regular 3uL drop
str(working_data)

#' ##Some things about the experiment

#+ working-tables
calibration <- filter(working_data, 
                      exp_status == "calibration", 
                      notes != "agar_stock_d") #data of the 3uL drop applied to agar instead of the raft
agar_stock <- filter(working_data, notes == "agar_stock_d") #data of the cfu in initial innoculum
experiment <- filter(working_data, exp_status != "calibration", cfu != "NA") #experimental data

#+ calibration-plots
calibplot <- ggplot(calibration, aes(sample_id, count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = as.factor(date)), position = "jitter")
calibplot #looking at the test 3uL drops
calibplot2 <- ggplot(agar_stock, aes(sample_id, cfu)) +
  geom_point(aes(colour = as.factor(date), size = dilution)) +
  scale_size_continuous(range = c(3,10))
calibplot2 #looking at the cfu in the intial innoculum

#+ plot-data, message=FALSE, fig.width=12, fig.height=10
scatterplot <- ggplot(experiment, aes(sample_id, cfu)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = as.factor(date), size = double_factor, shape = is_infected), position = 'jitter') +
  #  scale_size_discrete(range = c(3,5)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  facet_grid(~timepoint)
scatterplot #this plot shows the data as it is reported to get a sense of the results before data processing

#+ cfu-delivered
#create the old estimates to merge with the new ones
old_estimates <- data_frame(cfu_per_drop = c(15, 210, 2700),
                            sample_id = c("dil1", "dil2", "dil3"),
                            timepoint = 0,
                            date = as.factor("131206"),
                            sample_strain = "LAC",
                            notes = "old esitmates manually entered",
                            count = "",
                            dilution = "",
                            id_merge = c("0_dil1", "0_dil2", "0_dil3"),
                            double_factor = "",
                            is_infected = "",
                            cfu = "",
                            exp_status = "calibration")
#calculate the mean cfu delivered
agar_stock %>%
  mutate(cfu_per_drop = 3*cfu/1000) %>%
  rbind(old_estimates) %>%
  group_by(sample_id) %>%
  summarise(mean(cfu_per_drop))

#+ clean-n-save, eval=FALSE
# remove the "infected data points since they inflate the means as seen in scatterplot
working_data_add <- filter(working_data, is_infected != "infected")
# save for import into 160531_cfu_calibration_old.R
save(working_data_add, file = "Data/cfu_calibration_merged.Rdata")
