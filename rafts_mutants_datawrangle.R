#' ---
#' title: "merged raft cfu assays for publication, wrangle the data"
#'output:
#'  html_document:
#'    toc: true
#'    theme: united
#' ---

#/*    Data cleaning script for raft cfu assay*/
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
#import libraries for...
#plotting
require(gplots)
require(ggplot2)
require(ggthemes)
require(scales)
require(grid)
#modelling
require(car)
require(fitdistrplus)
require(lme4)
require(nlme)
#data wrangl'n
require(gdata)
require(reshape2)
require(data.table)
require(dplyr)
require(tidyr)
require(purrr)

sessionInfo() #for reproducibility

#+ functions, include=FALSE
#gets the AIC from a list of fitted data
get_AIC <- function(fit){
  AICs <- c()
  for (i in 1:length(fit)){
    AICs <- c(AICs, fit[[i]]$aic)
  }
  return(AICs)
}

#+ some-other-options
original_par <- par() #for resetting to original par after generating the plot in the null device

#' ##Data import and processing
#' 
#' import the premerged .csv file
#' add additional experiments as required
#' 

#+ import-data
#' import the data from csv
raw_data <- read.csv('Data/rafts_cfu_merged/160520_merged_KO.csv',comment.char = '#') #import the data from csv
str(raw_data)

#+ process-data
working_data <- raw_data %>%
  transform(timepoint = as.factor(timepoint)) %>% #make timepoint a factor
  mutate(id_merge = as.factor(paste(timepoint, sample_id, sep = '_')),
         exp_status = as.factor(ifelse(grepl("agar.*", as.character(notes)) == TRUE,
                                       "calibration",
                                       "")), #flag the experimental data from those used to QA the initial innoculum
         double_factor = as.factor(ifelse(grepl(".*double", as.character(notes)) == TRUE,
                                          "double",
                                          "")), #flag double drops, where the 3uL merges into one
         cfu = ifelse(double_factor == 'double',
                      10^dilution * 50 * count / 2, 
                      10^dilution * 50 * count)) #calculate cfu, where double is normalized to the regular 3uL drop
str(working_data)

#+ import-data2
#i import the data for the agr experiments (held seperately before)
raw_data <- read.csv('Data/rafts_cfu_merged/160520_agr_comp_merged.csv',comment.char = '#')
str(raw_data)

#+ process-data2
agr_data <- raw_data %>%
  transform(timepoint = as.factor(timepoint)) %>% #make timepoint a factor
  mutate(id_merge = as.factor(paste(timepoint, sample_id, sep = '_')),
         exp_status = as.factor(ifelse(grepl("agar.*", as.character(notes)) == TRUE,
                                       "calibration",
                                       "")), #flag the experimental data from those used to QA the initial innoculum
         double_factor = as.factor(ifelse(grepl(".*double", as.character(notes)) == TRUE,
                                          "double",
                                          "")), #flag double drops, where the 3uL merges into one
         cfu = ifelse(double_factor == 'double',
                      10^dilution * 50 * count / 2, 
                      10^dilution * 50 * count)) #calculate cfu, where double is normalized to the regular 3uL drop
str(agr_data)


#+ append
#append the agr data
working_data <- rbind(working_data, agr_data)

# append hla data from R.data file
load("Data/hla_tidy.RData")
working_data <- rbind(working_data, hla_tidy)
working_data$sample_id <- factor(working_data$sample_id)

#' ##Transform
#' some typical transformations which might prove useful

#+ transform-data
max_cfu <- max(working_data$cfu)
min_cfu <- min(working_data$cfu)

norm_summary <- working_data %>%
  filter(exp_status != "calibration") %>%
  filter(notes != "infected" & notes != "infected_double" & notes != "double_infected" & notes != "visible") %>%
  filter(sample_id != "agrA_C123F_comp-P2min") %>%
  filter(cfu != "NA") %>%
  mutate(cfu_log = log10(cfu), 
         cfu_scale = (cfu - min_cfu)/(max_cfu - min_cfu)) #use the log base 10 and scale to 0-1

#' ##Dirty group comparison
#' 
#' look at the the merged data sets
#' 
#+ da-plots, message=FALSE, fig.width=12, fig.height=10
summary1plot <- ggplot(norm_summary, aes(sample_id, cfu)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(aes(colour = date), position = "jitter") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  scale_y_log10()
summary1plot 

summary2plot <- ggplot(norm_data, aes(simple_id, OD)) +
  scale_shape_manual(values = c(1,16,7)) +
  scale_colour_brewer(palette = "YlOrRd") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(aes(shape = coating, colour = date), size = 3, position = "jitter") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  scale_y_log10()
summary2plot #OD by simple groupings

summary2plot <- ggplot(norm_data, aes(simple_id, OD)) +
  scale_shape_manual(values = c(1,16,7)) +
  scale_colour_brewer(palette = "YlOrRd") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(aes(shape = coating, colour = dilution), size = 3, position = "jitter") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  scale_y_log10()
summary2plot #OD by simple groupings, colored by dilution

#'
#' since there is some variability over time with the absolute OD of the values 
#' (also because I am including FBS and non FBS treated wells, and since I can't tell what the dilution factor plays into this)
#' I will be using the normalized values for all analyses.
#'
#+ da-plot-final, message=FALSE, fig.width=12, fig.height=10
summary3plot <- ggplot(norm_data, aes(simple_id, OD_adjusted)) +
  scale_shape_manual(values = c(1,16,7)) +
  scale_colour_brewer(palette = "YlOrRd") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(aes(shape = coating, colour = date), size = 3, position = "jitter") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  scale_y_log10()
summary3plot #adjusted OD by simple groupings

#' ###Some retrospective fixes
#' 
#' some deviation from 1 for wt is OK considering that there are sometimes two wt observations and drift between them...
#' 
#' but excessive deviation warrants some investigation.
#' It seems like the 160326 and 160324 experiments are driving this therefore the investiagtion will begin there
#' 
#' after this analysis was borne out, I made the changes by adding a drop column and removing the offending data from the larger set
#' 

#+ investigation, include=FALSE
unique(norm_data[simple_id == "wt_NA", .(date, plate, OD_backgroundsub, OD, wt_well_values, OD_adjusted)])
investigate <- working_data %>%
  filter(date == "160326" | date == "160324")
invest_platemap <- ggplot(investigate, aes(col, row, label = sample_id)) +
  geom_raster(aes(fill=OD)) +
  geom_text(fontface = "bold", size = 3, angle = -45) +
  facet_grid(plate~date)
invest_platemap
#'lots of variation in the empty wells
#'
#' of the plates in 160326 plate 1 has the highest disagreement between wt, then 2, 1
unique(norm_data[simple_id == "wt_NA" & date == "160326", .(date, plate, row, OD_backgroundsub, OD, wt_well_values, OD_adjusted)])
#' of the plates in 160324 plate 4 has the highest disagreement between wt, then 2, 1 also edge wells are responsible
unique(norm_data[simple_id == "wt_NA" & date == "160324", .(date, plate, row, OD_backgroundsub, OD, wt_well_values, OD_adjusted)])
#' in the end I think we will exclude the experiments of 160326 and pre average wt values from 160324

#+ effect-size
data_summary <- norm_data %>%
  group_by(simple_id) %>%
  summarise(
    mean = mean(OD_adjusted, na.rm = TRUE), # means comparison
    sdev = sd(OD_adjusted, na.rm = TRUE),
    ci_lower = t.test(OD_adjusted)$conf.int[1], #95% confidence intervals CANT DO IT CAUSE THE DATA?
    ci_upper = t.test(OD_adjusted)$conf.int[2])
data_summary

#+ effect-size-plot, message=FALSE, fig.width=12, fig.height=10
pos = position_dodge(width = 0.9)#for error bars to dodge dodging columns
data_summary_plot <- ggplot(data_summary, aes(simple_id, mean, ymin = ci_lower, ymax = ci_upper)) +
  geom_bar(aes(fill = simple_id), stat="identity", position = pos, width = 0.9) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  geom_errorbar(aes(fill = simple_id), width = 0.2, position = pos)
data_summary_plot

#' #Assumptions Testing
#' 
#' ##Explore distributions
#' 
#+ split-data, message=FALSE, results="hide"
OD_aggregate <- aggregate(cbind(OD, OD_adjusted, OD_log, OD_scale) ~ simple_id,
                          data = norm_data,
                          print) #the OD for each pariwise combinations, also required for Cliff's DELTA
xselect <- combn(OD_aggregate[["simple_id"]], 2)
names(OD_aggregate$OD_adjusted) <- OD_aggregate$simple_id

#+ explore-distributions, fig.width=12, fig.height=30
par(mfcol = c(5,3), mar = c(2,2,2,1))
distribution_explore <- lapply(OD_aggregate$OD_adjusted, function(x) descdist(x, boot = 1000))
par(original_par)

#+ fits
fit_norm <- lapply(OD_aggregate$OD_adjusted, function(z) fitdist(z, "norm", method = "mle"))
fit_lnorm <- lapply(OD_aggregate$OD_adjusted, function(z) fitdist(z, "lnorm", method = "mle"))
fit_exp <- lapply(OD_aggregate$OD_adjusted, function(z) fitdist(z, "exp", method = "mle"))
fit_gamma <- lapply(OD_aggregate$OD_adjusted, function(z) fitdist(z, "gamma", method = "mle"))
fit_weibull <- lapply(OD_aggregate$OD_adjusted, function(z) fitdist(z, "weibull", method = "mle"))
fit_unif <- lapply(OD_aggregate$OD_adjusted, function(z) fitdist(z, "unif", method = "mle"))
fits <- c(exp = fit_exp, gamma = fit_gamma, norm = fit_norm, lnorm = fit_lnorm, unif = fit_unif, weibull = fit_weibull) #fits using untransformed data are combined

#+ compare-AICs
AIC_extract <- list(names(fit_norm), c("exp", "gamma", "norm", "lnorm", "unif", "weibull"))
AIC_summary <- get_AIC(fits) %>%
  matrix(15,6, dimnames = AIC_extract)
AIC_summary

#+ set-up1
#the coordinates of the subset of the dat you want to look at to assess fit
#wt
index1 <- 15

#+ plot-fits
legend <- c("gamma", "normal", "lognormal", "weibull") #legend text
data_subset <- list(fit_gamma[[index1]], fit_norm[[index1]], fit_lnorm[[index1]], fit_weibull[[index1]]) #to compare the normal to lognormal and weibull
par(mfcol = c(2,2), mar = c(2,2,2,1))
denscomp(data_subset, legendtext = legend)
qqcomp(data_subset, legendtext = legend)
cdfcomp(data_subset, legendtext = legend)
ppcomp(data_subset, legendtext = legend)
par(original_par)

#+ set-up2
#the coordinates of the subset of the dat you want to look at to assess fit
#agrA_
index1 <- 1

#+ plot-fits2, ref.label="plot-fits"

#+ set-up3
#the coordinates of the subset of the dat you want to look at to assess fit
#agrA_pos_20
index1 <- 2

#+ plot-fits3, ref.label="plot-fits"

#+ set-up4
#the coordinates of the subset of the dat you want to look at to assess fit
#agrA_pos1_min
index1 <- 3

#+ plot-fits4, ref.label="plot-fits"

#+ set-up5
#the coordinates of the subset of the dat you want to look at to assess fit
#agrA_pos1_empty
index1 <- 4

#+ plot-fits5, ref.label="plot-fits"

#' ## Test Normality

#+ Q-Q-plot
qqnorm_data <- function(x){
  Q <- as.data.frame(qqnorm(x, plot = FALSE))
  names(Q) <- c("xq", substitute(x))
  Q
}
theoretical_qq <- norm_data %>%
  group_by(simple_id) %>%
  do(with(., qqnorm_data(OD_adjusted)))

ggplot(data = theoretical_qq, aes(x = xq, y = OD_adjusted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("Theoretical") +
  ylab("Sample") +
  facet_wrap(~simple_id)

#' null: the data are normally distributed
#+ shapiro-wilks
norm_data[, shapiro.test(.SD$OD_adjusted)$p.value, by = simple_id]

#' ##Homoskedasticity
#' 
#' the null hypotheises of these tests are: all group variances are equal

#+ homogeneity-of-variances
bartlett.test(OD_adjusted ~ simple_id, data=norm_data) # Bartlett Test of Homogeneity of Variances (parametric)
fligner.test(OD_adjusted ~ simple_id, data=norm_data) # Figner-Killeen Test of Homogeneity of Variances (nonparamatric)

#' #Output
#' 
#' save the outputs into the Data folder for visualization and stats

#+ save
save(norm_data, file="Data/crystal_violet_biofilm.RData")
