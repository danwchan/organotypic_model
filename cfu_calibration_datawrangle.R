#' ---
#' title: "CFU calibration analysis, data wrangle"
#'output:
#'  html_document:
#'    toc: true
#'    theme: united
#' ---

#/*    Data cleaning script for cfu calibration assay with wild type bacteria*/
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
#modelling
require(fitdistrplus)
#statistical tests
require(stats)
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

#+ functions, include=FALSE
#gets the AIC from a list of fitted data
get_AIC <- function(fit){
  AICs <- c()
  for (i in 1:length(fit)){
    AICs <- c(AICs, fit[[i]]$aic)
  }
  return(AICs)
}

#+ Session-info
sessionInfo() #for reproducibility

#+ some-other-options
original_par <- par() #for resetting to original par after generating the plot in the null device
#setwd("/home/danwchan/Documents/Organotypic_Model/160601_remake_old_graphs") #why is this bad again? set the working directory to this scripts location (don't separate me away)

#' ##Data processing: 
#' 
#+ import-proccess-data
#' copy and pasted woth modigfication from my very first R scripts :P
eval(parse(text=readLines("Data/calibration/131206_Wildtype_varied_innoculum.csv", n=1))) #extract the information in the old testnames comment
expdata <- read.table("Data/calibration/131206_Wildtype_varied_innoculum.csv", header=TRUE, sep=',', skip=1) #read in the csv
timepoints <- levels(as.factor(expdata$Time))
datacols <- names(expdata)
#### old formatting concerns
#calculate the amount of cells per 500uL of homogenate
#dilution factor in tablename adjusts to 10uL
#columns should be labelled test(i) for the counts and test(i).D for the dilution
####
for (i in 1:length(testnames)) {
  count <- paste('test', as.character(i), sep='')
  dilution <- paste('test', as.character(i), '.D', sep='')
  expdata[[testnames[i]]] <- expdata[[count]] * expdata[[dilution]] * 50
}
working_data <- melt(expdata[c('Time', testnames)], id.vars='Time', na.rm=TRUE, variable.name = "CFU_delivered", value.name = "cfu") #reshape the data

#+ process-data
# new-process
working_data <- working_data %>%
  transform(timepoint = as.factor(Time),
            CFU_delivered = as.factor(CFU_delivered)) %>% #make timepoint and CFU_delivered a factor
  mutate(date = as.factor("131206")) %>% #add a column with date as a factor
  filter(cfu != 'NA') %>% #remove NA's
  select(date,timepoint, CFU_delivered, cfu)
str(working_data)

#+ incorporate-new
# rename the CFU_delivered column levels
working_data$CFU_delivered <- plyr::mapvalues(working_data$CFU_delivered, c("2700_CFU", "210_CFU" ,"15_CFU"), c("dil3", "dil2", "dil1"))
#load the
load("Data/cfu_calibration_merged.Rdata")
working_data <- working_data_add %>%
  transform(CFU_delivered = sample_id) %>%
  filter(cfu != 'NA') %>%
  select(date,timepoint, CFU_delivered, cfu) %>%
  rbind(working_data) %>%
  mutate(id_merge = as.factor(paste(timepoint, CFU_delivered, sep = '_'))) %>% #add a column with date as a factor
  transform(CFU_delivered = factor(CFU_delivered),
            timepoint = factor(timepoint))

#+ transform-data
max_cfu <- max(working_data$cfu)
min_cfu <- min(working_data$cfu)

norm_data <- working_data %>%
  mutate(cfu_log = log10(cfu), 
         cfu_scale = (cfu - min_cfu)/(max_cfu - min_cfu)) #use the log base 10 and scale to 0-1

#' #Quick Overview of Data
#' 

#+ overview-statistics
data_summary <- working_data %>%
  group_by(CFU_delivered, timepoint) %>%
  summarise(mean = mean(cfu),
            sd = sd(cfu),
            ci_lower = t.test(cfu)$conf.int[1],
            ci_upper = t.test(cfu)$conf.int[2])
summary_plot <- ggplot(data_summary, aes(CFU_delivered, mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = abs(ci_lower), ymax = abs(ci_upper))) + #this abs correction is wrong for the t0 15CFU condition but I don't know why the other negaitve CI are the way they are
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  facet_grid(~timepoint)
summary_plot

#+ first-plot,  message=FALSE, fig.width=12, fig.height=10
scatterplot <- ggplot(working_data, aes(CFU_delivered, cfu)) +
  scale_color_brewer(palette = "YlOrBr") +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = date),position = 'jitter') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  facet_grid(.~timepoint)
scatterplot

#' #Assumptions Testing
#' 
#' I'm using log transformed data because of the notion that log transformed data better captures the biological error
#' 
#' given that bacteria grow exponentially so error at an early timepoint is magnified exponentailly at a latter timepoint
#'
#' ## Test Normality

#+ Q-Q-plot
qqnorm_data <- function(x){
  Q <- as.data.frame(qqnorm(x, plot = FALSE))
  names(Q) <- c("xq", substitute(x))
  Q
}
theoretical_qq <- norm_data %>%
  group_by(timepoint, CFU_delivered) %>%
  do(with(., qqnorm_data(cfu_log)))

ggplot(data = theoretical_qq, aes(x = xq, y = cfu_log)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("Theoretical") +
  ylab("Sample") +
  facet_grid(timepoint~CFU_delivered)

#' null: the data are normally distributed
#+ shapiro-wilks
aggregate(cfu_log ~ id_merge, data = norm_data, function(x) ifelse(length(x) >= 3, shapiro.test(x)$p.value, NA)) #shapiro wilks test

#' ##Homoskedasticity
#' 
#' the null hypotheises of these tests are: all group variances are equal

#+ homogeneity-of-variances
bartlett.test(cfu_log ~ id_merge, data=norm_data) # Bartlett Test of Homogeneity of Variances (parametric)
fligner.test(cfu_log ~ id_merge, data=norm_data) # Figner-Killeen Test of Homogeneity of Variances (nonparamatric)

#+ save, eval=FALSE
save(norm_data, file="Data/cfu_calibration.RData")