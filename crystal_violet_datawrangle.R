#' ---
#' title: "merged biofilm assays for publication, wrangle the data"
#'output:
#'  html_document:
#'    toc: true
#'    theme: united
#' ---

#/*    Data cleaning script for a crystal violet biofilm assay*/
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

#' ##Data import
#' 
#' import a variety of experimental done with non-standard annoation.
#' they are all done according to the protocol in the methods, but the conditions are frequently changed
#' 

#+ import-data
#' import the data from csv
raw_data1 <- read.csv('Data/crystal_violet_biofilm/160524_cv_merge.csv',comment.char = '#')
raw_data2 <- read.csv('Data/crystal_violet_biofilm/160316_cv_mutants_full1.csv',comment.char = '#')
raw_data2[,"date"] <- rep(160316, 240) # add dates for those instances it does not exist
raw_data3 <- read.csv('Data/crystal_violet_biofilm/160324_cv_mutants_full2.csv', comment.char = '#')
raw_data3[,"date"] <- rep(160324, 528)
raw_data3 <- filter(raw_data3, plate != 'p8', plate != 'p7', plate != 'p6', plate != 'p5') #those plates are diluted repeats
raw_data4 <- read.csv('Data/crystal_violet_biofilm/160326_cv_mutants_full3.csv', comment.char = '#')
raw_data4[,"date"] <- rep(160326, 264)
raw_data5 <- read.csv('Data/crystal_violet_biofilm/160505_cv_o2.csv', comment.char = '#')
raw_data5[,"date"] <- rep(160505, 264)
raw_data6 <- read.csv('Data/crystal_violet_biofilm/160414_cv_mutants_method.csv', comment.char = '#')
raw_data6[,"date"] <- rep(160414, 144)
raw_data7 <- read.csv('Data/crystal_violet_biofilm/160416_cv_mutants_method2.csv', comment.char = '#')
raw_data7[,"date"] <- rep(160416, 198)

#' ##Data processing
#' 
#' the raw data is melted and reshaped standardizing the column names and entries (somewhat) to make it easier to work with
#' 

#+ melt-merge
data_list <- list(raw_data1, raw_data2, raw_data3, raw_data4, raw_data5, raw_data6, raw_data7)
melt_pool <- melt(data_list[1], id.vars = c('sample_strain', 'sample_id', 'location', 'date'))
for(i in 2:length(data_list)){
  m <- melt(data_list[i], id.vars = c('sample_strain', 'sample_id', 'location', 'date'))
  melt_pool <- rbind(melt_pool, m)
  #  str(m) for debug
} #merge molten data
cast_data <- dcast(melt_pool, location + sample_strain + sample_id + date ~ variable, fill = "not_reported") #recast
cast_data <- data.table(cast_data) #first data.table useage!
cast_data[sample_strain == "", sample_strain := "LAC"] #update those rows with sample_strain=="" to "LAC"
cast_data[media == "not_reported", media := "TSB"] #absent media is TSB
cast_data[coating == "not_reported" | coating == "", coating := "none"] #absent coating is none
cast_data[O2 == "not_reported", O2 := "aerobic"] #absent O2 is TSB
cast_data[date == "160326", drop := 'true'] #remove these based on their lack of useable data (found when looking at the raw_platemap)
#EDIT after investigation remove all of 160326
cast_data[plate == "p3" & date == "160324", drop := 'true']
cast_data[sample_strain == "LAC" & sample_id == "empty", sample_strain := "empty"] #fix up this improper labelling
cast_data[drop == "not_reported", drop := "F"]
str(cast_data)

#' 
#' the data is now flagged and the types are changed to work better with downstream exploration and visulaization
#' 

#+ process-data
working_data <- cast_data %>%
  transform(date = as.factor(date),
            plate = as.factor(plate),
            coating = as.factor(coating),
            OD = as.numeric(OD)) %>% #convert date, plate to a factor and OD into numeric, merge the mutants together
  filter(sample_strain == 'LAC' | sample_strain == "empty",
         media == "TSB" | media == "0.5_glu",
         O2 == "aerobic") %>% #extract only the data for which there will be comparative rafts, also exclude media comparisons
  mutate(row = as.factor(ifelse(grepl("[a-g][0-9]{1,2}", as.character(well)) == TRUE,
                                sub("([a-g])[0-9]{1,2}", "\\1", as.character(well)),
                                "error")),
         col = as.factor(ifelse(grepl("[a-g][0-9]{1,2}", as.character(well)) == TRUE,
                                sub("[a-g]([0-9]{1,2})", "\\1", as.character(well)),
                                "error")), #process and flag data with positional information
         drop = as.logical(drop)) %>%
  separate(sample_id, c("genetic_background", "exogenous_genetic"), sep = "_[KOermCF123]{0,5}_?",
           extra = "merge", remove = FALSE) #separate the sample_id into components
str(working_data)

#'
#' from the separate function it's clear that there needs to be some reworking of the categories
#' queries are made and categories simplified/merged
#'

#+ checking-data, results="hide"
setkey(working_data, "exogenous_genetic")
unique(working_data[,.(exogenous_genetic)]) # we need to clean up this category
unique(working_data[NA_character_,.(exogenous_genetic, genetic_background, sample_strain, sample_id)]) #what are the entries which have NA?
unique(working_data[exogenous_genetic == "",.(exogenous_genetic, genetic_background, sample_strain, sample_id)])
unique(working_data["20",.(exogenous_genetic, genetic_background, sample_strain, sample_id)]) # 1st working through the list
working_data["20", exogenous_genetic := "pos1_agrA_20"]
unique(working_data[exogenous_genetic == "pos1_empty",.(exogenous_genetic, genetic_background, sample_strain, sample_id)])
unique(working_data[exogenous_genetic == "empty",.(exogenous_genetic, genetic_background, sample_strain, sample_id)])
working_data[exogenous_genetic == "empty", exogenous_genetic := "pos1_empty"]
unique(working_data[exogenous_genetic == "37" | exogenous_genetic == "77",.(exogenous_genetic, genetic_background, sample_strain, sample_id)])
working_data[exogenous_genetic == "37" | exogenous_genetic == "77", exogenous_genetic := ""]
unique(working_data[exogenous_genetic == "37_empty" | exogenous_genetic == "77_empty",.(exogenous_genetic, genetic_background, sample_strain, sample_id)])
working_data[exogenous_genetic == "37_empty" | exogenous_genetic == "77_empty", exogenous_genetic := "pos1_empty"]
unique(working_data[exogenous_genetic == "37_comp",.(exogenous_genetic, genetic_background, sample_strain, sample_id)])
working_data[exogenous_genetic == "37_comp", exogenous_genetic := "pos1_atl"]
unique(working_data[exogenous_genetic == "comp20" | exogenous_genetic == "77_pos1_20",.(exogenous_genetic, genetic_background, sample_strain, sample_id)])
working_data[exogenous_genetic == "comp20" | exogenous_genetic == "77_pos1_20", exogenous_genetic := "pos1_agrA_20"]
unique(working_data[exogenous_genetic == "comp",.(exogenous_genetic, genetic_background, sample_strain, sample_id)])
working_data[exogenous_genetic == "comp" & genetic_background == "atl", exogenous_genetic := "pos1_atl"
             ][exogenous_genetic == "comp" & genetic_background == "ica", exogenous_genetic := "pos1_ica"
               ][exogenous_genetic == "comp" & genetic_background == "srtA", exogenous_genetic := "pos1_srtA"]
unique(working_data[,.(genetic_background)]) # check the other separated column
unique(working_data[genetic_background == "",.(exogenous_genetic, genetic_background, sample_strain, sample_id)])
working_data[, simple_id := as.factor(paste(genetic_background, exogenous_genetic, sep = "_"))]

#' 
#' some last final cleaning
#' 

#+ remove-160324-col6
working_data[date == "160324" & col == "6", drop := TRUE]
working_data <- filter(working_data, drop == FALSE)

#' ##Quick visulaization
#' 
#' the plotted raw OD's in the style of Tecan output
#' 

#+ raw-OD-map, fig.width=15, fig.height=20
raw_platemap <- ggplot(working_data, aes(col, row, label = sample_id)) +
  geom_raster(aes(fill=OD)) +
  geom_text(aes(colour = sample_strain),fontface = "bold", size = 3, angle = -45) +
  facet_grid(plate~date)
raw_platemap

#' ##Normalize

#' because the "biofilms" are sensitive to washing I am normalizing them to the "empty" (background) well and "wt" (proportional constant) well in each row
#' if there are more than one wt samples of the same strain: the "wt" wells are averaged per row

#add the OD_adjusted column which takes the "empty" in the same row and subtracts it from the rest of the values then divde that by the adjusted "wt" value in the same row
#to test the function: x <- filter(working_data, sample_id =="wt")

#+ normalize-data
norm_data <- working_data[sample_strain == "empty", empty_well_values := OD
                          ][, OD_backgroundsub := OD - mean(empty_well_values, na.rm = TRUE), keyby = .(date, plate, row)
                            ][sample_id == "wt", wt_well_values := OD_backgroundsub
                              ][, OD_adjusted := OD_backgroundsub / mean(wt_well_values, na.rm = TRUE), keyby = .(date,plate,row)
                                ][, id_merge := as.factor(paste(date, genetic_background, "KO", exogenous_genetic, sep = '_')) ] # and make a merged id

#+ transform-data
norm_data <- norm_data[!OD_adjusted < 0,] #some anomolusly low OD_adjusted values
max_OD <- max(norm_data$OD_adjusted)
min_OD <- min(norm_data$OD_adjusted)
norm_data <- norm_data %>%
  mutate(OD_log = ifelse(OD_adjusted > 0, log10(OD_adjusted), 0),
         OD_scale = (OD_adjusted - min_OD)/(max_OD - min_OD))%>% #use the log base 10 and scale to 0-1
  filter(sample_strain != "empty", drop == FALSE)
norm_data <- norm_data[!OD_adjusted < 0] #get rid of some wierd under 0 points (what are they exactly?)

#+ norm-OD-map, fig.width=15, fig.height=20
norm_platemap <- ggplot(norm_data, aes(col, row, label = sample_id)) +
  geom_raster(aes(fill=OD_adjusted)) +
  geom_text(fontface = "bold", size = 3, angle = -45) +
  facet_grid(plate~date)
norm_platemap

#' ##Dirty group comparison
#' 
#' since I am going to exclude the anaerobic data in the later part of the this analysis here I will look at it
#' 
#+ da-plots, message=FALSE, fig.width=12, fig.height=10
summary1plot <- ggplot(norm_data, aes(sample_id, OD_adjusted)) +
  scale_colour_brewer(palette = "YlOrRd") +
  geom_boxplot(outlier.shape=NA) +
  geom_point(aes(shape = plate, colour = date), position = "jitter") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  scale_y_log10()
summary1plot #adjusted OD by detailed groupings

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
save(norm_data, file="Data/crystal_violet.RData")
