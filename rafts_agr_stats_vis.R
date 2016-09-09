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

norm_data$sample_id <- factor(norm_data$sample_id, levels(norm_data$sample_id)[c(5,1,9,7)])
norm_data$timepoint <- plyr::mapvalues(norm_data$timepoint, c("72", "120"), c("72 hours", "120 hours"))

#'
#' #Effect Size
#'
#' Effect size is calculated and plotted
#'

#+ calculate
orddom_cols <- c("comparison", "timepoint", "Var2", "1-alpha", "A X>Y", "A Y>X", "CI high", 
                 "CI low", "Cohen's d", "d CI high", "d CI low", "delta", "df", "H1 tails p/CI", 
                 "N #Y<X", "N #Y=X", "N #Y>X", "n in X", "n in Y", "NNT", "p", "PS X>Y", "PS Y>X", 
                 "s delta", "se delta", "type_title", "var d.i", "var delta", "var dij", "var dj.", 
                 "var1_X", "var2_Y", "z/t score") #column names for orddom ouput
orddom_summary <- norm_data %>%
  split(norm_data[[blocking_factor]], drop = TRUE) %>%
  map(
    function(blocked_data) {
      combinations <- combn(sort(unique(blocked_data[[test_factor1]])), 2)
      df <- data.frame()
      for (i in 1:dim(combinations)[[2]]) {
        id1 <- unlist(combinations[,i])[[1]]
        id2 <- unlist(combinations[,i])[[2]]
        data1 <- blocked_data[blocked_data[[test_factor1]] == id1,]
        data2 <- blocked_data[blocked_data[[test_factor1]] == id2,]
        comparison <- paste(id1, id2, sep = "-")
        #transformed data needs to know a position
        data <- cbind(melt(orddom(data1[[set_data]], data2[[set_data]], 
                                  alpha = 0.05, 
                                  symmetric = FALSE, 
                                  onetailed = FALSE,
                                  t.welch = TRUE)),
                      comparison)
        df <- rbind(data, df)}
      return(df)
    }) %>% #calculate orddom matrix by blocking factor
  at_depth(0, ~ bind_rows(.x, .id = "timepoint")) %>% #bind the data frames together by row
  dcast(comparison + timepoint + Var2 ~ Var1) %>% #merge intp final data frame
  map_at(orddom_cols[c(-1, -2, -3, -26,-31, -32)], as.numeric) %>%
  as.data.frame() %>%
  filter(Var2 == set_test)

#+ comparisons-label
l4 <- expression(atop(paste(italic(agrA[C123F]), " + pOS1 ", italic(agrA), " vs."), paste(italic(agrA[C123F]), " + pOS1 empty")))
l3 <- expression(atop("wild-type vs. ", paste(italic(agrA[C123F]), " + pOS1 ", italic(agrA))))
l2 <- expression(atop("wild-type vs. ", paste(italic(agrA[C123F]), " + pOS1 empty")))
l1 <- expression(atop("wild-type vs. ", italic(agrA[C123F])))

#+ filtered-plot
orddom_sliced <- slice(orddom_summary, c(1,2,7:12)) #filter
orddom_sliced$timepoint <- factor(orddom_sliced$timepoint, levels(orddom_sliced$timepoint)[c(2,1)]) #reorder factor

effsize_plot <- ggplot(orddom_sliced, aes(delta, comparison)) +
  geom_vline(xintercept =  0, linetype = 2, alpha = 0.5) +
  geom_point() +
  geom_errorbarh(height = 0.2, aes(xmin = CI.low, 
                                   xmax = CI.high)) +
  labs(x = expression(paste("Cohen's ", italic(d))), y = "Comparison") +
  scale_y_discrete(labels = c(l4, l3, l2, l1)) +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  #  facet_grid(.~timepoint) +
  facet_grid(reformulate(blocking_factor, ".")) +
  theme_mod
effsize_plot #plot

#+ save-graph, eval=FALSE
ggsave("Figures/rafts_agr_effsize.tiff", plot = effsize_plot, width = 30, height = 10, units = "cm", dpi = 1200) #this code is only evaluate when the script us run ourside of knitr
