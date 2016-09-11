#' ---
#' title: "merged biofilm assays for publication, analyize and visualize the data"
#'output:
#'  html_document:
#'    toc: true
#'    theme: united
#' ---

#/*    Analysis script for statistics and visualization of crystal violet biofilm assay*/
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
require(lmerTest)
require(multcomp)
require(dunn.test)
require(orddom)
require(equivalence)
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
load("Data/crystal_violet_biofilm.RData")

date_index <- norm_data[simple_id %in% c("atl_", "icaA_", "srtA_"),unique(date)] #find all experiments with the agr sample
norm_data <- norm_data[date %in% date_index][simple_id %in% c("wt_NA", "atl_", "icaA_", "srtA_")] #extract them from the data with thier matched wt

#' 
#' #Initial Visualize
#' 

#+ dirty-visualze
#+ data-summary
#data_summary <- norm_data %>%
#  group_by(simple_id) %>%
#  summarise(
#    mean = mean(OD_adjusted, na.rm = TRUE), # means comparison
#    sdev = sd(OD_adjusted, na.rm = TRUE),
#    ci_lower = t.test(OD_adjusted)$conf.int[1], #95% confidence intervals CANT DO IT CAUSE THE DATA?
#    ci_upper = t.test(OD_adjusted)$conf.int[2])
#data_summary

#+ summary-plot, message=FALSE, fig.width=12, fig.height=10
#pos = position_dodge(width = 0.9)#for error bars to dodge dodging columns
#data_summary_plot <- ggplot(data_summary, aes(simple_id, mean, ymin = ci_lower, ymax = ci_upper)) +
#  geom_bar(aes(fill = simple_id), stat="identity", position = pos, width = 0.9) +
#  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
#  geom_errorbar(aes(fill = simple_id), width = 0.2, position = pos)
#data_summary_plot

#' 
#' #Set values
#'
#' These variables affect the coming scripts.
#' blocking_factor is the column (data frame variable) in which all experimental conditions are tested in
#' test_factor1 is the column (data frame variable) which contains the experimental conditions
#'

#+ set
transformed_data <- c("OD", "OD_adjusted") #different transformations/normalizations
alpha_level <- 0.05
blocking_factor <- "drop" #this is not a blocking factor, it includes all data
test_factor1 <- "simple_id"
set_data <- transformed_data[[2]]
set_test <- "ordinal" #ordinal <- non-paramtric, cliff's D, metric <- parametric, cohen's D

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
      combinations <- combn(unique(blocked_data[[test_factor1]]), 2)
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
l4 <- expression(atop("wild-type vs. ", paste("srtA", ":", ":", "erm")))
l5 <- expression(atop("wild-type vs. ", paste("icaA", ":", ":", "erm")))
l6 <- expression(atop("wild-type vs. ", italic(Delta*atl)))

#+ filtered-plot
orddom_sliced <- slice(orddom_summary, c(4:6)) #filter
effsize_plot <- ggplot(orddom_sliced, aes(delta, comparison)) +
  geom_vline(xintercept =  0, linetype = 2, alpha = 0.5) +
  geom_point() +
  geom_errorbarh(height = 0.2, aes(xmin = CI.low, 
                                   xmax = CI.high)) +
  labs(x = expression(paste("Cliff's ", Delta)), y = "Comparison") +
  scale_y_discrete(labels = c(l4, l5, l6)) +
  coord_cartesian(xlim = c(-1, 1)) +
#  facet_grid(.~timepoint) +
#  facet_grid(reformulate(".", facet_type)) +
  theme_mod
effsize_plot #plot

#+ save-graph, eval=FALSE
ggsave("Figures/crystal_violet_supp_effsize.tiff", plot = effsize_plot, width = 30, height = 10, units = "cm", dpi = 1200) #this code is only evaluate when the script us run ourside of knitr

#' #NHST
#' 
#' null hypothesis: the cfu from mutants are not different
#' since the fomer analysis was done with the non-parametric cliff's delta we are using the kruskall-wallis with dunn's post hoc here
#' 

#non parametric tests
#+ Kruskall-Wallis
kruskal.test(OD_adjusted ~ simple_id, data = norm_data) #Kruskall-Wallis test for a group with stochastic dominance

#+ Post-hoc-pairwise, message=FALSE
dunns_test <- dunn.test(norm_data$OD_adjusted, norm_data$simple_id, method = 'bonferroni') #pairwise comparison if Kruskall-Wallis is rejected using Dunn's Z statistic
sig_matrix1 <- matrix(c(dunns_test$Z, dunns_test$P, dunns_test$P.adjusted), 6, 3, dimnames = list(dunns_test$comparisons, c("Z", "pvalue", "padjusted")))
sig_matrix1 <- data.table(sig_matrix1, keep.rownames = TRUE)
sig_matrix1[,plabel := numeric_to_label(padjusted, prefix = "p=="),]

#'
#' #Main figure
#' 
#' Labels
#'

#+ labels
l1 <- "wild-type"
l3 <- expression(italic(Delta*atl))
l4 <- expression(italic(paste("icaA", ":", ":", "erm")))
l5 <- expression(italic(paste("srt", ":", ":", "erm")))
facet1 <- data.frame(x = 1:4, y = 1:4) #for the overlay later to allow for drawing stats comparison paths
norm_data$simple_id <- factor(norm_data$simple_id, levels(norm_data$simple_id)[c(17,6,11,14)])
comparisons <- list(c(1,2), c(1,3), c(1,4)) #state the comparisons

#+ overview-plot, fig.width=7, fig.height=7
main_fig <- ggplot(norm_data, aes(simple_id, OD_adjusted)) +
  geom_hline(yintercept =  1, linetype = 2, alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.5) +
  geom_point(position = position_jitter(width = 0.25), size = 0.3) +
  coord_cartesian(ylim = c(0,4.5)) +
  scale_x_discrete(labels = c(l1, l3, l4, l5)) +
  labs(x = " ", y = expression(normalized~OD[595])) +
  theme_mod + 
  geom_path(aes(x=rep(comparisons[[1]], each = 2),y=c(3,3.1,3.1,3)), data = facet1) +
  geom_text(aes(x=median(comparisons[[1]]),y=3.2,label=sig_matrix1[[1,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[2]], each = 2),y=c(3.6,3.7,3.7,3.6)), data = facet1) +
  geom_text(aes(x=median(comparisons[[2]]),y=3.8,label=sig_matrix1[[2,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[3]], each = 2),y=c(4.2,4.3,4.3,4.2)), data = facet1) +
  geom_text(aes(x=median(comparisons[[3]]),y=4.4,label=sig_matrix1[[4,5]]), data = facet1, parse = TRUE)
main_fig

#+ save-graph2, eval=FALSE
ggsave("Figures/crystal_violet_supp_biofilm.tiff", plot = main_fig, width = 25, height = 25, units = "cm", dpi = 1200) #this code is only evaluate when the script us run ourside of knitr

