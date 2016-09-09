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

#' #NHST
#' 
#' null hypothesis: the cfu from mutants are not different
#' we are using the kruskall-wallis with dunn's post hoc on the untrasnformed data to understand which group are different from each other for a given timepoint
#' 
#+split
splitted <- split(norm_data, norm_data[[blocking_factor]], drop = TRUE)

#non parametric tests
#+ Kruskall-Wallis
map(splitted, function(x) kruskal.test(cfu ~ sample_id, data = x)) #Kruskall-Wallis test for a group with stochastic dominance

#+ Post-hoc-pairwise, message=FALSE
dunns_test <- dunn.test(splitted[[1]]$cfu, splitted[[1]]$sample_id, method = 'bonferroni') #pairwise comparison if Kruskall-Wallis is rejected using Dunn's Z statistic
sig_matrix1 <- matrix(c(dunns_test$Z, dunns_test$P, dunns_test$P.adjusted), 6, 3, dimnames = list(dunns_test$comparisons, c("Z", "pvalue", "padjusted")))
sig_matrix1 <- data.table(sig_matrix1, keep.rownames = TRUE)
sig_matrix1[,plabel := numeric_to_label(padjusted, prefix = "p=="),]

dunns_test <- dunn.test(splitted[[2]]$cfu, splitted[[2]]$sample_id, method = 'bonferroni') #pairwise comparison if Kruskall-Wallis is rejected using Dunn's Z statistic
sig_matrix2 <- matrix(c(dunns_test$Z, dunns_test$P, dunns_test$P.adjusted), 6, 3, dimnames = list(dunns_test$comparisons, c("Z", "pvalue", "padjusted")))
sig_matrix2 <- data.table(sig_matrix2, keep.rownames = TRUE)
sig_matrix2[,plabel := numeric_to_label(padjusted, prefix = "p=="),]


#'
#' #Main figure
#' 
#' Labels
#'

#+ labels
l1 <- "wild-type"
l2 <- expression(italic(agrA[C123F]))
l3 <- expression(atop(italic(agrA[C123F]), "+ pOS1 empty"))
l4 <- expression(atop(italic(agrA[C123F]), paste("+ pOS1 ", italic(agrA))))


facet1 <- data.frame(x = 1:4, y = 1:4, timepoint = "72 hours") #for the overlay later to allow for drawing stats comparison paths
facet2 <- data.frame(x = 1:4, y = 1:4, timepoint = "120 hours") #for the overlay later to allow for drawing stats comparison paths
comparisons <- list(c(1,2), c(1,3), c(1,4), c(3,4)) #state the comparisons
p <- logscale_sigbars_generator(5e10, 2e9, 3, text_spacing = 1.5) # to calculate the postions for statisitical significance bars on a log scale

#+ overview-plot, fig.width=7, fig.height=7
main_fig <- ggplot(norm_data, aes(sample_id, cfu)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.5) +
  geom_point(position = position_jitter(width = 0.25), size = 0.3) +
  scale_x_discrete(labels = c(l1, l2, l3, l4)) +
  labs(x = " ", y = "CFU / raft") +
  coord_cartesian(ylim = c(1e5, 1e11)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  facet_wrap(~timepoint) +
  theme_mod + 
  #facet 1, labels are still manually added as is the level which subsets p in each path-text pair
  geom_path(aes(x=rep(comparisons[[1]], each = 2),y=p[1,1:4]), data = facet1) +
  geom_text(aes(x=median(comparisons[[1]]),y=p[1,5],label=sig_matrix1[[1,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[2]], each = 2),y=p[2,1:4]), data = facet1) +
  geom_text(aes(x=median(comparisons[[2]]),y=p[2,5],label=sig_matrix1[[2,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[3]], each = 2),y=p[3,1:4]), data = facet1) +
  geom_text(aes(x=median(comparisons[[3]]),y=p[3,5],label=sig_matrix1[[4,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[4]], each = 2),y=p[1,1:4]), data = facet1) +
  geom_text(aes(x=median(comparisons[[4]]),y=p[1,5],label=sig_matrix1[[6,5]]), data = facet1, parse = TRUE) +
  # facet 2
  geom_path(aes(x=rep(comparisons[[1]], each = 2),y=p[1,1:4]), data = facet2) +
  geom_text(aes(x=median(comparisons[[1]]),y=p[1,5],label=sig_matrix2[[1,5]]), data = facet2, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[2]], each = 2),y=p[2,1:4]), data = facet2) +
  geom_text(aes(x=median(comparisons[[2]]),y=p[2,5],label=sig_matrix2[[2,5]]), data = facet2, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[3]], each = 2),y=p[3,1:4]), data = facet2) +
  geom_text(aes(x=median(comparisons[[3]]),y=p[3,5],label=sig_matrix2[[4,5]]), data = facet2, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[4]], each = 2),y=p[1,1:4]), data = facet2) +
  geom_text(aes(x=median(comparisons[[4]]),y=p[1,5],label=sig_matrix2[[6,5]]), data = facet2, parse = TRUE)
main_fig

#+ save-graph2, eval=FALSE
ggsave("Figures/rafts_agr_cfu.tiff", plot = main_fig, width = 30, height = 25, units = "cm", dpi = 1200) #this code is only evaluate when the script us run ourside of knitr
