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

#+ functions, include=FALSE
logscale_sigbars_generator <- function (max_draw_dim, min_draw_dim, number_bar_levels = 1, tick_size = 0.01, default_step = 1.5, text_spacing = 2) {
  #someday it'll be nice to have some input verification
  print("positions generated:", quote = FALSE)
  print("the levels are counted from the bottom to top",quote = FALSE) 
  print("p[level, 1:4] are the positions", quote = FALSE)
  print("p[level, 5] contains the text position", quote = FALSE) #some guidance
  range <- log(max_draw_dim) - log(min_draw_dim) #the range that the bars will be plotted in
  tick_size_log <- log(max_draw_dim) * tick_size # the size of the downturned ticks
  step <- ifelse((range / number_bar_levels) < default_step, (range / number_bar_levels), default_step) # the spacing between bars
  p <- matrix(0,number_bar_levels, 5) # the matrix of the results
  for (i in 1:number_bar_levels) {
    bar_position <- log(min_draw_dim) + (step*i)
    tick_postion <- bar_position - tick_size_log
    text_position <- bar_position + (text_spacing * tick_size_log)
    p[i,] <- as.numeric(c(exp(tick_postion), exp(bar_position), exp(bar_position),exp(tick_postion), exp(text_position)))
  } #make it
  return(p)
}

numeric_to_label <- function (numeric, prefix = "", precision = 3, round_method = floor) {
  value <- as.integer(numeric)
  exponent <- round_method(log10(numeric))
  number <- signif(numeric * 10^(exponent * -1), digits = precision)
  label <- paste0(prefix, number, "%*%10^{", exponent, "}")
  return(label)
}

#+ style, inlcude=FALSE
theme_mod <- theme_bw() +
  theme(text = element_text(size = 16),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())
#and some xtable formatting options too
italic <- function(x){
  paste0('{\\emph{', x, '}}')
}
bold <- function(x){
  paste0('{\\bfseries ', x, '}')
}

#+ Session-info
sessionInfo() #for reproducibility

#+ some-other-options
original_par <- par() #for resetting to original par after generating the plot in the null device

#+ load-data
load("Data/crystal_violet_biofilm.RData")
norm_data <- norm_data[simple_id %in% c("wt_NA", "agrA_", "atl_", "icaA_", "srtA_", "agrA_pos1_agrA_20", "agrA_pos1_empty")] #filter

#' 
#' #Initial Visualize
#' 

#+ dirty-visualze
#+ data-summary
data_summary <- norm_data %>%
  group_by(simple_id) %>%
  summarise(
    mean = mean(OD_adjusted, na.rm = TRUE), # means comparison
    sdev = sd(OD_adjusted, na.rm = TRUE),
    ci_lower = t.test(OD_adjusted)$conf.int[1], #95% confidence intervals CANT DO IT CAUSE THE DATA?
    ci_upper = t.test(OD_adjusted)$conf.int[2])
data_summary

#+ summary-plot, message=FALSE, fig.width=12, fig.height=10
pos = position_dodge(width = 0.9)#for error bars to dodge dodging columns
data_summary_plot <- ggplot(data_summary, aes(simple_id, mean, ymin = ci_lower, ymax = ci_upper)) +
  geom_bar(aes(fill = simple_id), stat="identity", position = pos, width = 0.9) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  geom_errorbar(aes(fill = simple_id), width = 0.2, position = pos)
data_summary_plot

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
l1 <- expression(atop(paste(italic(agrA[C123F]), " + pOS1 ", italic(agrA), " vs."), paste(italic(agrA[C123F]), " + pOS1 empty")))
l2 <- expression(atop("wild-type vs. ", paste(italic(agrA[C123F]), " + pOS1 empty")))
l3 <- expression(atop("wild-type vs. ", paste(italic(agrA[C123F]), " + pOS1 ", italic(agrA))))
l4 <- expression(atop("wild-type vs. ", paste("srtA", ":", ":", "erm")))
l5 <- expression(atop("wild-type vs. ", paste("icaA", ":", ":", "erm")))
l6 <- expression(atop("wild-type vs. ", italic(Delta*atl)))
l7 <- expression(atop("wild-type vs. ", italic(agrA[C123F])))

#+ filtered-plot
orddom_sliced <- slice(orddom_summary, c(1,16:21)) #filter
effsize_plot <- ggplot(orddom_sliced, aes(delta, comparison)) +
  geom_vline(xintercept =  0, linetype = 2, alpha = 0.5) +
  geom_point() +
  geom_errorbarh(height = 0.2, aes(xmin = CI.low, 
                                   xmax = CI.high)) +
  labs(x = expression(paste("Cliff's ", Delta)), y = "Comparison") +
  scale_y_discrete(labels = c(l7, l6, l5, l4, l3, l2, l1)) +
  coord_cartesian(xlim = c(-1, 1)) +
#  facet_grid(.~timepoint) +
#  facet_grid(reformulate(".", facet_type)) +
  theme_mod
effsize_plot #plot

#+ save-graph, eval=FALSE
ggsave("Figures/crystal_violet_effsize.tiff", plot = effsize_plot, width = 30, height = 15, units = "cm", dpi = 1200) #this code is only evaluate when the script us run ourside of knitr

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
sig_matrix1 <- matrix(c(dunns_test$Z, dunns_test$P, dunns_test$P.adjusted), 21, 3, dimnames = list(dunns_test$comparisons, c("Z", "pvalue", "padjusted")))
sig_matrix1 <- data.table(sig_matrix1, keep.rownames = TRUE)
sig_matrix1[,plabel := numeric_to_label(padjusted, prefix = "p=="),]

#'
#' #Main figure
#' 
#' Labels
#'

#+ labels
l1 <- "wild-type"
l2 <- expression(italic(agrA[C123F]))
l3 <- expression(italic(Delta*atl))
l4 <- expression(italic(paste("icaA", ":", ":", "erm")))
l5 <- expression(italic(paste("srt", ":", ":", "erm")))
l6 <- expression(atop(italic(agrA[C123F]), paste("+ pOS1 ", italic(agrA))))
l7 <- expression(atop(italic(agrA[C123F]), "+ pOS1 empty"))
facet1 <- data.frame(x = 1:4, y = 1:4) #for the overlay later to allow for drawing stats comparison paths
norm_data$simple_id <- factor(norm_data$simple_id, c("wt_NA", "srtA_", "icaA_", "atl_", "agrA_", "agrA_pos1_agrA_20", "agrA_pos1_empty")) #reorder factor
comparisons <- list(c(1,2), c(1,3), c(1,4), c(1,4.9), c(1,6), c(1,7), c(6,7)) #state the comparisons

#+ overview-plot, fig.width=7, fig.height=7
main_fig <- ggplot(norm_data, aes(simple_id, OD_adjusted)) +
  geom_hline(yintercept =  1, linetype = 2, alpha = 0.5) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.5) +
  geom_point(position = position_jitter(width = 0.25), size = 0.3) +
  coord_cartesian(ylim = c(0,12)) +
  scale_x_discrete(labels = c(l1, l5, l4, l3, l2, l6, l7)) +
  labs(x = " ", y = expression(normalized~OD[595])) +
  theme_mod + 
  geom_path(aes(x=rep(comparisons[[1]], each = 2),y=c(6,6.2,6.2,6)), data = facet1) +
  geom_text(aes(x=median(comparisons[[1]]),y=6.4,label=sig_matrix1[[1,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[2]], each = 2),y=c(7,7.2,7.2,7)), data = facet1) +
  geom_text(aes(x=median(comparisons[[2]]),y=7.4,label=sig_matrix1[[2,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[3]], each = 2),y=c(8,8.2,8.2,8)), data = facet1) +
  geom_text(aes(x=median(comparisons[[3]]),y=8.4,label=sig_matrix1[[4,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[4]], each = 2),y=c(9,9.2,9.2,9.2)), data = facet1) + #modified to not cross the whisker
  geom_text(aes(x=median(comparisons[[4]]),y=9.4,label=sig_matrix1[[7,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[5]], each = 2),y=c(10,10.2,10.2,10)), data = facet1) +
  geom_text(aes(x=median(comparisons[[5]]),y=10.4,label=sig_matrix1[[11,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[6]], each = 2),y=c(11,11.2,11.2,11)), data = facet1) +
  geom_text(aes(x=median(comparisons[[6]]),y=11.4,label=sig_matrix1[[16,5]]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[7]], each = 2),y=c(6,6.2,6.2,6)), data = facet1) +
  geom_text(aes(x=median(comparisons[[7]]),y=6.4,label=sig_matrix1[[21,5]]), data = facet1, parse = TRUE)
main_fig

#+ save-graph2, eval=FALSE
ggsave("Figures/crystal_violet_biofilm.tiff", plot = main_fig, width = 30, height = 15, units = "cm", dpi = 1200) #this code is only evaluate when the script us run ourside of knitr

