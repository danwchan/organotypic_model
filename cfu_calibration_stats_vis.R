#' ---
#' title: "Statistics and visualization of CFU caibration assay"
#'output:
#'  html_document:
#'    toc: true
#'    theme: united
#' ---

#/*    Analysis script for statistics and visualization of cfu calibration assay with wild type bacteria*/
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
load("Data/cfu_calibration.RData")

#' 
#' #Initial Visualize
#' 

#+ dirty-visualze
#+ data-summary
data_summary <- norm_data %>%
  group_by(id_merge) %>%
  summarise(
    mean = mean(cfu_log, na.rm = TRUE), # means comparison
    sdev = sd(cfu_log, na.rm = TRUE),
    ci_lower = t.test(cfu_log)$conf.int[1],
    ci_upper = t.test(cfu_log)$conf.int[2])
data_summary

#+ summary-plot, message=FALSE, fig.width=12, fig.height=10
pos = position_dodge(width = 0.9)#for error bars to dodge dodging columns
data_summary_plot <- ggplot(data_summary, aes(id_merge, mean, ymin = ci_lower, ymax = ci_upper)) +
  geom_bar(aes(fill = id_merge), stat="identity", position = pos, width = 0.9) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  geom_errorbar(aes(fill = id_merge), width = 0.2, position = pos)
data_summary_plot

#' #Effect Size
#' 
#+ expermental-variables
#getting more generalized
blocking_factor <- "timepoint"  #the blocking factor is timepoint (operations are performed within blocks), required for ANOVA + orddom
test_factor1 <- "CFU_delivered"
transformed_data <- c("cfu", "cfu_log", "cfu_scale")
alpha_level <- 0.05

#norm_data$timepoint <- factor(norm_data$timepoint, levels(norm_data$timepoint)[1,3,4,2]) #reorder blocks
norm_data$timepoint <- plyr::mapvalues(norm_data$timepoint, c("0", "24", "72", "120"), c("3.5 hours", "24 hours", "72 hours", "120 hours")) #rename factors

#' 
#' ##Ordinal transformed effect size (non-parametric)
#' 
#+ set-variables
set_data <- transformed_data[[2]]
set_test <- "metric"
#' the following code chunks are run using the transformed data 
{{as.name(set_data)}} 
#' to obtain 
{{as.name(set_test)}} 
#' effect
#' 
#' sizes calculated by the package orddom 
{{packageVersion("orddom")}}
#' 

#+ pairwise-effsize-via-orddom
#rewrote this function so I don't have to create the pairwise data table, function (split) iterates as required over the split
#data and compares the transformed cfu's indicated by the argument choice (see the experimental-varibles chunk)
#
orddom_cols <- c("comparison", "timepoint", "Var2", "1-alpha", "A X>Y", "A Y>X", "CI high", 
                 "CI low", "Cohen's d", "d CI high", "d CI low", "delta", "df", "H1 tails p/CI", 
                 "N #Y<X", "N #Y=X", "N #Y>X", "n in X", "n in Y", "NNT", "p", "PS X>Y", "PS Y>X", 
                 "s delta", "se delta", "type_title", "var d.i", "var delta", "var dij", "var dj.", 
                 "var1_X", "var2_Y", "z/t score") #column names for orddom ouput
orddom_summary <- norm_data %>%
  split(norm_data[[blocking_factor]]) %>%
  map(
    function(split) {
      combinations <- combn(unique(split[[test_factor1]]), 2)
      df <- data.frame()
      for (i in 1:dim(combinations)[[2]]) {
        id1 <- unlist(combinations[,i])[[1]]
        id2 <- unlist(combinations[,i])[[2]]
        data1 <- split[split[[test_factor1]] == id1,]
        data2 <- split[split[[test_factor1]] == id2,]
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

#+ labels-effsize
label1 <- expression(atop("13 vs 156", "CFU / raft"))
label2 <- expression(atop("13 vs 1680", "CFU / raft"))
label3 <- expression(atop("156 vs 1680", "CFU / raft"))
orddom_summary$timepoint <- factor(orddom_summary$timepoint, levels(orddom_summary$timepoint)[c(3,2,4,1)]) #reorder factor

#+ effsize-plot
ordinal_effsize <- ggplot(orddom_summary, aes(delta, comparison)) + #plot the data
  geom_vline(xintercept =  0, linetype = 2, alpha = 0.5) +
  geom_point() +
  geom_errorbarh(height = 0.2, aes(xmin = CI.low, 
                                   xmax = CI.high)) +
  labs(x = expression(paste("Cohen's ", italic(d))), y = "Comparison") +
  scale_y_discrete(labels = c(label3, label2, label1)) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  facet_grid(~timepoint)
ordinal_effsize

#+ save-graph, eval=FALSE
ggsave("Figures/cfu_calibration_effsize.tiff", plot = ordinal_effsize, width = 30, height = 15, units = "cm", dpi = 1200) #this code is only evaluate when the script us run ourside of knitr

#' #NHST
#' 
#' null hypothesis: the cfu from mutants are not different
#' 
#' ##parametric tests
#' 
#' (using log transformed values, comparing means as a central tendency) robust(?) against non-normality and heterodiascity
#' 

#+ ANOVA-old
anova_fit <- aov(cfu_log ~ CFU_delivered*timepoint, data = norm_data)

#+ ANOVA-diagnostic-plots
par(mfcol = c(2,2), mar = c(2,2,2,1))
plot(anova_fit)
par(original_par)
summary(anova_fit) #Type I SS, F-test
drop1(anova_fit,~.,test="F") # type III SS and F Tests

#+ Tukey-HSD
#for when the anova did not include the error factor of date
tukey_anova_fit <- TukeyHSD(anova_fit, "CFU_delivered") #tukey's HSD pairwise if ANOVA null is rejected
par(mar = c(5,9,3,1))
plot(tukey_anova_fit, las = 1) #visulaize
par(original_par)

#+ t-test
norm_data %>%
  split(norm_data[[blocking_factor]]) %>%
  map(function(x) pairwise.t.test(x$cfu_log, x$CFU_delivered, p.adjust.method = 'bonferroni', pool.sd = FALSE)) #t-tests if the ANOVA null is rejected

#' ##non-parametric tests
#' (using ranks, comparing medians as central tendency)

#+ Kruskall-Wallis, message=FALSE
kruskal_sig <- norm_data %>%
  split(norm_data[[blocking_factor]]) %>%
  map(function(x) kruskal.test(cfu ~ CFU_delivered, data = x)) %>% #Kruskall-Wallis test for a group with stochastic dominance
  map(function(x) matrix(c(x$statistic, x$parameter, x$p.value), 1, 3, dimnames = list(x[[blocking_factor]], c("chi-squared", "df", "p-value")))) %>% 
  print() %>% #summarize in matrix and print to stdout
  lapply(function(x) x[[3]] <= alpha_level) %>% #test for significance
  as.logical()

#+ Post-hoc-pairwise, message=FALSE
norm_data %>%
  split(norm_data[[blocking_factor]]) %>%
  discard(!kruskal_sig) %>% #discard those that did not pass the Kruskall wallis significance test
  map(function(x) pairwise.wilcox.test(x$cfu, x$CFU_delivered, p.adjust.method = 'bonferroni')) %>% #pairwise comparison if Kruskall-Wallis is rejected using Wilcoxon distrubution
  map(function(x) x["p.value"]) #print wilcox tests

sig_matrix1 <- norm_data %>%
  split(norm_data[[blocking_factor]]) %>%
  discard(!kruskal_sig) %>% #discard those that did not pass the Kruskall wallis significance test
  map(function(x) dunn.test(x$cfu, x$CFU_delivered, method = 'bonferroni', table = FALSE, kw = FALSE)) %>% #pairwise comparison if Kruskall-Wallis is rejected using Dunn's Z statistic
  map(function(x) matrix(c(x$Z, x$P, x$P.adjusted), 3, 3, dimnames = list(x$comparisons, c("Z", "pvalue", "padjusted")))) %>%
  map(function (x) as.data.frame(x)) %>%
  map(function (x) cbind(x, plabel = numeric_to_label(x$padjusted, prefix = "p==")))

#+ labels-setup
label1 <- "13 CFU"
label2 <- "156 CFU"
label3 <- "1680 CFU"
#for the overlay layer to allow for drawing stats comparison paths, one df for each facet
facet1<- data.frame(x = 1:4, y = 1:4, timepoint = "3.5 hours") 
facet2 <- data.frame(x = 1:4, y = 1:4, timepoint = "24 hours") 
facet3 <- data.frame(x = 1:4, y = 1:4, timepoint = "72 hours")
facet4 <- data.frame(x = 1:4, y = 1:4, timepoint = "120 hours")
#state the comparisons
comparisons <- list(c(1,1.9), c(2.1,3), c(1,3))
# to calculate the postions for statisitical significance bars o a log scale
p <- logscale_sigbars_generator(1e11, 9e8, 2)

#+ overview-plot, fig.width=7, fig.height=7
cfu_calibration <- ggplot(norm_data, aes(CFU_delivered, cfu)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.5) +
  geom_point(position = position_jitter(width = 0.25), size = 0.3) +
  scale_x_discrete(labels = c(label1, label2, label3)) +
  labs(x = "mean CFU delivered / raft", y = "CFU / raft") +
  coord_cartesian(ylim = c(1, 1e11)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  facet_grid(~timepoint) +
  theme_mod +
  #facet 1
  geom_path(aes(x=rep(comparisons[[1]], each = 2),y=p[1,1:4]), data = facet1) +
  geom_text(aes(x=median(comparisons[[1]]),y=p[1,5],label=sig_matrix1[['3.5 hours']][1,4]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[2]], each = 2),y=p[1,1:4]), data = facet1) +
  geom_text(aes(x=median(comparisons[[2]]),y=p[1,5],label=sig_matrix1[['3.5 hours']][3,4]), data = facet1, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[3]], each = 2),y=p[2,1:4]), data = facet1) +
  geom_text(aes(x=median(comparisons[[3]]),y=p[2,5],label=sig_matrix1[['3.5 hours']][2,4]), data = facet1, parse = TRUE) +
  # facet 2
  geom_path(aes(x=rep(comparisons[[1]], each = 2),y=p[1,1:4]), data = facet2) +
  geom_text(aes(x=median(comparisons[[1]]),y=p[1,5],label=sig_matrix1[['24 hours']][1,4]), data = facet2, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[2]], each = 2),y=p[1,1:4]), data = facet2) +
  geom_text(aes(x=median(comparisons[[2]]),y=p[1,5],label=sig_matrix1[['24 hours']][3,4]), data = facet2, parse = TRUE) +
  geom_path(aes(x=rep(comparisons[[3]], each = 2),y=p[2,1:4]), data = facet2) +
  geom_text(aes(x=median(comparisons[[3]]),y=p[2,5],label=sig_matrix1[['24 hours']][2,4]), data = facet2, parse = TRUE)
cfu_calibration

#+ save-graph2, eval=FALSE
ggsave("Figures/cfu_calibration.tiff", plot = cfu_calibration, width = 30, height = 15, units = "cm", dpi = 1200) #this code is only evaluate when the script us run ourside of knitr

