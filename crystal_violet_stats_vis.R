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
logscale_sigbars_generator <- function (max_draw_dim, min_draw_dim, number_bar_levels = 1, tick_size = 0.01, default_step =1.5) {
  #someday it'll be nice to have some input verification
  print("positions generated:", quote = FALSE)
  print("the levels are counted from the bottom to top",quote = FALSE) 
  print("p[level, 1:4] are the positions", quote = FALSE)
  print("p[level, 5] contains the text position", quote = FALSE) #some guidance
  range <- log(max_draw_dim) - log(min_draw_dim) #the range that the bars will be plotted in
  tick_size_log <- log(max_draw_dim) * tick_size # the size of the downturned ticks
  step <- ifelse((range / number_bar_levels) < 1.5, (range / number_bar_levels), 1.5) # the spacing between bars
  p <- matrix(0,number_bar_levels, 5) # the matrix of the results
  for (i in 1:number_bar_levels) {
    bar_position <- log(min_draw_dim) + (step*i)
    tick_postion <- bar_position - tick_size_log
    text_position <- bar_position + (2 * tick_size_log)
    p[i,] <- as.numeric(c(exp(tick_postion), exp(bar_position), exp(bar_position),exp(tick_postion), exp(text_position)))
  } #make it
  return(p)
}
effect_size_plot <- function(orddom_table) {
  cliffd_plot <- ggplot(orddom_table, aes(delta, comparison)) +
    geom_vline(xintercept =  0, linetype = 2, alpha = 0.5) +
    geom_point() +
    geom_errorbarh(height = 0.2, aes(xmin = CI.low, 
                                     xmax = CI.high)) +
    labs(x = "Cliff's Delta", y = "Comparison") +
    coord_cartesian(xlim = c(-1, 1)) +
    facet_grid(.~timepoint) +
    theme_mod
  return(cliffd_plot)
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
norm_data <- norm_data[simple_id %in% c("wt_NA", "agrA_", "atl_", "icaA_", "srtA_", "agrA_pos1_agrA_20", "agrA_pos1_empty")]

#' 
#' #Initial Visualize
#' 

#+ dirty-visualze
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

#' 
#' ##Set values
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
#' ##Effect Size
#'
#' Effect size is calculated and plotted
#'

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

#+ orddom-effect-size
orddom_cols <- c("comparison", "timepoint", "Var2", "1-alpha", "A X>Y", "A Y>X", "CI high", 
                 "CI low", "Cohen's d", "d CI high", "d CI low", "delta", "df", "H1 tails p/CI", 
                 "N #Y<X", "N #Y=X", "N #Y>X", "n in X", "n in Y", "NNT", "p", "PS X>Y", "PS Y>X", 
                 "s delta", "se delta", "type_title", "var d.i", "var delta", "var dij", "var dj.", 
                 "var1_X", "var2_Y", "z/t score") #column names for orddom ouput
orddom <- pairwise_data_table %>%
  at_depth(1, ~ orddom(.x$cfu1, .x$cfu2)) %>%
  at_depth(1, melt) %>% #manipulate to data frame
  at_depth(0, ~ bind_rows(.x, .id = "comparison")) %>% #bind the data frames together by row
  dcast(comparison + Var2 ~ Var1) %>% #merge intp final data frame
  map_at(orddom_cols[c(-1, -2, -3, -26,-31, -32)], as.numeric) #convert to numeric data
#make it back into a dataframe and filter to keep only ordinal information
orddom_ordinal <- as.data.frame(orddom) %>%
  filter(Var2 == "ordinal")
#make it back into a dataframe and filter to keep only metric information
orddom_metric <- as.data.frame(orddom) %>%
  filter(Var2 == "metric")

#+ save-KO_subsets, results="hide"
save(KO_subset, file = "KO_subset.RData") #the subsetted norm data
save(data_summary, file = "KO_subset_summary.RData") #the mean, sdev and t.test computed CI
save(orddom_ordinal, file = "KO_subset_cliffsd.RData") #cliff's d: effectsizes
save(orddom_metric, file = "KO_subset_metric.RData") #metric effectsizes

#+ plot-ordinal-effsize
ordinal_effsize <- ggplot(orddom_ordinal, aes(delta, comparison)) +
  geom_point() +
  geom_errorbarh(height = 0.2, aes(xmin = CI.low, 
                                   xmax = CI.high))
ordinal_effsize

#+ plot-metric-effsize
metric_effsize <- ggplot(orddom_metric, aes(delta, comparison)) +
  geom_point() +
  geom_errorbarh(height = 0.2, aes(xmin = CI.low, 
                                   xmax = CI.high))
metric_effsize

#' ##NHST
#' 
#' null hypothesis: the cfu from mutants are not different
#' 
#+ ANOVA-new
anova_fit <- lmer(OD_adjusted ~ simple_id*coating + (1|date), data = KO_subset, REML = TRUE) #making the mixed effects model 2 factors one random error

#+ ANOVA-new-stats
#some difficulty getting all of the plots into one window...
plot(anova_fit, type = c("p","smooth")) #fitted v residual
plot(anova_fit, sqrt(abs(resid(.)))~fitted(.), type = c("p","smooth")) #scale location
qqnorm(resid(anova_fit))
qqline(resid(anova_fit)) #qq
anova(anova_fit) #test if the factors make a difference
#'both timepoint and sample_id have a significant effect on cfu_log

#+ model-based-effect-sizes
effect_size <- difflsmeans(anova_fit, test.effs = "simple_id")
names(effect_size$diffs.lsmeans.table)[c(5,6)] <-c("Lower.CI", "Upper.CI") #remove the spaces from the names
effect_size_plot <- ggplot(effect_size$diffs.lsmeans.table, aes(Estimate, dimnames(effect_size$diffs.lsmeans.table)[[1]])) +
  geom_point() +
  geom_errorbarh(height = 0.2, aes(xmin = Lower.CI, xmax = Upper.CI))
effect_size_plot

#+ Tukey-test
posthoc <- glht(anova_fit, linfct = mcp(simple_id = "Tukey"))
tukey_anova_fit <- summary(posthoc, test = adjusted("single-step"))
tukey_anova_fit

#non parametric tests
#+ Kruskall-Wallis
kruskal.test(OD_adjusted ~ simple_id, data = KO_subset) #Kruskall-Wallis test for a group with stochastic dominance

#+ Post-hoc-pairwise, message=FALSE
wilcox_tests <- pairwise.wilcox.test(KO_subset$OD_adjusted, KO_subset$simple_id, p.adjust.method = 'bonferroni') #pairwise comparison if Kruskall-Wallis is rejected using Wilcoxon distrubution
wilcox_tests #print wilcox tests
dunns_test <- dunn.test(KO_subset$OD_adjusted, KO_subset$simple_id, method = 'bonferroni') #pairwise comparison if Kruskall-Wallis is rejected using Dunn's Z statistic
matrix(c(dunns_test$Z, dunns_test$P, dunns_test$P.adjusted), 10, 3, dimnames = list(dunns_test$comparisons, c("Z", "p-value", "p-adjusted")))

#' #Complemented mutants: agrA
#' 
#+ subset-data-agrA, results="hide"
agr_subset_index <- unique(norm_data[simple_id == "agrA_pos1_empty" | simple_id == "agrA_pos1_agrA_20", .(date, plate)])
agr_subset <- norm_data[agr_subset_index][simple_id == "agrA_pos1_empty" | simple_id == "agrA_pos1_agrA_20" | simple_id == "wt_NA"]
agr_subset$simple_id <- droplevels(agr_subset$simple_id)
OD_aggregate <- aggregate(cbind(OD, OD_adjusted, OD_log, OD_scale) ~ simple_id,
                          data = agr_subset,
                          print) #the OD for each pariwise combinations, also required for Cliff's DELTA
xselect <- combn(OD_aggregate[["simple_id"]], 2)
names(OD_aggregate$OD_adjusted) <- OD_aggregate$simple_id

#+ effect-size-agr
#rewrite effect size in data.table
data_summary <- agr_subset[simple_id != "wt_NA", summarise(.SD,
                                                          mean = mean(OD_adjusted),
                                                          sdev = sd(OD_adjusted),
                                                          ci_lower = t.test(OD_adjusted)$conf.int[1],
                                                          ci_upper = t.test(OD_adjusted)$conf.int[2]), by = .(simple_id)]
data_summary <- rbindlist(list(data_summary, 
                               agr_subset[simple_id == "wt_NA", summarise(.SD, mean = mean(OD_adjusted),
                                                                         sdev = sd(OD_adjusted)), by = .(simple_id)]),
                          fill = TRUE)

#+ effect-size-plot3, ref.label="effect-size-plot"

#+foda-se-amazing, message=FALSE
#HARDCODED someday you might want to make this nicer, I think it's close to being general
pairwise_data_table <- data.frame(pair = I(list(cfu1 = 1, cfu2 = 2)),
                                  pair = I(list(cfu1 = 1, cfu2 = 2)),
                                  pair = I(list(cfu1 = 1, cfu2 = 2)))
#make the empty datatable you want

#edited to remove the subframes
fetch <- function(x){
  names_list <- c()
  for (i in 1:ncol(xselect)){
    index <- as.character(xselect[,i])
    names_list <- c(names_list, paste(index[1], index[2], sep ="-"))
    pairwise_data_table[[i]][[1]] <- x[[index[1]]]
    pairwise_data_table[[i]][[2]] <- x[[index[2]]]
  }
  colnames(pairwise_data_table) <- names_list
  return(pairwise_data_table)
}

#populate the list
z <- OD_aggregate[["OD_adjusted"]]
pairwise_data_table <- fetch(z)
str(pairwise_data_table)
######foda-se

#+ orddom-effect-size2, ref.label="orddom-effect-size"

#+ plot-ordinal-effsize2, ref.label="plot-ordinal-effsize"

#+ plot-metric-effsize2, ref.label="plot-metric-effsize"

#+ save-agr_subsets, results="hide"
save(agr_subset, file = "agr_subset.RData") #the subsetted norm data
save(data_summary, file = "agr_subset_summary.RData") #the mean, sdev and t.test computed CI
save(orddom_ordinal, file = "agr_subset_cliffsd.RData") #cliff's d and effectsizes
save(orddom_metric, file = "agr_subset_metric.RData") #metric effectsizes

#' ##NHST
#' 
#' null hypothesis: the cfu from mutants are not different
#' 
#+ ANOVA-new2
anova_fit <- lmer(OD_adjusted ~ simple_id*coating + (1|date), data = agr_subset, REML = TRUE) #making the mixed effects model 2 factors one random error

#+ ANOVA-new-stats2, ref.label="ANOVA-new-stats"

#+ model-based-effect-sizes2, ref.label="model-based-effect-sizes"

#+ Tukey-test2, ref.label="Tukey-test"

#non parametric tests
#+ Kruskall-Wallis2
kruskal.test(OD_adjusted ~ simple_id, data = agr_subset) #Kruskall-Wallis test for a group with stochastic dominance

#+ Post-hoc-pairwise2, message=FALSE
wilcox_tests <- pairwise.wilcox.test(agr_subset$OD_adjusted, agr_subset$simple_id, p.adjust.method = 'bonferroni') #pairwise comparison if Kruskall-Wallis is rejected using Wilcoxon distrubution
wilcox_tests #print wilcox tests
dunns_test <- dunn.test(agr_subset$OD_adjusted, agr_subset$simple_id, method = 'bonferroni') #pairwise comparison if Kruskall-Wallis is rejected using Dunn's Z statistic
matrix(c(dunns_test$Z, dunns_test$P, dunns_test$P.adjusted), 3, 3, dimnames = list(dunns_test$comparisons, c("Z", "p-value", "p-adjusted")))

#' #Complemented mutants: atl
#' 
#+ subset-data-atl, results="hide"
atl_subset_index <- unique(norm_data[simple_id == "atl_pos1_empty" | simple_id == "atlA_pos1_atl", .(date, plate)])
atl_subset <- norm_data[agr_subset_index][simple_id == "atl_pos1_empty" | simple_id == "atl_pos1_atl" | simple_id == "wt_NA"]
atl_subset$simple_id <- droplevels(atl_subset$simple_id)
OD_aggregate <- aggregate(cbind(OD, OD_adjusted, OD_log, OD_scale) ~ simple_id,
                          data = atl_subset,
                          print) #the OD for each pariwise combinations, also required for Cliff's DELTA
xselect <- combn(OD_aggregate[["simple_id"]], 2)
names(OD_aggregate$OD_adjusted) <- OD_aggregate$simple_id

#+ effect-size-atl
#rewrite effect size in data.table
data_summary <- atl_subset[simple_id != "wt_NA", summarise(.SD,
                                                           mean = mean(OD_adjusted),
                                                           sdev = sd(OD_adjusted),
                                                           ci_lower = t.test(OD_adjusted)$conf.int[1],
                                                           ci_upper = t.test(OD_adjusted)$conf.int[2]), by = .(simple_id)]
data_summary <- rbindlist(list(data_summary, 
                               atl_subset[simple_id == "wt_NA", summarise(.SD, mean = mean(OD_adjusted),
                                                                          sdev = sd(OD_adjusted)), by = .(simple_id)]),
                          fill = TRUE)

#+ effect-size-plot4, ref.label="effect-size-plot"

#+foda-se-amazing2, ref.label="foda-se-amazing"

#+ orddom-effect-size3, ref.label="orddom-effect-size"

#+ save-atl_subsets, results="hide"
save(atl_subset, file = "atl_subset.RData") #the subsetted norm data
save(data_summary, file = "atl_subset_summary.RData") #the mean, sdev and t.test computed CI
save(orddom_ordinal, file = "atl_subset_cliffsd.RData") #cliff's d and effectsizes
save(orddom_metric, file = "atl_subset_metric.RData") #metric effectsizes

#+ plot-ordinal-effsize3, ref.label="plot-ordinal-effsize"

#+ plot-metric-effsize3, ref.label="plot-metric-effsize"

#' ##NHST
#' 
#' null hypothesis: the cfu from mutants are not different
#' 
#+ ANOVA-new3
anova_fit <- lmer(OD_adjusted ~ simple_id*coating + (1|date), data = atl_subset, REML = TRUE) #making the mixed effects model 2 factors one random error

#+ ANOVA-new-stats3, ref.label="ANOVA-new-stats"

#+ model-based-effect-sizes3, ref.label="model-based-effect-sizes"

#+ Tukey-test3, ref.label="Tukey-test"

#non parametric tests
#+ Kruskall-Wallis3
kruskal.test(OD_adjusted ~ simple_id, data = atl_subset) #Kruskall-Wallis test for a group with stochastic dominance

#+ Post-hoc-pairwise3, message=FALSE
wilcox_tests <- pairwise.wilcox.test(atl_subset$OD_adjusted, atl_subset$simple_id, p.adjust.method = 'bonferroni') #pairwise comparison if Kruskall-Wallis is rejected using Wilcoxon distrubution
wilcox_tests #print wilcox tests
dunns_test <- dunn.test(atl_subset$OD_adjusted, atl_subset$simple_id, method = 'bonferroni') #pairwise comparison if Kruskall-Wallis is rejected using Dunn's Z statistic
matrix(c(dunns_test$Z, dunns_test$P, dunns_test$P.adjusted), 3, 3, dimnames = list(dunns_test$comparisons, c("Z", "p-value", "p-adjusted")))
