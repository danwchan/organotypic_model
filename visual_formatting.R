
#my preferred theme at the moment
theme_mod <- theme_bw() +
  theme(text = element_text(size = 16),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())

#and some formatting options too, orginally made to play with xtable
italic <- function(x){
  paste0('{\\emph{', x, '}}')
}

bold <- function(x){
  paste0('{\\bfseries ', x, '}')
}