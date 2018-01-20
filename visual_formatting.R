
#my preferred theme at the moment
theme_mod <- theme_bw() +
  theme(text = element_text(size = 22),
        axis.title.x = element_text(margin = unit(c(20,0,0,0), "pt")),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())

#and some formatting options too, orginally made to play with xtable
italic <- function(x){
  paste0('{\\emph{', x, '}}')
}

bold <- function(x){
  paste0('{\\bfseries ', x, '}')
}