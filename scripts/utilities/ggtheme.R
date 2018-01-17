
######################
# Create your own ggplot2 theme or see ggthemes package
# https://cran.r-project.org/web/packages/ggthemes/vignettes/ggthemes.html
# http://sape.inf.usi.ch/quick-reference/ggplot2/themes
# The following is modified from:
# https://stackoverflow.com/questions/31404433/is-there-an-elegant-way-of-having-uniform-font-size-for-the-whole-plot-in-ggplot
# If re-using functions it is much better to copy this to a new script and run here as:
# source('my_ggplot_theme.R')
theme_my <- function(base_size = 14, base_family = "Times") {
  normal_text <- element_text(size = as.numeric(base_size), colour = "black", face = "plain")
  large_text <- element_text(size = as.numeric(base_size + 1), colour = "black", face = "plain")
  bold_text <- element_text(size = as.numeric(base_size + 1), colour = "black", face = "bold")
  axis_text <- element_text(size = as.numeric(base_size - 1), colour = "black", face = "plain")
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(legend.key = element_blank(),
          strip.background = element_blank(),
          text = normal_text,
          plot.title = bold_text,
          axis.title = large_text,
          axis.text = axis_text,
          legend.title = bold_text,
          legend.text = normal_text,
          #plot.margin = unit(c(1,1,1,1),"mm")
          plot.margin = grid::unit(c(1,1,1,1), "mm") 
          )
  }
theme_my()
######################


