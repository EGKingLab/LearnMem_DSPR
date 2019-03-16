font_size <- 10
my_theme <- theme(
  axis.text = element_text(size = font_size),
  axis.title = element_text(size = font_size),
  legend.text = element_text(size = font_size),
  legend.title=element_text(size = font_size),
  plot.title = element_text(size = font_size + 1))

font_size <- 8
my_theme_sm <- theme(
  axis.text = element_text(size = font_size),
  axis.title = element_text(size = font_size),
  legend.text = element_text(size = font_size),
  legend.title=element_text(size = font_size),
  plot.title = element_text(size = font_size + 1))



library(scales)
LM_cols <- hue_pal()(2)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}