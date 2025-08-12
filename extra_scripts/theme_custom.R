# library(extrafont)
# extrafont::font_import()
# loadfonts(device = "win")
sysfonts::font_add_google("Open Sans")
library(tidyverse)

theme_custom <- function () { 
  theme_bw(base_size=12, base_family="Open Sans") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="white", color=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size=18, family="Open Sans"),
      axis.text = element_text(size=14),
      legend.text = element_text(size=14),
      legend.title = element_blank()
    )
}

