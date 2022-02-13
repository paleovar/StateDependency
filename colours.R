COLS <- list()
COLS["LGM"] <- "#DC267F"
COLS["PI"] <- "grey30"
COLS <- lapply(COLS, function(x){adjustcolor(x, alpha.f = 0.9)})

COLS_a <- lapply(COLS, function(x){adjustcolor(x, alpha.f = 0.2)})
  
pal <- c(
  "#009E73", #green # "tropics", 
  "#56B4E9", #blue "mid_lats_n", 
  "#E69F00", #light orange "mid_lats_s", 
  "#0072B2",# dark blue "hight_lats_n", 
  "#D55E00" #dark orange "hight_lats_s"
) 

pal <-  adjustcolor(pal, 0.9)

pal_a <- adjustcolor(pal, 0.2)

#barplot(seq(1,(length(COLS)+length(pal))), col=c(unlist(COLS), pal))
#barplot(seq(1,(length(COLS_a))), col=c(unlist(COLS_a)))

IPCCColPal <- list("BlPu"=c(rgb(237,248,251,maxColorValue=255),
                            rgb(179,205,227,maxColorValue=255),
                            rgb(140,150,198,maxColorValue=255),
                            rgb(136, 86,167,maxColorValue=255),
                            rgb(129, 15,124,maxColorValue=255)),
                   "YlRd"=c(rgb(255,255,178,maxColorValue=255),
                            rgb(254,204, 92,maxColorValue=255),
                            rgb(253,141, 60,maxColorValue=255),
                            rgb(240, 59, 32,maxColorValue=255),
                            rgb(189,  0, 38,maxColorValue=255)),
                   "LineShade"=c(rgb(91,  174, 178,maxColorValue=255),
                                 rgb(204, 174, 113,maxColorValue=255),
                                 rgb(191, 191, 191,maxColorValue=255),
                                 rgb(67,  147, 195,maxColorValue=255),
                                 rgb(223, 237, 195,maxColorValue=255)),
                   "Precip"=c(rgb(84,48,5,maxColorValue = 255),
                              rgb(140,81,10,maxColorValue = 255),
                              rgb(191,129,45,maxColorValue = 255),
                              rgb(223,194,125,maxColorValue = 255),
                              rgb(246,232,195,maxColorValue = 255),
                              rgb(245,245,245,maxColorValue = 255),
                              rgb(199,234,229,maxColorValue = 255),
                              rgb(128,205,193,maxColorValue = 255),
                              rgb(53,151,143,maxColorValue = 255),
                              rgb(1,102,94,maxColorValue = 255),
                              rgb(0,60,48,maxColorValue = 255)),
                   "Temp"=c(rgb(103,0,31,maxColorValue = 255),
                            rgb(178,24,43,maxColorValue = 255),
                            rgb(214,96,77,maxColorValue = 255),
                            rgb(244,165,130,maxColorValue = 255),
                            rgb(253,219,199,maxColorValue = 255),
                            rgb(247,247,247,maxColorValue = 255),
                            rgb(209,229,240,maxColorValue = 255),
                            rgb(146,197,222,maxColorValue = 255),
                            rgb(67,147,195,maxColorValue = 255),
                            rgb(33,102,172,maxColorValue = 255),
                            rgb(5,48,97,maxColorValue = 255)))