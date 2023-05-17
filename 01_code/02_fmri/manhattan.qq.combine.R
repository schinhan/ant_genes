#!/usr/bin/env Rscript

# ======================================
# === create combined manhattan plot ===
# ======================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=6) {
  stop(paste0('expected 6 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits = args[1] # traits="gap_gm,gap_wm,gap_gwm"
manhattanPlots = args[2] # manhattanPlots="/.../03_data_visualization/01_fMRI/Manhattan_QQ_Plots/manhattan.png"
qqPlots = args[3] # qqPlots="/.../03_data_visualization/01_fMRI/Manhattan_QQ_Plots/qqplot.png"
outputFile = args[4] # outputFile="/.../03_data_visualization/01_fMRI/Manhattan_QQ_Plots/qqplot.manhattan.png"
width = as.numeric(args[5]) # width = 12
height = as.numeric(args[6]) # height = 14

message(paste0('\n--- Settings for qqplot.manhattan Figure ---',
               '\ntraits: ', traits,
               '\nmanhattanPlots: ', manhattanPlots,
               '\nqqPlots: ', qqPlots,
               '\noutputFile: ', outputFile,
               '\nwidth: ', width,
               '\nheight: ', height,'\n'))

# attach required packages
for (pkg in c('ggpubr','ggplot2','patchwork','magick','stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# transform variables
traits = str_split(traits, ',')[[1]]
manhattanPlots = str_split(manhattanPlots, ',')[[1]]
qqPlots = str_split(qqPlots, ',')[[1]]

# load images
for (i in 1:length(traits)) {
  
  # load manhattan plots
  message(paste0('[', i, '/', length(traits), '] Loading ', manhattanPlots[i], ' and ', qqPlots[i]))
  manh = image_read(manhattanPlots[i])
  manh.width = image_info(manh)$width
  manh.height = image_info(manh)$height
  manh = ggplot() +
  background_image(manh) + coord_fixed(ratio = image_info(manh)$height/image_info(manh)$width)
  assign(paste0(traits[i],"_manh"), manh)
  
  # load qq-plots
  qq = image_read(qqPlots[i])
  qq.width = image_info(qq)$width
  qq.height = image_info(qq)$height
  qq = ggplot() +
  background_image(qq) + coord_fixed(ratio = image_info(qq)$height/image_info(qq)$width)
  assign(paste0(traits[i],"_qq"), qq)
  
  # collect plots
  tmp = manh + qq + plot_layout(widths = c(manh.width,qq.width))
  if (i == 1) { pl = tmp } else { pl = pl / tmp }
}

# set layout
if (length(traits) > 1) {
  pl = pl + plot_annotation(tag_levels = list(c(rbind(letters[1:length(traits)],letters[(1+length(traits)):(2*length(traits))])))) & # list(c(rbind(letters,' ')))
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 22))
}

# save file
png(width = width, height = height, units = "in", res = 300, filename = outputFile)
pl
invisible(dev.off())
