
# library -----------------------------------------------------------------


library(magrittr)
library(ggplot2)


# set path ----------------------------------------------------------------

setwd('/workspace/xiegy/EVAtool/')

# load trimmed file -------------------------------------------------------

trimm_file <- readr::read_tsv('./tmp_result/SRR8185773.freq.stat', col_names = c('seq_len','count', 'ratio'), col_types = "cdd")


# stat --------------------------------------------------------------------
trimm_file %>% 
  dplyr::filter(seq_len != "ok") %>% 
  dplyr::arrange(count) -> .data_for_plot


# plot --------------------------------------------------------------------

ggplot(data = .data_for_plot, mapping = aes(x = seq_len, y = ratio, group = 1)) +
  geom_line()+
  theme(
    axis.text = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    plot.title = element_text(hjust = 0.5),
    legend.position = 'none',
    panel.grid = element_blank(),
    panel.background=element_rect(fill="white"),
    axis.line = element_line(size = 0.5, colour = 'black')
  ) +
  labs(
    x = 'Read length',
    y = 'Gene count ratio',
    title = 'Distribution of reads length')+
  scale_y_continuous(labels = scales::percent) -> .d_plot

# save image --------------------------------------------------------------

width <- 10
height <- 7
png_name ='tmp_result/distribution_of_read_len.png'
ggsave(filename = png_name, plot = .d_plot, device = "png", width = width, height = height)
pdf_name = 'tmp_result/distribution_of_read_len.pdf'
ggsave(filename = pdf_name, plot = .d_plot, device = "pdf", width = width, height = height)


# imge --------------------------------------------------------------------

save.image('/workspace/xiegy/EVAtool/data/rds/02-reads-len-distribution.rds')
load('/workspace/xiegy/EVAtool/data/rds/02-reads-len-distribution.rds')
  





