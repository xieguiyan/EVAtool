
# library --------------------------------------------------------------------

library(magrittr)
library(ggplot2)


# load data ---------------------------------------------------------------

rna_type <- readr::read_tsv('tmp_result/SRR8185773.stat', comment = "#", col_names = c('type', 'count','percentage'))


# plot --------------------------------------------------------------------

.ncrna_type <- ggplot(data = rna_type, mapping = aes(x = reorder(type, count), y = count, fill = type)) +
  geom_bar(stat = 'identity')+
  geom_text(mapping = aes(label = percentage), vjust= -0.5)+
  scale_fill_brewer(palette = 'Set2', name = 'type')+
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
    x = 'ncRNA type',
    y = 'Reads count',
    title = 'Number of ncRNA in different types')+
  scale_y_log10()


# save image --------------------------------------------------------------

width <- 10
height <- 7
png_name ='tmp_result/distribution_of_ncRNA_type.png'
ggsave(filename = png_name, plot = .ncrna_type, device = "png", width = width, height = height)
pdf_name = 'tmp_result/distribution_of_ncRNA_type.pdf'
ggsave(filename = pdf_name, plot = .ncrna_type, device = "pdf", width = width, height = height)


# imge --------------------------------------------------------------------

save.image('/workspace/xiegy/EVAtool/data/rds/02-reads-type-distribution.rds')
load('/workspace/xiegy/EVAtool/data/rds/02-reads-type-distribution.rds')
