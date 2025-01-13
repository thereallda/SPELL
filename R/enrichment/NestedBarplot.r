# Nested barplot for enrichment results visualization
library(tidyverse)
library(paintingr)
library(gprofiler2)
goM1 <- gost(aov_sig_ls$M.3$geneID,
     organism = 'hsapiens', exclude_iea = TRUE, 
     correction_method = 'fdr', evcodes = TRUE)
gostplot(goM1)

godat_M1 <- goM1$result %>% 
  filter(source %in% c('KEGG', "REAC")) %>% 
  group_by(source) %>% 
  slice_min(n=3, order_by = p_value) %>% 
  mutate(inter=paste(term_name, source, sep='.'))

godat_M1 %>% 
  ggplot(aes(x = weave_factors(inter), y = -log10(p_value), fill = source)) +
  geom_bar(stat = 'identity', width = 0.8, color = 'black') +
  coord_flip() +
  scale_x_discrete(guide = "axis_nested") +
  theme_minimal() +
  theme(axis.text = element_text(color='black'),
        axis.line = element_line(),
        ggh4x.axis.nestline.y = element_line(size=1),
        ggh4x.axis.nesttext.y = element_text(color='black', face='bold', size=12),
        panel.grid = element_blank(),
        legend.position = 'top'
  ) +
  scale_fill_manual(values=paint_palette("Villeneuve", length(unique(godat_M1$source)), 'continuous')) +
  scale_y_continuous(expand = expansion(0)) +
  labs(x='', title='signatures', fill='')
ggsave('results/signatures.png', width=8, height=3, bg='white')