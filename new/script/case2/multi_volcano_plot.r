library(dplyr)
library(ggrepel)

dat <- read.csv("new/result/case2/data.csv")
head(dat)

color.arr <- NULL
onlyAnnotateUp <- F
log2Foldchang <- 0.25
adjp <- 0.05
max_overlaps <- 10
width <- 0.9
if (is.null(color.arr)) {
  len <- length(unique(dat$cluster))
  color.arr <- scales::hue_pal()(len)
}

dat.plot <- dat %>% mutate(
  "significance" = case_when(
    p_val_adj < adjp & avg_log2FC >= log2Foldchang ~ "Up" ,
    TRUE ~ "None"
  )
)
head(dat.plot)
tbl <- table(dat.plot$significance)
print(tbl)
background.dat <- data.frame(
  dat.plot %>% group_by(cluster) %>% filter(avg_log2FC > 0) %>% summarise("y.localup" = max(avg_log2FC)),
  x.local = seq(1:length(unique(dat.plot$cluster)))
)
head(background.dat)
x.number <- background.dat %>% select(cluster, x.local)
dat.plot <- dat.plot %>% left_join(x.number, by = "cluster")
dat.marked.up <- dat.plot %>%
  filter(significance == "Up" & isTR == 1) %>%
  group_by(cluster) %>% 
  arrange(avg_log2FC) %>% 
  top_n(5,abs(avg_log2FC))

dat.marked <- dat.marked.up 
dat.marked[which(dat.marked$gene=="E2F1")[[1]],]
dat.infor <- background.dat %>%
  mutate("y.infor" = rep(0, length(cluster)))

vol.plot <- ggplot() +
  geom_col(background.dat,
    mapping = aes(x.local, y.localup),
    fill = "grey80", alpha = 0.2, width = 0.9, just = 0.5
  ) +
  geom_jitter(dat.plot,
    mapping = aes(x.local, avg_log2FC,
      color = significance
    ),
    size = 0.8, width = 0.4, alpha = 1
  ) +
  scale_color_manual(
    name = "significance",
    breaks = c("Up", "None"),
    values = c("#d56e5e", "#cccccc", "#5390b5")
  ) +
  geom_tile(dat.infor,
    mapping = aes(x.local, y.infor),
    height = log2Foldchang * 2,
    fill = color.arr[1:length(unique(dat.plot$cluster))],
    alpha = 0.5,
    width = width
  ) +
  labs(x = NULL, y = "log2 Fold change") +
  geom_text(dat.infor, mapping = aes(x.local, y.infor, label = cluster)) +
  ggrepel::geom_label_repel(
    data = if (onlyAnnotateUp) dat.marked.up else dat.marked,
    mapping = aes(x = x.local, y = avg_log2FC, label = gene),
    force = 2,
    max.overlaps = max_overlaps,
    label.size = 0,
    fill = "#00000000",
    seed = 233,
    min.segment.length = 0,
    force_pull = 2,
    box.padding = 0.1,
    segment.linetype = 3,
    hjust = 0.5
  ) +
  annotate("text",
    x = 1.5, y = max(background.dat$y.localup) + 1,
    label = paste0("log2FC>", log2Foldchang, " & FDR<", adjp)
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(size = 13, color = "black"),
    axis.text = element_text(size = 15, color = "black"),
    axis.line.y = element_line(color = "black", size = 0.8),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.background = element_blank(),
    legend.box = "horizontal",
    legend.position = c(0.13, 0.77), legend.justification = c(1, 0)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5), title = "Change")
  )
vol.plot

ggsave(filename = "new/result/case2/mult_volcano_plot.svg", plot = vol.plot, device = "svg", width = 16, height = 6)
