# > head(DEG)
#              p_val avg_log2FC pct.1 pct.2     p_val_adj     cluster  gene label
# RPS12 1.273332e-143  0.7298951 1.000 0.991 1.746248e-139 Naive CD4 T RPS12    Up
# RPS6  6.817653e-143  0.6870694 1.000 0.995 9.349729e-139 Naive CD4 T  RPS6    Up
library(dplyr)
library(ggrepel)

color.pals <- c(
  "#DC143C", "#0000FF", "#20B2AA", "#FFA500", "#9370DB", "#98FB98", "#F08080", "#1E90FF", "#7CFC00", "#FFFF00",
  "#808000", "#FF00FF", "#FA8072", "#7B68EE", "#9400D3", "#800080", "#A0522D", "#D2B48C", "#D2691E", "#87CEEB", "#40E0D0", "#5F9EA0",
  "#FF1493", "#0000CD", "#008B8B", "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", "#4682B4", "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE",
  "#FF6347", "#6A5ACD", "#9932CC", "#8B008B", "#8B4513", "#DEB887"
)

#
#' multi volcano plot for scRNA-seq
#' @version 0.2 change legend order
#' @version 0.3 add max_overlaps for annotation
#'
#' @param dat Seurat FindAllMarkers returns, must set only.pos = F;
#' @param color.arr color list, default same as Seurat
#' @param onlyAnnotateUp only annote gene symbols for up genes
#' @param log2Foldchang threshold for annotation
#' @param adjp  threshold for annotation
#' @param top_marker gene number for annotation
#' @param max_overlaps annotation label overlapping
#'
#' @return ggplot2 obj
#' @export
#'
#' @examples
multiVolcanoPlot <- function(dat, color.arr = NULL, onlyAnnotateUp = T,
                             log2Foldchang = 0.58, adjp = 0.05, top_marker = 5,
                             max_overlaps = 10, width = 0.9) {
  if (is.null(color.arr)) {
    len <- length(unique(dat$cluster))
    color.arr <- scales::hue_pal()(len)
  }

  dat.plot <- dat %>% mutate(
    "significance" = case_when(
      p_val_adj < adjp & avg_log2FC >= log2Foldchang ~ "Up",
      p_val_adj < adjp & avg_log2FC <= -log2Foldchang ~ "Down",
      TRUE ~ "None"
    )
  )
  head(dat.plot)
  tbl <- table(dat.plot$significance)
  print(tbl)
  background.dat <- data.frame(
    dat.plot %>% group_by(cluster) %>% filter(avg_log2FC > 0) %>% summarise("y.localup" = max(avg_log2FC)),
    dat.plot %>% group_by(cluster) %>% filter(avg_log2FC <= 0) %>% summarise("y.localdown" = min(avg_log2FC)),
    x.local = seq(1:length(unique(dat.plot$cluster)))
  ) %>% select(-cluster.1)
  x.number <- background.dat %>% select(cluster, x.local)
  dat.plot <- dat.plot %>% left_join(x.number, by = "cluster")
  dat.marked.up <- dat.plot %>%
    filter(significance == "Up") %>%
    group_by(cluster) %>%
    arrange(-avg_log2FC) %>%
    top_n(top_marker, abs(avg_log2FC))
  dat.marked.down <- dat.plot %>%
    filter(significance == "Down") %>%
    group_by(cluster) %>%
    arrange(avg_log2FC) %>%
    top_n(top_marker, abs(avg_log2FC))
  dat.marked <- dat.marked.up %>% bind_rows(dat.marked.down)
  dat.infor <- background.dat %>%
    mutate("y.infor" = rep(0, length(cluster)))
  vol.plot <- ggplot() +
    geom_col(background.dat,
      mapping = aes(x.local, y.localup),
      fill = "grey80", alpha = 0.2, width = 0.9, just = 0.5
    ) +
    geom_col(background.dat,
      mapping = aes(x.local, y.localdown),
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
      breaks = c("Up", "None", "Down"),
      values = c("#d56e5e", "#cccccc", "#5390b5")
    ) +
    geom_tile(dat.infor,
      mapping = aes(x.local, y.infor),
      height = log2Foldchang * 1.3,
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
      label = paste0("|log2FC|>=", log2Foldchang, " & FDR<", adjp)
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
  ggsave(filename = "new/result/2.18/mult_volcano_plot.svg", plot = vol.plot, device = "svg", width = 20, height = 10)
}
# multiVolcanoPlot(DEG, color.pals)
# multiVolcanoPlot(scObj.markers.time)

data <- read.csv("new/result/case2/data.csv")
head(data)
multiVolcanoPlot(data, onlyAnnotateUp = F)
