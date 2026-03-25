rm(list = ls())
source("R/utilities.R")
source("R/Data_functions.R")
library(ggplot2)
library(gridExtra)
library(grid)
Data_plot <- function(dataset, TITLE = "", interval = NULL) {
  cont_size <- nrow(dataset)
  
  # --- Base Data ---
  plot_data <- data.frame(
    x_values = rep(0:(cont_size - 1), 2),
    y_values = c(dataset[, 1], dataset[, 2]),
    legend = factor(rep(c("Group 1", "Group 2"), each = cont_size))
  )
  
  diff_data <- data.frame(
    x_values = 0:(cont_size - 1),
    y_values = dataset[, 2] - dataset[, 1],
    legend = factor(rep("Effect", cont_size))
  )
  
  combined_data <- rbind(plot_data, diff_data)
  
  # --- Colors and linetypes (use plain names for keys) ---
  color_values <- c(
    "Group 1" = "tomato",
    "Group 2" = "cadetblue",
    "Effect"  = "black",
    "D1"      = "gold",   # will map to \u211C_1 in legend
    "D1p"     = "darkviolet"   # will map to \u211C_1' in legend
  )
  
  linetype_values <- c("Group 1" = "solid", "Group 2" = "solid", "Effect" = "dotted")
  
  # --- Base plot ---
  p <- ggplot(combined_data, aes(x = x_values, y = y_values, color = legend, linetype = legend)) +
    geom_line(linewidth = 1) +
    labs(title = TITLE, x = "Domain", y = "Value") +
    scale_color_manual(
      values = color_values,
      breaks = c("Group 1", "Group 2", "Effect", "D1", "D1p"),
      labels = c("Group 1", "Group 2", "Effect", expression(D[1]), expression(D[1]^"'"))
    ) +
    scale_linetype_manual(values = linetype_values) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      axis.title.x = element_text(size = 12),
    )
  
  # --- Compute contiguous ranges where Effect != 0 ---
  effect_mask <- diff_data$y_values != 0
  runs <- rle(effect_mask)
  ends <- cumsum(runs$lengths)
  starts <- c(1, head(ends, -1) + 1)
  if (any(runs$values)) {
    intervals_effect <- data.frame(start = starts[runs$values], end = ends[runs$values])
  } else {
    intervals_effect <- data.frame(start = numeric(0), end = numeric(0))
  }
  
  # --- Horizontal bar positions below x-axis ---
  y_axis_min <- min(diff_data$y_values)
  # offset1 <- y_axis_min - (max(abs(y_axis_min), 1) * 0.15)  # relative spacing
  # offset2 <- y_axis_min - (max(abs(y_axis_min), 1) * 0.30)
   offset1 <- -20  # relative spacing
   offset2 <- -30
  
  # --- Add Effect ≠ 0 bars (mapped to "D1") ---
  if (nrow(intervals_effect) > 0) {
    df_D1 <- data.frame(
      x = intervals_effect$start - 1,
      xend = intervals_effect$end - 1,
      y = offset1,
      yend = offset1
    )
    p <- p + geom_segment(
      data = df_D1,
      inherit.aes = FALSE,
      aes(x = x, xend = xend, y = y, yend = yend, color = "D1"),
      linewidth = 3,
      lineend = "round"
    ) + guides(
      color = guide_legend(override.aes = list(
        linetype = ifelse(is.null(interval),c("solid", "solid", "dotted", "solid"),
                          c("solid", "solid", "dotted", "solid", "solid"))
      )),
      linetype = "none")
  }
  
  # --- Add user-defined intervals (mapped to "D1p") ---
  if (!is.null(interval)) {
    if (is.matrix(interval) || is.data.frame(interval)) {
      interval_df <- as.data.frame(interval, check.names = FALSE)
      names(interval_df) <- c("xmin", "xmax")
    } else if (is.vector(interval) && length(interval) == 2) {
      interval_df <- data.frame(xmin = interval[1], xmax = interval[2])
    } else {
      stop("interval must be a vector of length 2, or a matrix/data frame with 2 columns")
    }
    df_D1p <- data.frame(
      x = interval_df$xmin,
      xend = interval_df$xmax,
      y = offset2,
      yend = offset2
    )
    p <- p + geom_segment(
      data = df_D1p,
      inherit.aes = FALSE,
      aes(x = x, xend = xend, y = y, yend = yend, color = "D1p"),
      linewidth = 3,
      lineend = "round"
    ) + guides(
      color = guide_legend(override.aes = list(
        linetype = c("solid", "solid", "dotted", "solid", "solid"),
        size = c(1, 1, 1, 3, 3)
      )),
      linetype = "none")
  }
  
  return(p)
}





p1 <- Data_plot(MF_data("both")[,c(2,1)], TITLE="Real Mean Functions", interval=c(32,54)) + ylim(-30, 200)


Simulation_data <- square_pulse(start_end_pulse = c(32,54), start_height= 0,
                                pulse_height = max(MF_data("both")[,1] -
                                  MF_data("both")[,2]),
                                continuum_size=101) 
Simulation_data <- cbind(rep(0,101),Simulation_data)

p2 <- Data_plot(Simulation_data, TITLE="Simulated Mean Functions", interval= c(32,54)) + ylim(-30, 200)



# Ensure both plots have the same y-limits and breaks
y_limits <- c(-30, 200)
y_breaks <- seq(0, 200, by = 50)

p1_clean <- p1 + 
  ylab(NULL) + 
  scale_y_continuous(limits = y_limits, breaks = y_breaks) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none")

p2_clean <- p2 + 
  ylab(NULL) +
  scale_y_continuous(limits = y_limits, breaks = y_breaks) +
  theme(legend.position = "none")




# Create a legend using one of the plots
legend_plot <-   p1
# Convert the plot to a grob object
grob_legend <- ggplotGrob(legend_plot)

# Extract legend component from the gtable
legend_index <- which(sapply(grob_legend$grobs, function(x) x$name) == "guide-box")
shared_legend <- grob_legend$grobs[[legend_index]]

combined_plot <- grid.arrange(p2_clean, p1_clean, ncol = 2)

final_plot <- cowplot::plot_grid(combined_plot, cowplot::plot_grid(shared_legend),
                                 ncol = 1, rel_heights = c(0.9, 0.1))
# Display the final plot


# Add shared y-axis on the left
plot <- grid.arrange(
  final_plot,
  left = textGrob("Value", rot = 90, gp = gpar(fontsize = 12), hjust = .1)
)

# Save the final plot, enable if needed
ggsave("outputs/Figure1.jpeg", plot, width = 190, height = 70, units = "mm", dpi = 1200)


