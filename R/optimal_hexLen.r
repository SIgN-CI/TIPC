#' Determine the optimal hexagon length
#'
#' Determine the optimal hexagon length based on
#' (i) the smallest difference in Q3-Q1 across the 6 TIPC spatial parameters, and
#' (ii) the smallest variance/sd across 5 directions within each TIPC spatial parameters
#'
#' @param data A data.frame containing the 6 normalized TIPC spatial parameters (in %)
#'   in all 5 directions across all hexLen output from
#'   \code{\link[TIPC]{trend_plot_hexLen}}.
#' @param output_dir A character string indicating the output folder name; the directory the output pdf will be saved.
#' @export
#' @importFrom grDevices pdf
#' @importFrom stats quantile
#' @importFrom dplyr summarise group_by
#' @import ggplot2
#' @return A character indicating the optimal hexLen.
optimal_hexLen <- function(data, output_dir){
  ## ======================
  ## collapse all samples by TIPC params, hexLen, and directions
  ## ======================
  mean_df <- data %>% dplyr::group_by(TIPC.cat, hexLen, direction) %>% dplyr::summarise(avg_perc = mean(counts))
  pl <- ggplot(data = mean_df, aes(x=hexLen, y=avg_perc, group=hexLen, color=TIPC.cat)) + geom_boxplot() + theme_classic()+geom_point()
  ggsave(filename = file.path(output_dir, paste0('Avg_TIPCmetric_vs_hexLen',paste0(unique(data$hexLen), collapse = '_'),'_normalized.pdf')), plot = pl)

  ## ======================
  ## collapse all directions by TIPC params and hexLen
  ## ======================
  mean_df2 <- mean_df %>% dplyr::group_by(hexLen, TIPC.cat) %>%
    summarise(sd=sd(avg_perc),max=min(avg_perc),min=min(avg_perc), avg=mean(avg_perc),
              Q3_Q1 = quantile(avg_perc, 0.75) - quantile(avg_perc, 0.25))

  ## ======================
  ## collapse all TIPC params by hexLen and compute Q3-Q1
  ## ======================
  Q3_Q1 <- mean_df2 %>% dplyr::group_by(hexLen) %>%
    dplyr::summarise(Q3_Q1 = quantile(avg, 0.75) - quantile(avg, 0.25))
  minN <- min(4, nrow(Q3_Q1))
  min_Q3_Q1 <- Q3_Q1[order(Q3_Q1$Q3_Q1, decreasing = F)[1:minN],]
  ## ======================
  ## stability/variability of hexagonal shifting in 5 directions within each TIPC spatial parameters
  ## ======================
  sd_df <- mean_df2 %>% dplyr::group_by(hexLen) %>% dplyr::summarise(avg_sd = mean(sd))
  opt_dir_sd <- sd_df[order(sd_df$avg_sd, decreasing = F)[1:minN],]

  opt_hexLen <- intersect(min_Q3_Q1$hexLen, opt_dir_sd$hexLen)
  if(length(opt_hexLen)>1){
    opt_hexLen <- intersect(Q3_Q1$hexLen[which.min(Q3_Q1$Q3_Q1)], opt_dir_sd$hexLen[which.min(opt_dir_sd$avg_sd)])
  }
   return(opt_hexLen)
}
