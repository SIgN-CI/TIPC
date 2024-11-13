#' Clustering of TIPC metrics
#'
#' Clustering of TIPC metrics using
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#'
#' @param root_dir A directory path containing the TIPC_metrics, i.e. normalized
#'   TIPC metrics of all 5 directions output from
#'   \code{\link[TIPC]{normalize_metrics}}; nrow = total. of tumors, ncol = 6
#'   (TIPC metrics) x 5 (directions).
#' @param min_k An integer indicating the minimum number of clusters.
#' @param max_k An integer indicating the maximum number of clusters.
#' @param distance A character value indicating the distance method used in
#' \code{\link[ConsensusClusterPlus]{ConsensusClusterPlus}}
#' @param seed An integer seed for randomizing case ordering and for ConsensusClusterPlus.
#' @param output_bnm A character string appended to output folder name;
#'   sub-folders are created for different k from min_k to max_k.
#' @export
#' @importFrom grDevices pdf
#' @importFrom utils write.csv
#' @importFrom grDevices rainbow
#' @importFrom graphics hist legend lines
#' @importFrom utils assignInNamespace
consensus_clustering <- function(min_k = 2, max_k = 6,distance='pearson',
                                 root_dir = NULL, output_bnm = 'test', seed = 999){

  ## modified CDF function to output delta values
  CDF <- function (ml, breaks = 100) {
    plot(c(0), xlim = c(0, 1), ylim = c(0, 1), col = "white",
         bg = "white", xlab = "consensus index", ylab = "CDF",
         main = "consensus CDF", las = 2)
    k = length(ml)
    this_colors = rainbow(k - 1)
    areaK = c()
    for (i in 2:length(ml)) {
      v = ConsensusClusterPlus:::triangle(ml[[i]], mode = 1)
      h = hist(v, plot = FALSE, breaks = seq(0, 1, by = 1/breaks))
      h$counts = cumsum(h$counts)/sum(h$counts)
      thisArea = 0
      for (bi in 1:(length(h$breaks) - 1)) {
        thisArea = thisArea + h$counts[bi] * (h$breaks[bi +
                                                         1] - h$breaks[bi])
        bi = bi + 1
      }
      areaK = c(areaK, thisArea)
      lines(h$mids, h$counts, col = this_colors[i - 1], lwd = 2,
            type = "l")
    }
    legend(0.8, 0.5, legend = paste(rep("", k - 1), seq(2, k,
                                                        by = 1), sep = ""), fill = this_colors)
    deltaK = areaK[1]
    for (i in 2:(length(areaK))) {
      deltaK = c(deltaK, (areaK[i] - areaK[i - 1])/areaK[i -
                                                           1])
    }
    save(deltaK, file = file.path(res_subdir,'deltaK.RData'))
    par(mar = c(5, 8, 3, 3))
    plot(1 + (1:length(deltaK)), y = deltaK, xlab = "k", ylab = "Relative change in\narea under CDF curve",
         main = "", type = "b", cex.main=2, cex.lab=2, cex.axis=2)
  }

  assignInNamespace("CDF",CDF,ns="ConsensusClusterPlus")

  ## ======================
  ## root dir check
  ## ======================
  if(is.null(root_dir)) stop('No root directory is provided!\n')
  ## ======================
  ## load TIPC_metrics
  ## ======================
  TIPC_metrics_holder <- load(file.path(root_dir,'TIPC_metrics.Rda'))
  TIPC_metrics = get(TIPC_metrics_holder)

  ## ======================
  ## create output directory
  ## ======================
  output_dir_bnm <- paste0(output_bnm)
  res_subdir <- file.path(root_dir, output_dir_bnm)
  dir.create(res_subdir)

  ## ======================
  ## data scaling
  ## ======================
  df <- as.data.frame(scale(TIPC_metrics))

  ## ======================
  ## randomize sample ordering
  ## ======================
  set.seed(seed)
  rand_order <- sample(x = nrow(df), size = nrow(df))
  df <- df[rand_order,]

  ## ======================
  ## exclude NAs columns
  ## ======================
  NA_spatial_param <- apply(df, MARGIN = 2, FUN=function(z){sum(is.na(z))})
  df <- df[,NA_spatial_param==0]
  ## ======================
  ## caling ConsensusClusterPlus
  ## ======================
  clustering_res = ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(t(df)),maxK=max_k,reps=50,
                                                              pItem=0.8,pFeature=1,
                                                              title= res_subdir,
                                                              distance=distance,seed=seed,plot="pdf")

  setwd(res_subdir)
  icl = ConsensusClusterPlus::calcICL(clustering_res,title='cluster_item_consensus_plots',plot="pdf")

  for (k in c(min_k:max_k)){
    res_k <- data.frame('tumor'=names(clustering_res[[k]][["consensusClass"]]),
                        'cluster_no'=clustering_res[[k]][["consensusClass"]])

    setwd(res_subdir)
    sub_folderNm <- paste0('k',k)
    dir.create(sub_folderNm)
    setwd(sub_folderNm)
    write.csv(res_k, file=paste0('cluster_no_k',k,'.csv'), row.names = FALSE)

  }


}
