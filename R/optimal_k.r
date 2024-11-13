#' Determine the optimal cluster number k
#'
#' Determine the optimal cluster number k based on
#' (i) first finding the shoulder point in cumulative CDF plot representing the smallest k with stability, and followed up
#' (ii) finding the largest (most granular) k in tracking plot that gives the highest information.
#'
#' @param data A data.frame containing the 6 normalized TIPC spatial parameters (in %)
#'   in all 5 directions across all hexLen output from
#'   \code{\link[TIPC]{trend_plot_hexLen}}.
#' @param optHex_dir A character string indicating the folder name corresponding the selected optimal hexLen
#' which can be determined using \code{\link[TIPC]{optimal_hexLen}}.
#' @export
#' @import infotheo
#' @importFrom plyr mapvalues
#' @return A character indicating the optimal cluster k.
optimal_k <- function(optHex_dir){

  ## ======================
  ## Identify shoulder point in cumulative CDF
  ## ======================
  delta <- get(load(file = file.path(optHex_dir,'deltaK.RData')))
  delta_df <- data.frame(k=c(1:length(delta))+1, delta=delta, stringsAsFactors = FALSE)
  # Compute differences between consecutive k
  diff_values <- diff(delta_df$delta)
  # Identify the shoulder point
  shoulder_point <- which.min(diff_values)+1
  CDF_stableK <- delta_df$k[shoulder_point]
  cat('Smallest stable cluster (shoulder point)=', CDF_stableK, '\n')
  ## ======================
  ## Identify the largest (most granular) k in tracking plot
  ## ======================
  ## collect cluster assignment from all k folders
  k_folders <- list.files(path = optHex_dir, pattern = '^k[0-9]+$')
  cluster_across_k <- c()
  ## get first k
  ii=1; ff=k_folders[ii]
  k_ff <- gsub(x=ff, pattern = 'k', replacement = '')
  cluster_ff <- read.csv(file = file.path(optHex_dir, ff, paste0('cluster_no_k',k_ff,'.csv')))
  cluster_df <- data.frame(cluster_ff)
  colnames(cluster_df) <- plyr::mapvalues(x=colnames(cluster_df), from='cluster_no', to=ff)
  cluster_across_k <- cluster_df
  ## loop over and collect all k's
  for(ff in k_folders[-1]){
    k_ff <- gsub(x=ff, pattern = 'k', replacement = '')
    cluster_ff <- read.csv(file = file.path(optHex_dir, ff, paste0('cluster_no_k',k_ff,'.csv')))
    cluster_df <- data.frame(cluster_ff)
    colnames(cluster_df) <- mapvalues(x=colnames(cluster_df), from='cluster_no', to=ff)
    cluster_across_k <- merge(cluster_across_k, cluster_df, by='tumor')

  }#end all k folders
  ## clean up the tumor ids
  rownames(cluster_across_k) <- cluster_across_k$tumor
  cluster_across_k$tumor <- NULL
  ## reorder from the data based on increasing k
  col_kID <- colnames(cluster_across_k)
  col_kID <- gsub(x=col_kID, pattern = 'k', replacement = '')
  col_kID <- as.integer(col_kID)
  reorder_col <- order(col_kID)
  cluster_across_k <- cluster_across_k[, reorder_col]
  ## Function to compute normalized mutual information
  compute_nmi <- function(x, y) {
    mi <- mutinformation(x, y)
    hx <- entropy(x)
    hy <- entropy(y)
    nmi <- mi / sqrt(hx * hy)
    return(nmi)
  }
  ## Compute NMI for each pair of columns
  n <- ncol(cluster_across_k)
  nmi_out <- c()
  mi_out <- c()
  for (i in 1:(n-1)) {
    j <- i+1
    if(j>n)next
    nmi_ij <- compute_nmi(cluster_across_k[, i], cluster_across_k[, j])
    #cat(i, ' vs ', j, ' nmi=', nmi_ij , '\n')
    nmi_out <- c(nmi_out,nmi_ij)
    names(nmi_out)[i] <- colnames(cluster_across_k)[j]
    mi_ij <- mutinformation(cluster_across_k[, i], cluster_across_k[, j])
    mi_out <- c(mi_out, mi_ij)
    names(mi_out)[i] <- paste0('k',j)
  }
  ## identify the most informative cluster
  nmi_out_df <- data.frame(k=gsub(x=names(nmi_out),pattern = 'k', replacement = ''),
                           nmi=nmi_out, stringsAsFactors = FALSE)
  nmi_out_df$k <- as.integer(nmi_out_df$k)
  ## satisify first criteria
  nmi_out_df <- nmi_out_df[nmi_out_df$k > CDF_stableK,]
  nmi_out_df$delta_nmi <- c(Inf, abs(diff(nmi_out_df$nmi)))
  #max_k <- which.max(nmi_out_df$nmi)
  #cat('Largest most informative cluster=', nmi_out_df[max_k-1,"k"], '\n')

  max_stable_k <- which.min(nmi_out_df$delta_nmi)
  cat('Largest most informative & stable cluster=', nmi_out_df[max_stable_k-2,"k"], '\n')
  return(nmi_out_df[max_stable_k-2,'k'])
}
