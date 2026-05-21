## MERS posterior predictive check and probability mass contours.

rm(list = ls())

code_dir <- "codes"
source(file.path(code_dir, "load_data.R"))
source(file.path(code_dir, "funcs_util.R"))
source(file.path(code_dir, "funcs_solveH.R"))

seed <- 1

mass_levels   <- c(0.50, 0.80, 0.95) 
out_plot_pdf    <- "./results/Fig2_ppc_mass_contour.pdf"
out_summary_csv <- "./results/ppc_panel_summary.csv"

draws <- read.csv("./results/posterior_draws.csv", check.names = FALSE)

inv_logit <- function(x) 1 / (1 + exp(-x))

get_draw_params <- function(df_row){
  r     <- exp(as.numeric(df_row[["log_r"]]))
  k     <- exp(as.numeric(df_row[["log_k"]]))
  S     <- inv_logit(as.numeric(df_row[["logit_S"]]))
  fracC <- inv_logit(as.numeric(df_row[["logit_fracC"]]))
  
  list(r=r, k=k, S=S, fracC=fracC)
}

## 
M <- matrix(c(4.3, 3.0,
              1.3, 2.7), nrow = 2, byrow = TRUE)
M <- M / max(eigen(M)$values)

## 
dmax <- max(df_summ$n_reports, na.rm = TRUE)
N <- max(2^5, 2^(ceiling(log2(dmax)) + 1))

cat("dmax =", dmax, "  N =", N, "\n")

## 
set.seed(seed)
n_all <- nrow(draws)
idx_use <- 1:n_all

## accumulate mean pmf matrices ---
pmf1_sum <- matrix(0, nrow = dmax + 1, ncol = dmax + 1)  # index = child
pmf2_sum <- matrix(0, nrow = dmax + 1, ncol = dmax + 1)  # index = adult
pmfM_sum <- matrix(0, nrow = dmax + 1, ncol = dmax + 1)  # mixture (for index unknown)

mass1_vec <- numeric(0)
mass2_vec <- numeric(0)
massM_vec <- numeric(0)

ok_count <- 0L

for(tt in seq_along(idx_use)){
  i <- idx_use[tt]
  par <- get_draw_params(draws[i, , drop = FALSE])
  
  r     <- par$r
  k     <- par$k
  S     <- par$S
  fracC <- par$fracC
  
  MM <- diag(c(1, S)) %*% M / max(eigen(diag(c(1, S)) %*% M)$values)
  ngm <- r * MM
  
  pars <- offspring_params(ngm)
  R <- pars$R
  p_list <- pars$p_list
  
  ## extinction prob
  if(r <= 1){
    Q <- c(1,1)
  } else {
    Q <- extinction_prob(R, p_list, k)
  }
  
  ## radius for solveH
  rad <- choose_rad_solveH(Q, R, p_list, k,
                           rad_cap = 0.95,
                           target_rho = 0.975)
  
  ## compute FFT tables
  tables <- tryCatch(
    precompute_tables(N, R, p_list, k, rad),
    error = function(e) NULL
  )
  if(is.null(tables)) next
  
  ## unconditioned pmf
  pmf1_u <- final_size_table_uncond(dmax, index = 1, tables = tables)
  pmf2_u <- final_size_table_uncond(dmax, index = 2, tables = tables)
  
  pmf1 <- pmf1_u
  pmf2 <- pmf2_u
  pmfM <- fracC * pmf1_u + (1 - fracC) * pmf2_u
  
  pmf1_sum <- pmf1_sum + pmf1
  pmf2_sum <- pmf2_sum + pmf2
  pmfM_sum <- pmfM_sum + pmfM
  
  mass1_vec <- c(mass1_vec, sum(pmf1))
  mass2_vec <- c(mass2_vec, sum(pmf2))
  massM_vec <- c(massM_vec, sum(pmfM))
  
  ok_count <- ok_count + 1L
  
  if(tt %% 50 == 0){
    cat("processed", tt, "/", length(idx_use), " draws (ok =", ok_count, ")\n")
  }
}

if(ok_count == 0) stop("failure in precompute_tables().")

pmf1_mean <- pmf1_sum / ok_count
pmf2_mean <- pmf2_sum / ok_count
pmfM_mean <- pmfM_sum / ok_count


# Threshold calculation for HPD (mass) regions
hpd_thresh_unnormalized <- function(pmf, alpha = 0.95){
  v <- as.vector(pmf)
  v <- v[is.finite(v) & v > 0]
  if(length(v) == 0) return(NA_real_)
  
  o <- order(v, decreasing = TRUE)
  cs <- cumsum(v[o])
  target <- alpha * sum(v[o])  # Target mass relative to the finite portion
  v[o][which(cs >= target)[1]]
}

get_contour_thresholds <- function(pmf, mass_levels){
  sapply(mass_levels, function(a) hpd_thresh_unnormalized(pmf, a))
}

panel_summary_df <- data.frame(
  panel = c("child_index", "adult_index", "index_unknown"),
  pmf_sum = c(sum(pmf1_mean), sum(pmf2_mean), sum(pmfM_mean)),
  mass50 = c(
    hpd_thresh_unnormalized(pmf1_mean, 0.50),
    hpd_thresh_unnormalized(pmf2_mean, 0.50),
    hpd_thresh_unnormalized(pmfM_mean, 0.50)
  ),
  mass80 = c(
    hpd_thresh_unnormalized(pmf1_mean, 0.80),
    hpd_thresh_unnormalized(pmf2_mean, 0.80),
    hpd_thresh_unnormalized(pmfM_mean, 0.80)
  ),
  mass95 = c(
    hpd_thresh_unnormalized(pmf1_mean, 0.95),
    hpd_thresh_unnormalized(pmf2_mean, 0.95),
    hpd_thresh_unnormalized(pmfM_mean, 0.95)
  ),
  ok_draws = ok_count,
  total_draws = n_all,
  stringsAsFactors = FALSE
)

write.csv(panel_summary_df, out_summary_csv, row.names = FALSE)

df_obs <- df_summ
df_obs$idx_unknown <- df_obs$n_reports - df_obs$n_child - df_obs$n_adult

complete_child <- subset(df_obs, index_category == "child" & idx_unknown == 0)
complete_adult <- subset(df_obs, index_category == "adult" & idx_unknown == 0)

partial_child  <- subset(df_obs, index_category == "child" & idx_unknown > 0)
partial_adult  <- subset(df_obs, index_category == "adult" & idx_unknown > 0)

noindex_all    <- subset(df_obs, index_category == "NR")   # index unknown (may have idx_unknown 0 or >0)

############################################################
## ggplot version: continuous heatmap + mass contours + overlay
############################################################

library(ggplot2)

eps <- 1e-30


# Convert PMF matrix to data.frame for ggplot ---
pmf_to_df <- function(pmf, min_C = 0, min_A = 0){
  dmax <- nrow(pmf) - 1
  df <- expand.grid(n_child = 0:dmax, n_adult = 0:dmax)
  df$p <- as.vector(pmf)
  
  # Calculate true log probability for accurate contour boundaries
  df$logp_true <- log10(df$p + eps)
  
  # Clip extremely small probabilities to prevent washing out the color scale
  min_p <- 1e-6
  df$p_plot <- pmax(df$p, min_p)
  
  # Set structurally impossible regions to NA to render as strictly blank (white)
  is_impossible <- (df$n_child < min_C) | (df$n_adult < min_A) | (df$n_child == 0 & df$n_adult == 0)
  df$p_plot[is_impossible] <- NA
  
  df
}

# Convert points and segments
make_pts_df <- function(df, min_C = 0, min_A = 0){
  if(is.null(df) || nrow(df) == 0) return(NULL)
  
  res <- data.frame(
    x = df$n_child,
    y = df$n_adult
  )
  # Filter points that meet minimum required constraints
  keep <- res$x >= min_C & res$y >= min_A & !(res$x == 0 & res$y == 0)
  res <- res[keep, , drop = FALSE]
  
  if(nrow(res) == 0) return(NULL)
  return(res)
}

make_seg_df <- function(df, min_C = 0, min_A = 0){
  if(is.null(df) || nrow(df) == 0) return(NULL)
  
  res <- data.frame(
    x0 = df$n_child,
    y0 = df$n_adult + df$idx_unknown,
    x1 = df$n_child + df$idx_unknown,
    y1 = df$n_adult
  )
  
  # Total counts are constant along the segment
  total <- res$x0 + res$y0
  
  # Clip by minimum Child condition (x >= min_C)
  res$x0 <- pmax(res$x0, min_C)
  res$y0 <- total - res$x0
  
  # Clip by minimum Adult condition (y >= min_A)
  res$y1 <- pmax(res$y1, min_A)
  res$x1 <- total - res$y1
  
  # Verify valid segments (start <= end)
  keep_valid <- res$x0 <= res$x1 & total > 0
  res <- res[keep_valid, , drop = FALSE]
  
  if(nrow(res) == 0) return(NULL)
  return(res)
}

# Plotting function 
plot_pmf_gg <- function(pmf, title,
                        pts_complete = NULL,
                        seg_partial  = NULL,
                        seg_noindex  = NULL,
                        min_C = 0,
                        min_A = 0){
  
  dfp <- pmf_to_df(pmf, min_C, min_A)
  
  thr <- get_contour_thresholds(pmf, mass_levels)
  thr <- pmax(thr, eps)
  thr_log <- log10(thr)
  
  g <- ggplot(dfp, aes(x = n_child, y = n_adult)) +
    geom_raster(aes(fill = p_plot), interpolate = FALSE) +
    geom_contour(
      aes(z = logp_true),
      breaks = thr_log,
      colour = "black",
      linewidth = 0.5
    ) +
    coord_equal(expand = FALSE) +
    scale_fill_gradient(
      trans = "log10",
      low  = "#9ecae1",
      high = "#081d58",
      name = "Probability",
      limits = c(1e-6, 1),
      breaks = c(1e-6, 1e-4, 1e-2, 1),
      labels = expression(10^{-6}, 10^{-4}, 10^{-2}, 1),
      na.value = "white"
    ) +
    labs(
      title = title,
      x = "C (child)",
      y = "A (adult)"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
  
  obs_color <- "red"
  
  add_segments_and_points <- function(gg, seg_df, line_type = "solid") {
    if (is.null(seg_df) || nrow(seg_df) == 0) return(gg)
    
    is_pt <- seg_df$x0 == seg_df$x1 & seg_df$y0 == seg_df$y1
    
    if (any(is_pt)) {
      gg <- gg + geom_point(
        data = seg_df[is_pt, ],
        aes(x = x0, y = y0),
        colour = obs_color,
        shape = 16,
        size = 2.0,
        inherit.aes = FALSE
      )
    }
    
    if (any(!is_pt)) {
      gg <- gg + geom_segment(
        data = seg_df[!is_pt, ],
        aes(x = x0, y = y0, xend = x1, yend = y1),
        colour = obs_color,
        linewidth = 0.7,
        linetype = line_type,
        lineend = "round",
        inherit.aes = FALSE
      )
    }
    
    gg
  }
  
  pts <- make_pts_df(pts_complete, min_C, min_A)
  if (!is.null(pts) && nrow(pts) > 0) {
    g <- g + geom_point(
      data = pts,
      aes(x = x, y = y),
      colour = obs_color,
      shape = 16,
      size = 2.0,
      inherit.aes = FALSE
    )
  }
  
  seg1 <- make_seg_df(seg_partial, min_C, min_A)
  g <- add_segments_and_points(g, seg1, line_type = "solid")
  
  seg2 <- make_seg_df(seg_noindex, min_C, min_A)
  g <- add_segments_and_points(g, seg2, line_type = "dashed")
  
  g
}

p_child <- plot_pmf_gg(
  pmf1_mean,
  title = "Index = child",
  pts_complete = complete_child,
  seg_partial  = partial_child,
  seg_noindex  = noindex_all,
  min_C = 1,
  min_A = 0
)

p_adult <- plot_pmf_gg(
  pmf2_mean,
  title = "Index = adult",
  pts_complete = complete_adult,
  seg_partial  = partial_adult,
  seg_noindex  = noindex_all,
  min_C = 0,
  min_A = 1
)

p_mix <- plot_pmf_gg(
  pmfM_mean,
  title = "Index unknown",
  pts_complete = NULL,
  seg_partial  = NULL,
  seg_noindex  = noindex_all,
  min_C = 0,
  min_A = 0
)

library(patchwork)
p_all <- p_child + p_adult + p_mix + patchwork::plot_layout(ncol = 3)
ggsave(out_plot_pdf, p_all, width = 16, height = 5)
