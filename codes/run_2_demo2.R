## Demo 2: likelihood surface for a simulated two-type dataset.

rm(list = ls())

code_dir <- "codes"
source(file.path(code_dir, "funcs_util.R"))
source(file.path(code_dir, "funcs_solveH.R"))
source(file.path(code_dir, "funcs_demo2.R"))

library(ggplot2)

out_dir <- "./results/demo2_surface_paper_sub075_k010_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_index <- c(20L, 20L)
max_gen <- 200L
max_total <- 120L
max_resim <- 200L

scenario <- list(
  rho = 0.75,
  k = 0.10,
  seed = 1234L
)

## Code convention: rows are offspring types and columns are parent types.
M_contact <- matrix(
  c(4.3, 3.0,
    1.3, 2.7),
  nrow = 2,
  byrow = TRUE
)
ngm_base <- normalize_by_spectral_radius(M_contact)

opts_mle <- list(
  algorithm = "NLOPT_GN_DIRECT_L",
  maxeval = 2000,
  print_level = 0
)

lb_par <- c(log(0.05), log(1e-3))
ub_par <- c(log(5.00), log(50.0))

starts <- rbind(
  c(log(0.70), log(0.10)),
  c(log(0.70), log(0.30)),
  c(log(1.00), log(0.10)),
  c(log(0.90), log(1.00))
)

ngrid_rho <- 31L
ngrid_k <- 31L
lrt_levels <- c(0.50, 0.80, 0.95, 0.99)
pmf_mass_levels <- c(0.50, 0.80, 0.95, 0.99)

fill_cols <- c(
  "#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1",
  "#6BAED6", "#4292C6", "#2171B5", "#084594"
)

rho_true <- scenario$rho
k_true <- scenario$k
seed0 <- scenario$seed

cat("\nDemo 2: simulated dataset, likelihood surface, and PMF at MLE\n")
cat(sprintf("true rho = %.3f, true k = %.3f\n", rho_true, k_true))

sim_ans <- simulate_no_trunc(
  rho_true = rho_true,
  k_true = k_true,
  seed0 = seed0,
  n_index = n_index,
  ngm_base = ngm_base,
  max_gen = max_gen,
  max_total = max_total,
  max_resim = max_resim
)

sim_df <- sim_ans$sim_df
n_trunc <- sim_ans$n_trunc
seed_use <- sim_ans$seed
n_attempt <- sim_ans$attempt

cat("truncated clusters:", n_trunc, "/", nrow(sim_df),
    "| seed used:", seed_use,
    "| attempts:", n_attempt, "\n")

lists <- sim_df_to_complete_lists(sim_df)
listC <- lists$list_complete_C
listA <- lists$list_complete_A

grid_size <- choose_dmax_from_lists(listC, listA)
dmax_obs <- grid_size$dmax
N <- grid_size$N
cat("dmax_obs =", dmax_obs, "| N =", N, "\n")

loglik_fun <- make_loglik_uncond(
  listC = listC,
  listA = listA,
  N = N,
  ngm_base = ngm_base
)

ll_true <- loglik_safe(loglik_fun, c(log(rho_true), log(k_true)))
cat(sprintf("logLik at true parameter: %.2f\n", ll_true))

fit <- fit_mle_nloptr(
  loglik_fun = loglik_fun,
  starts = starts,
  lb = lb_par,
  ub = ub_par,
  opts = opts_mle
)

if (!isTRUE(fit$ok)) {
  stop("MLE failed for the demo 2 dataset.")
}

cat(sprintf("MLE: rho = %.3f, k = %.3f | logLik = %.2f | status = %s\n",
            fit$rho, fit$k, fit$logLik, as.character(fit$status)))

logrho_grid <- seq(-1.5, 0.5, length.out = ngrid_rho)
logk_grid <- seq(-4.0, 1.0, length.out = ngrid_k)
cat("grid evaluations:", length(logrho_grid) * length(logk_grid), "\n")

df_grid <- eval_ll_grid(logrho_grid, logk_grid, loglik_fun)

p_ll <- plot_ll_surface_2d(
  df_grid = df_grid,
  mle = fit,
  lrt_levels = lrt_levels,
  fill_cols = fill_cols
)

pmf_fit <- make_pmf_tables_at_fit(
  fit = fit,
  ngm_base = ngm_base,
  N = N,
  dmax = dmax_obs
)

p_pmf_child <- plot_pmf_overlay(
  pmf = pmf_fit$pmf_child,
  pts = listC,
  index_label = "child",
  mass_levels = pmf_mass_levels,
  fill_cols = fill_cols
)

p_pmf_adult <- plot_pmf_overlay(
  pmf = pmf_fit$pmf_adult,
  pts = listA,
  index_label = "adult",
  mass_levels = pmf_mass_levels,
  fill_cols = fill_cols
)

file_three <- file.path(out_dir, "demo2_three_panel.pdf")
file_sumcsv <- file.path(out_dir, "demo2_surface_summary.csv")

save_three_panel_pdf(
  p1 = p_ll,
  p2 = p_pmf_child,
  p3 = p_pmf_adult,
  file = file_three,
  width = 15.0,
  height = 5.0
)

summary_df <- data.frame(
  rho_true = rho_true,
  k_true = k_true,
  rho_mle = fit$rho,
  k_mle = fit$k,
  logLik_true = ll_true,
  logLik_mle = fit$logLik,
  Q1_hat = pmf_fit$Q[1],
  Q2_hat = pmf_fit$Q[2],
  dmax_obs = dmax_obs,
  N = N,
  n_truncated = n_trunc,
  seed_base = seed0,
  seed_used = seed_use,
  n_resim_attempts = n_attempt,
  stringsAsFactors = FALSE
)

write.csv(summary_df, file_sumcsv, row.names = FALSE)
print(summary_df)
