## Demo 1: proposed method versus Lagrange-Good calculation.

rm(list = ls())

code_dir <-"codes"
source(file.path(code_dir, "funcs_util.R"))
source(file.path(code_dir, "funcs_solveH.R"))
source(file.path(code_dir, "funcs_demo1.R"))

compare_one_dmax <- function(res_LG, res_H, dmax, index) {
  key <- paste0("dmax", dmax)
  A <- res_LG$tables[[key]][[paste0("tab_index", index)]]
  B <- res_H$tables[[key]][[paste0("tab_index", index)]]
  D <- A - B
  
  c(
    max_abs = max(abs(D)),
    sum_LG = sum(A),
    sum_H = sum(B),
    target = res_LG$Q[index]
  )
}

run_compare_one_scenario <- function(rho, r1, p11, p22, k,
                                     dmaxs = c(5, 10, 20, 40),
                                     eps_rad_LG = 1e-2,
                                     rad_cap_LG = 0.9) {
  res_LG <- run_scenario_LG(
    rho = rho,
    r1 = r1,
    p11 = p11,
    p22 = p22,
    k = k,
    dmaxs = dmaxs,
    eps_rad = eps_rad_LG,
    rad_cap = rad_cap_LG
  )
  
  res_H <- run_stepA_scenario_solveH(
    rho = rho,
    r1 = r1,
    p11 = p11,
    p22 = p22,
    k = k,
    dmaxs = dmaxs
  )
  
  out <- NULL
  for (dmax in dmaxs) {
    for (index in 1:2) {
      v <- compare_one_dmax(res_LG, res_H, dmax, index)
      out <- rbind(out, c(dmax = dmax, index = index, v))
    }
  }
  
  out <- as.data.frame(out)
  out$rho <- rho
  out$r1 <- r1
  out$p11 <- p11
  out$p22 <- p22
  out$k <- k
  out$Q1 <- res_LG$Q[1]
  out$Q2 <- res_LG$Q[2]
  out$rad_LG1 <- res_LG$rad[1]
  out$rad_LG2 <- res_LG$rad[2]
  out$rad_H1 <- res_H$rad[1]
  out$rad_H2 <- res_H$rad[2]
  
  out$N_LG <- NA
  out$N_H <- NA
  for (ii in 1:nrow(out)) {
    key <- paste0("dmax", out$dmax[ii])
    out$N_LG[ii] <- res_LG$tables[[key]]$N
    out$N_H[ii] <- res_H$tables[[key]]$N
  }
  
  list(table = out, res_LG = res_LG, res_H = res_H)
}

run_compare_batch <- function(scen_df,
                              dmaxs = c(20, 40),
                              eps_rad_LG = 1e-2,
                              rad_cap_LG = 0.95) {
  all_tab <- NULL
  all_res <- vector("list", nrow(scen_df))
  
  for (s in 1:nrow(scen_df)) {
    ans <- run_compare_one_scenario(
      rho = scen_df$rho[s],
      r1 = scen_df$r1[s],
      p11 = scen_df$p11[s],
      p22 = scen_df$p22[s],
      k = scen_df$k[s],
      dmaxs = dmaxs,
      eps_rad_LG = eps_rad_LG,
      rad_cap_LG = rad_cap_LG
    )
    
    all_tab <- rbind(all_tab, ans$table)
    all_res[[s]] <- ans
  }
  
  list(table = all_tab, results = all_res)
}

scen_df <- expand.grid(
  rho = c(0.5, 1.25),
  r1 = c(0.25),
  p11 = c(0.25, 0.75),
  k = c(0.1, 1)
)
scen_df$p22 <- scen_df$p11

batch <- run_compare_batch(scen_df, dmaxs = c(20, 40))
batch$table

metric_df_output <- function(df, digits = 3) {
  out <- df
  num <- sapply(out, is.numeric)
  out[num] <- lapply(out[num], function(x) signif(x, digits))
  out
}

tab_dmax40 <- subset(batch$table, dmax == 40)
tab_dmax40_metric <- metric_df_output(tab_dmax40, digits = 4)
write.csv(tab_dmax40_metric, file = "./results/comparison_dmax40.csv", row.names = FALSE)
