rm(list = ls())

code_dir <- "codes"
source(file.path(code_dir, "load_data.R"))
source(file.path(code_dir, "funcs_util.R"))
source(file.path(code_dir, "funcs_solveH.R"))

library(BayesianTools)

M <- matrix(
  c(4.3, 3.0, 1.3, 2.7), nrow = 2, byrow = TRUE)
M <- M / max(eigen(M)$values)
M

df_summ %>%
  filter(index_category == "adult" & n_reports == n_child + n_adult) %>% 
  mutate(count = 1) %>%
  group_by(n_child, n_adult) %>%
  summarize(count = sum(count)) -> list_complete_A

df_summ %>%
  filter(index_category == "child" & n_reports == n_child + n_adult) %>% 
  mutate(count = 1) %>%
  group_by(n_child, n_adult) %>%
  summarize(count = sum(count)) -> list_complete_C

df_summ %>%
  filter(index_category == "child" & n_reports != n_child + n_adult) %>% 
  mutate(count = 1) %>%
  group_by(n_child, n_adult, n_reports) %>%
  summarize(count = sum(count)) -> list_C_index_only

df_summ %>%
  filter(index_category == "adult" & n_reports != n_child + n_adult) %>%   
  mutate(count = 1) %>%
  group_by(n_child, n_adult, n_reports) %>%
  summarize(count = sum(count)) -> list_A_index_only

df_summ %>%
  filter(index_category == "NR") %>% 
  mutate(count = 1) %>%
  group_by(n_child, n_adult, n_reports) %>%
  summarize(count = sum(count)) -> list_no_index_info

n_known_child_index = sum(df_summ$index_category == "child")
n_known_adult_index = sum(df_summ$index_category == "adult")


loglik_2G = function(par){

  r <- exp(par[1])
  k <- exp(par[2])
  logit_S <- par[3]
  logit_fracC <- par[4]

  S <- 1 / (1+exp(-logit_S))
  fracC <- 1 / (1+exp(-logit_fracC))
  
  MM <- diag(c(1, S)) %*% M / max(eigen(diag(c(1, S))%*%M)$values)
  
  ngm <- r*MM
  
  offspring_params <- offspring_params(ngm)
  
  
  
  if(r <=1){
    Q <- c(1,1)
  }else{
    
    Q <- extinction_prob(offspring_params$R, offspring_params$p_list, k)
    
  }
  
  rad <- choose_rad_solveH(Q, offspring_params$R, offspring_params$p_list, k,
                             rad_cap = 0.95,
                             target_rho = 0.975)
  
  
  max_deg <- max(
    list_complete_A$n_child + list_complete_A$n_adult,
    list_complete_C$n_child + list_complete_C$n_adult,
    list_A_index_only$n_reports,
    list_C_index_only$n_reports,
    list_no_index_info$n_reports,
    na.rm = TRUE
  )
  
  N <- max(2^5, 2^(ceiling(log2(max_deg))+1))
  
  tables <- precompute_tables(N, offspring_params$R, offspring_params$p_list, k, rad)
  
  
  ### ============
  ### LIKELIHOODS
  ### ============
  
  llk_total <- 0
  
  if (nrow(list_complete_C) > 0){
    
    for (ii in 1:nrow(list_complete_C)){
      
      nC <- list_complete_C$n_child[ii]
      nA <- list_complete_C$n_adult[ii]
      p  <- coef_nm(1, nC, nA, tables) 
      
      llk_total <- llk_total + log(p) * list_complete_C$count[ii]
    }
  }
  
  if (nrow(list_complete_A) > 0){
    
    for (ii in 1:nrow(list_complete_A)){
      
      nC <- list_complete_A$n_child[ii]
      nA <- list_complete_A$n_adult[ii]
      p  <- coef_nm(2, nC, nA, tables) 
      
      llk_total <- llk_total + log(p) * list_complete_A$count[ii]
    }
  }
  
  
  if (nrow(list_C_index_only) > 0){
    
    for (ii in 1:nrow(list_C_index_only)){
      nC_obs <- list_C_index_only$n_child[ii]
      nA_obs <- list_C_index_only$n_adult[ii]
      
      nTOT  <- list_C_index_only$n_reports[ii]
      idx  <- nTOT - nC_obs - nA_obs
      
      p_tmp <- 0
      
      for (j in 0 : idx){
        p_tmp <- p_tmp + coef_nm(1, nC_obs + j, nA_obs + (idx-j), tables) 
      }
      
      llk_total <- llk_total + log(p_tmp) * list_C_index_only$count[ii]
    }
  }
  
  if (nrow(list_A_index_only) > 0){
    for (ii in 1:nrow(list_A_index_only)){
      nC_obs <- list_A_index_only$n_child[ii]
      nA_obs <- list_A_index_only$n_adult[ii]
      
      nTOT<- list_A_index_only$n_reports[ii]
      idx  <- nTOT - nC_obs - nA_obs
      
      p_tmp <- 0
      
      for (j in 0 : idx){
        p_tmp <- p_tmp + coef_nm(2, nC_obs + j, nA_obs + (idx-j), tables) 
      }
      
      llk_total <- llk_total + log(p_tmp) * list_A_index_only$count[ii]
    }
  }
  
  if (nrow(list_no_index_info) > 0){
    for (ii in 1:nrow(list_no_index_info)){
      
      nC_obs <- list_no_index_info$n_child[ii]
      nA_obs <- list_no_index_info$n_adult[ii]
      
      nTOT<- list_no_index_info$n_reports[ii]
      idx  <- nTOT - nC_obs - nA_obs
      
      p_tmp <- 0
      
      for (j in 0 : idx){
        p_tmp <- p_tmp + 
          fracC * coef_nm(1, nC_obs + j, nA_obs + (idx-j), tables)  + 
          (1 - fracC) * coef_nm(2, nC_obs + j, nA_obs + (idx-j), tables) 
      }
      
      llk_total <- llk_total + log(p_tmp) * list_no_index_info$count[ii]
    }
  }
  
  
  
  llk_total <- llk_total  +
    n_known_child_index * log(fracC) +
    n_known_adult_index * log1p(-fracC)
  
  return(llk_total)
}


loglik_2G_safe <- function(par){
  ll <- tryCatch(loglik_2G(par), error = function(e) NA_real_)
  if(!is.finite(ll)) return(-Inf)
  ll
}

loglike_vec <- function(par){
  if (is.null(dim(par))) {
    return(loglik_2G_safe(par)) 
  } else {
    apply(par, 1, function(row) loglik_2G_safe(row)) 
  }
}

prior <- createPrior(
  density = function(par){
    log_r <- par[1]; log_k <- par[2]; logit_S <- par[3]; logit_f <- par[4]
    dnorm(log_r,  0, 1, log=TRUE) +
      dnorm(log_k,  0, 2, log=TRUE) +
      dnorm(logit_S, 0, 5, log=TRUE) +
      dnorm(logit_f, 0, 5, log=TRUE)
  },
  sampler = function(n = 1){
    cbind(
      rnorm(n, 0, 1),
      rnorm(n, 0, 2),
      rnorm(n, 0, 5),
      rnorm(n, 0, 5)
    )
  },
  lower = rep(-Inf, 4), upper = rep(Inf, 4)
)


setup <- createBayesianSetup(
  likelihood = loglike_vec,
  prior      = prior,
  parallel   = 20 # parallel::detectCores(logical = FALSE) - 2L
) 


# MCMC
set.seed(123)
settings_dezs <- list(iterations = 6e3, nrChains = 4, burnin = 3e3)
out_dezs <- runMCMC(setup, sampler = "DEzs", settings = settings_dezs)

summary(out_dezs)
# > summary(out_dezs)
# # # # # # # # # # # # # # # # # # # # # # # # # # 
# ## MCMC chain summary ## 
# # # # # # # # # # # # # # # # # # # # # # # # # # 
# 
# # MCMC sampler:  DEzs 
# # Nr. Chains:  12 
# # Iterations per chain:  1001 
# # Rejection rate:  0.735 
# # Effective sample size:  840 
# # Runtime:  1453.703  sec. 
# 
# # Parameters
# psf    MAP   2.5% median  97.5%
#   par 1 1.008 -0.272 -0.564 -0.258  0.076
# par 2 1.007 -1.154 -2.260 -1.190  0.186
# par 3 1.033  3.952  2.122  5.140 11.465
# par 4 1.011 -2.706 -4.687 -2.938 -1.788
# 
# ## DIC:  Inf 
# ## Convergence 
# Gelman Rubin multivariate psrf:  1.036 
# 
# ## Correlations 
# par 1  par 2  par 3  par 4
# par 1  1.000 -0.050  0.024  0.024
# par 2 -0.050  1.000 -0.021  0.002
# par 3  0.024 -0.021  1.000 -0.022
# par 4  0.024  0.002 -0.022  1.000

par_names <- c("log_r", "log_k", "logit_S", "logit_fracC")
post_chains <- list()

for (rr in seq_along(out_dezs)) {
  raw_chains <- out_dezs[[rr]]$chain
  
  for (cc in seq_along(raw_chains)) {
    mat <- as.matrix(raw_chains[[cc]])
    
    ## Discard initial state
    mat <- mat[-1L, 1:4, drop = FALSE]
    colnames(mat) <- par_names
    
    nm <- paste0("run", rr, "_internal", cc)
    post_chains[[nm]] <- coda::mcmc(mat, start = 1001L, thin = 1L)
  }
}

post_coda <- coda::mcmc.list(post_chains)
full_mat <- do.call(rbind, lapply(post_chains, as.matrix))
saveRDS(post_coda, "./results/posterior_chains.rds")

## output summaries and diagnostics
print(summary(post_coda))
ess <- coda::effectiveSize(post_coda)
gr <- coda::gelman.diag(post_coda, autoburnin = FALSE)
# print(ess)
# print(gr)
# print(cor(full_mat))
# > print(ess)
# log_r       log_k     logit_S logit_fracC 
# 765.2107    868.9346    759.1528    793.8884 
# > print(gr)
# Potential scale reduction factors:
#   
#   Point est. Upper C.I.
# log_r             1.01       1.01
# log_k             1.01       1.03
# logit_S           1.02       1.03
# logit_fracC       1.01       1.03
# 
# Multivariate psrf
# 
# 1.03
# > print(cor(full_mat))
# log_r        log_k     logit_S  logit_fracC
# log_r        1.00000000 -0.041939946  0.02535417  0.016219369
# log_k       -0.04193995  1.000000000 -0.02015713 -0.003337052
# logit_S      0.02535417 -0.020157133  1.00000000 -0.017266793
# logit_fracC  0.01621937 -0.003337052 -0.01726679  1.000000000

write.csv(data.frame(parameter = par_names,
                     ESS = unname(ess),
                     Rhat = gr$psrf[, 1],
                     Rhat_upper = gr$psrf[, 2]),
          "./results/posterior_diagnostics.csv", row.names = FALSE)

r_full     <- exp(full_mat[, 1])
k_full     <- exp(full_mat[, 2])
S_full     <- plogis(full_mat[, 3])
fracC_full <- plogis(full_mat[, 4])

natural_mat <- cbind(R = r_full, k = k_full, S = S_full, pi = fracC_full)
parameter_summary <- t(apply(natural_mat, 2, quantile,
                             probs = c(.025, .5, .975)))
print(parameter_summary)
# > print(parameter_summary)
# 2.5%        50%     97.5%
# R  0.569919709 0.77266479 1.0790005
# k  0.104310617 0.30419923 1.1988240
# S  0.893813205 0.99418597 0.9999895
# pi 0.009225222 0.05029934 0.1424487
write.csv(parameter_summary, "./results/posterior_parameter_summary.csv")


### plots
## Trace/density plots use the full, cleaned chains.
pdf("./results/SFig1_trace_plot.pdf", width = 8, height = 6)
plot(post_coda)
dev.off()

pdf("./results/SFig3_corrplot_posterior.pdf", width = 8, height = 6)
BayesianTools::correlationPlot(full_mat, thin = 1)
dev.off()

## Preserve the original prior-versus-posterior marginal comparison.
set.seed(124)
prior_mat <- out_dezs[[1]]$setup$prior$sampler(10000)
colnames(prior_mat) <- par_names
pdf("./results/SFig2_marginal_plot.pdf", width = 8, height = 6)
BayesianTools::marginalPlot(full_mat, prior = prior_mat, thin = 1)
dev.off()

## retain 100 entire rows per internal chain.
last_mat <- do.call(rbind, lapply(post_chains, function(x) {
  as.matrix(x)[seq.int(10L, 1000L, by = 10L), , drop = FALSE]
}))
post_draws <- as.data.frame(last_mat)
colnames(post_draws) <- par_names
readr::write_csv(post_draws, "./results/posterior_draws.csv")
