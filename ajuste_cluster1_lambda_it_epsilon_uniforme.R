# ==============================================================================
# ajuste_cluster1_epsilon_Uniforme.R
# Modelo Bayesiano para microrregiões do Cluster 1 (MELHOR reportação)
#
# CONFIRMADO:
#   cluster_ids = 1 → MELHOR reportação (BH, etc.)
#   cluster_ids = 4 → PIOR reportação (Araçuaí, etc.)
#
# epsilon[t] ~ U(0.95, 1) i.i.d.
# lambda[i,t] com FFBS por região
# ==============================================================================

inicio_global <- Sys.time()



pkgs <- c("nimble", "coda", "dplyr", "ggplot2", "tidyr", "readr", "stringr", "tibble")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# ==============================================================================
# 1. DADOS
# ==============================================================================

source("_dataCaseStudy.r")

E_norm    <- E / mean(E)
N_regions <- nrow(Y_mat)
n_times   <- ncol(Y_mat)
p         <- dim(x)[3]

h_mat       <- data$hAI
cluster_ids <- apply(h_mat, 1, sum)  # 1 = melhor, 4 = pior

# Selecionar apenas cluster 1
cluster1_idx <- which(cluster_ids == 1)
N_c1 <- length(cluster1_idx)

cat(sprintf("Cluster 1: %d regiões (BH está aqui? %s)\n", 
            N_c1, any(grepl("BELO_HORIZONTE", rownames(Y_mat)[cluster1_idx]))))

# Filtrar dados
Y_c1 <- Y_mat[cluster1_idx, , drop = FALSE]
E_c1 <- E_norm[cluster1_idx, , drop = FALSE]
x_c1 <- x[cluster1_idx, , , drop = FALSE]

region_names_c1 <- rownames(Y_c1)
anos_label <- colnames(Y_mat)

# ==============================================================================
# 2. SUBGRAFO DE ADJACÊNCIA
# ==============================================================================

old_to_new <- setNames(seq_along(cluster1_idx), cluster1_idx)

adj_c1 <- c()
num_c1 <- integer(N_c1)
weights_c1 <- c()

for (i_new in seq_len(N_c1)) {
  i_old <- cluster1_idx[i_new]
  
  if (data$num[i_old] > 0) {
    start_idx <- if (i_old == 1) 1 else sum(data$num[1:(i_old-1)]) + 1
    end_idx   <- start_idx + data$num[i_old] - 1
    neighbors_old <- data$adj[start_idx:end_idx]
    
    neighbors_new <- old_to_new[as.character(neighbors_old)]
    neighbors_new <- neighbors_new[!is.na(neighbors_new)]
    
    num_c1[i_new] <- length(neighbors_new)
    if (length(neighbors_new) > 0) {
      adj_c1 <- c(adj_c1, neighbors_new)
      weights_c1 <- c(weights_c1, rep(1, length(neighbors_new)))
    }
  }
}

n_adj_c1 <- length(adj_c1)

# ==============================================================================
# 3. REGIÕES PARA MONITORAR
# ==============================================================================

# BH e mais algumas
bh_idx_c1 <- grep("BELO_HORIZONTE_NOVA_LIMA_CAETE", region_names_c1)[1]

set.seed(42)
REGIONS_INTEREST <- c(bh_idx_c1, sample(setdiff(seq_len(N_c1), bh_idx_c1), 11))
REGIONS_INTEREST <- sort(REGIONS_INTEREST)

cat("\nRegiões monitoradas:\n")
print(data.frame(Nome = region_names_c1[REGIONS_INTEREST]))

# ==============================================================================
# 4. CONSTANTES E DADOS
# ==============================================================================

a0 <- 1.0; b0 <- 1.0; w <- 0.85
eps_lower <- 0.95; eps_upper <- 1.00

constants <- list(
  n_regions = N_c1, n_times = n_times, p = p,
  a0 = a0, b0 = b0, w = w,
  eps_lower = eps_lower, eps_upper = eps_upper,
  adj = adj_c1, num = num_c1, weights = weights_c1, n_adj = n_adj_c1
)

data_nimble <- list(Y = Y_c1, E = E_c1, x = x_c1)

# ==============================================================================
# 5. INICIALIZAÇÕES
# ==============================================================================

set.seed(123)
lambda_init <- matrix(1.0, nrow = N_c1, ncol = n_times)

inits_1 <- list(
  lambda  = lambda_init,
  beta    = rep(0, p),
  epsilon = rep(0.975, n_times),
  sigma_s = 0.5,
  s       = rep(0, N_c1)
)
inits_2 <- list(
  lambda  = lambda_init * matrix(runif(N_c1 * n_times, 0.8, 1.2), N_c1, n_times),
  beta    = rnorm(p, 0, 0.3),
  epsilon = runif(n_times, 0.96, 0.99),
  sigma_s = 1.0,
  s       = rep(0, N_c1)
)

# ==============================================================================
# 6. MODELO NIMBLE
# ==============================================================================

code <- nimbleCode({
  
  for (j in 1:p) { beta[j] ~ dnorm(0, sd = 1) }
  
  for (t in 1:n_times) {
    epsilon[t] ~ dunif(eps_lower, eps_upper)
  }
  
  sigma_s ~ T(dt(0, 1, 1), 0, )
  tau_s   <- 1 / (sigma_s^2)
  s[1:n_regions] ~ dcar_normal(
    adj[1:n_adj], weights[1:n_adj], num[1:n_regions], tau_s, zero_mean = 1
  )
  
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      lambda[i, t] ~ dgamma(a0, b0)
    }
  }
  
  for (i in 1:n_regions) {
    for (t in 1:n_times) {
      mu[i, t] <- lambda[i, t] * E[i, t] * epsilon[t] *
        exp(inprod(beta[1:p], x[i, t, 1:p]) + s[i])
      Y[i, t]        ~ dpois(mu[i, t])
      logLik_Y[i, t] <- dpois(Y[i, t], mu[i, t], log = TRUE)
    }
  }
})

# ==============================================================================
# 7. FFBS PARA lambda[i,t]
# ==============================================================================

ffbs_lambda <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    n_times  <- control$n_times; p <- control$p
    a0 <- control$a0; b0 <- control$b0; w <- control$w
    region_i <- control$region_i
    at_buf <- nimNumeric(n_times + 1, 0)
    bt_buf <- nimNumeric(n_times + 1, 0)
    calcNodes   <- model$getDependencies(target, self = FALSE)
    targetNodes <- model$expandNodeNames(target)
    setupOutputs(at_buf, bt_buf)
  },
  run = function() {
    declare(t, integer()); declare(t_idx, integer())
    declare(t_back, integer()); declare(k, integer())
    declare(prod_val, double()); declare(g_it, double()); declare(nu, double())
    
    at_buf[1] <<- a0; bt_buf[1] <<- b0
    
    for (t in 1:n_times) {
      prod_val <- 0.0
      for (k in 1:p) prod_val <- prod_val + model$x[region_i, t, k] * model$beta[k]
      g_it <- model$E[region_i, t] * model$epsilon[t] * 
        exp(prod_val + model$s[region_i])
      at_buf[t + 1] <<- w * at_buf[t] + model$Y[region_i, t]
      bt_buf[t + 1] <<- w * bt_buf[t] + g_it
    }
    
    model$lambda[region_i, n_times] <<- rgamma(1,
                                               shape = at_buf[n_times + 1], rate = max(bt_buf[n_times + 1], 1e-10))
    
    for (t_idx in 1:(n_times - 1)) {
      t_back <- n_times - t_idx
      nu <- rgamma(1, shape = (1 - w) * at_buf[t_back + 1],
                   rate = max(bt_buf[t_back + 1], 1e-10))
      model$lambda[region_i, t_back] <<- nu + w * model$lambda[region_i, t_back + 1]
    }
    
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = targetNodes, logProb = TRUE)
  },
  methods = list(reset = function() {})
)

# ==============================================================================
# 8. COMPILAR E CONFIGURAR
# ==============================================================================

model <- nimbleModel(code = code, constants = constants,
                     data = data_nimble, inits = inits_1, check = FALSE)
Cmodel <- compileNimble(model)

conf <- configureMCMC(model)

# FFBS para cada região
conf$removeSamplers("lambda")
for (i in seq_len(N_c1)) {
  conf$addSampler(
    target = paste0("lambda[", i, ", 1:", n_times, "]"),
    type   = ffbs_lambda,
    control = list(n_times = n_times, p = p, a0 = a0, b0 = b0, w = w, region_i = i)
  )
}

# Slice para epsilon
conf$removeSamplers("epsilon")
for (t in seq_len(n_times)) {
  conf$addSampler(target = paste0("epsilon[", t, "]"), type = "slice")
}

# Monitores
mu_monitors <- unlist(lapply(REGIONS_INTEREST,
                             function(r) paste0("mu[", r, ", ", seq_len(n_times), "]")))

conf$addMonitors(c("beta", "epsilon", "lambda", "logLik_Y", "s", "sigma_s", "tau_s"))
conf$addMonitors(mu_monitors)

Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = model)

# ==============================================================================
# 9. EXECUTAR MCMC
# ==============================================================================

niter <- 50000; nburnin <- 10000; nchains <- 2; thin <- 10

cat(sprintf("niter=%d | nburnin=%d | thin=%d\n", niter, nburnin, thin))

samples <- runMCMC(Cmcmc, niter = niter, nburnin = nburnin, nchains = nchains,
                   thin = thin, inits = list(inits_1, inits_2),
                   samplesAsCodaMCMC = TRUE, summary = FALSE, WAIC = FALSE)

output_dir <- "C:/Users/vlara/OneDrive/Estatistica UFMG/Mestrado/Pesquisa/Aplicação/resultados_cluster1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(samples, file.path(output_dir, "samples.rds"))

# ==============================================================================
# 10. SUMÁRIOS E GRÁFICOS
# ==============================================================================

samples_mat <- as.matrix(samples)
mcmc_list   <- mcmc.list(lapply(1:nchains, function(ch) as.mcmc(samples[[ch]])))

safe_hpd <- function(sv) {
  if (var(sv) < 1e-12) return(c(NA_real_, NA_real_))
  as.numeric(HPDinterval(as.mcmc(sv), prob = 0.95))
}

make_label <- function(names) {
  sapply(names, function(n) {
    n_clean <- gsub("_", " ", n)
    if (nchar(n_clean) > 25) paste0(substr(n_clean, 1, 22), "...") else n_clean
  })
}

# ── epsilon[t] ─────────────────────────────────────────────────────────────

epsilon_names <- paste0("epsilon[", 1:n_times, "]")
eps_summary <- do.call(rbind, lapply(1:n_times, function(t) {
  sv <- samples_mat[, epsilon_names[t]]; hpd <- safe_hpd(sv)
  tibble(Time = t, Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
}))

ggplot(eps_summary, aes(x = Time)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkorange", alpha = 0.3) +
  geom_line(aes(y = Mean), color = "darkorange", linewidth = 0.8) +
  geom_hline(yintercept = c(0.95, 1), linetype = "dashed", color = "grey50") +
  scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Epsilon[t] ~ U(0.95, 1) - Cluster 1", x = "Ano", y = expression(epsilon[t]))
ggsave(file.path(output_dir, "painel_epsilon.png"), width = 10, height = 5)

# ── beta ───────────────────────────────────────────────────────────────────

beta_names <- paste0("beta[", 1:p, "]")
beta_summary <- do.call(rbind, lapply(beta_names, function(nm) {
  sv <- samples_mat[, nm]; hpd <- safe_hpd(sv)
  tibble(Parameter = nm, Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
}))
write_csv(beta_summary, file.path(output_dir, "beta_summary.csv"))

ggplot(beta_summary, aes(x = Parameter, y = Mean)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw(base_size = 12) +
  labs(title = "Coeficientes beta - Cluster 1", y = expression(beta))
ggsave(file.path(output_dir, "beta_posterior.png"), width = 8, height = 5)

# ── lambda[i,t] ────────────────────────────────────────────────────────────

lambda_names <- outer(1:N_c1, 1:n_times, function(i, t) paste0("lambda[", i, ", ", t, "]"))
lambda_sel <- do.call(rbind, lapply(REGIONS_INTEREST, function(i) {
  do.call(rbind, lapply(1:n_times, function(t) {
    sv <- samples_mat[, lambda_names[i, t]]; hpd <- safe_hpd(sv)
    tibble(Region = i, Nome = region_names_c1[i], Time = t,
           Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
  }))
}))
lambda_sel$Label <- make_label(lambda_sel$Nome)

ggplot(lambda_sel, aes(x = Time)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
  geom_line(aes(y = Mean), color = "black", linewidth = 0.8) +
  facet_wrap(~ Label, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "lambda[i,t] - Cluster 1", x = "Ano", y = expression(lambda[it]))
ggsave(file.path(output_dir, "painel_lambda.png"), width = 16, height = 12)

# ── mu[i,t] ────────────────────────────────────────────────────────────────

mu_names <- outer(REGIONS_INTEREST, 1:n_times, function(r, t) paste0("mu[", r, ", ", t, "]"))
mu_sel <- do.call(rbind, lapply(seq_along(REGIONS_INTEREST), function(idx) {
  i <- REGIONS_INTEREST[idx]
  do.call(rbind, lapply(1:n_times, function(t) {
    sv <- samples_mat[, mu_names[idx, t]]; hpd <- safe_hpd(sv)
    tibble(Region = i, Nome = region_names_c1[i], Time = t,
           Mean = mean(sv), Lower = hpd[1], Upper = hpd[2])
  }))
}))
mu_sel$Label <- make_label(mu_sel$Nome)

ggplot(mu_sel, aes(x = Time)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
  geom_line(aes(y = Mean), color = "steelblue", linewidth = 0.8) +
  facet_wrap(~ Label, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = 1:n_times, labels = anos_label) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "mu[i,t] - Cluster 1", x = "Ano", y = expression(mu[it]))
ggsave(file.path(output_dir, "painel_mu.png"), width = 16, height = 12)

# ── Efeito espacial s[i] ───────────────────────────────────────────────────

s_names <- paste0("s[", 1:N_c1, "]")
s_summary <- do.call(rbind, lapply(1:N_c1, function(i) {
  sv <- samples_mat[, s_names[i]]; hpd <- safe_hpd(sv)
  tibble(Region = i, Nome = region_names_c1[i], Mean = mean(sv),
         Lower = hpd[1], Upper = hpd[2])
}))
s_summary$Highlight <- ifelse(grepl("BELO_HORIZONTE_NOVA_LIMA_CAETE", s_summary$Nome), "BH", "Outras")

ggplot(s_summary, aes(x = reorder(Nome, Mean), y = Mean, color = Highlight)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.4, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("BH" = "red", "Outras" = "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "none") +
  labs(title = "Efeito espacial s[i] - Cluster 1", subtitle = "Vermelho: BH",
       y = "s[i]", x = "")
ggsave(file.path(output_dir, "s_posterior.png"), width = 14, height = 6)

# ── WAIC e LPML ───────────────────────────────────────────────────────────

loglik_names <- grep("logLik_Y", colnames(samples_mat), value = TRUE)
lm     <- samples_mat[, loglik_names, drop = FALSE]
lppd   <- sum(apply(lm, 2, function(x) { mx <- max(x); mx + log(mean(exp(x - mx))) }))
p_waic <- sum(apply(lm, 2, var))
waic   <- -2 * (lppd - p_waic)
LPML   <- sum(log(1 / apply(lm, 2, function(x) mean(exp(-x)))))

cat(sprintf("\nWAIC = %.2f | LPML = %.2f\n", waic, LPML))
write_csv(tibble(WAIC = waic, LPML = LPML, lppd = lppd, pWAIC = p_waic),
          file.path(output_dir, "criteria.csv"))

# ── Tempo total ────────────────────────────────────────────────────────────

cat("\nTempo total de execução:\n")
print(Sys.time() - inicio_global)