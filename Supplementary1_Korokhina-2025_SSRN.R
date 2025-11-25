### Supplementary Material 1. R-script supporting the workflow

### of the unpublished paper: Korokhina et al. Better Than Chance? Modeling Ceramic Group Distributions in Small Samples (Through the Case-study of Byzantine Amphora From Eastern Europe)

### --- PREPARATORY STAGE ---

## --- Install and load required packages ---

# List of required packages
packages <- c(
  "dplyr", "tidyr", "reshape2",  "tibble", # Data manipulation and transformation
  "ggplot2", "scatterpie", "maps", "ggraph",  "bayesplot", "RColorBrewer", # Data visualisation
  "igraph", # Network analysis
  "vegan", # Multivariate analysis
  "caret", "randomForest", "rstan", "MCMCpack" # Statistical modeling and machine learning
)

# Install missing packages
installed_packages <- rownames(installed.packages())
for(pkg in packages){
  if(!pkg %in% installed_packages)
    install.packages(pkg, dependencies = TRUE)
}

# Load all packages
lapply(packages, library, character.only = TRUE)

# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## --- Load data and select three main amphora groups ---

a <- data.frame(
  row.names = c("Dardanelles", "Chalcis", "S-Crimea"),
  Kyiv = c(36, 15, 3),
  Chernihiv = c(3, 4, 0),
  Pereiaslav = c(7, 13, 0),
  Volodymyr = c(4, 0, 5),
  Czermno = c(15, 0, 1),
  Krakow = c(10, 2, 0)
) # Define amphora group counts
groups_sel <- c("Dardanelles", "Chalcis", "S-Crimea") # select three main groups
a_sel <- a[groups_sel, ]

### --- Compute distances between sites ---

d <- data.frame(
  Site = c("Kyiv", "Chernihiv", "Pereiaslav", "Volodymyr", "Czermno", "Krakow"),
  X = c(50.4501, 51.4982, 50.0656, 50.8499, 50.6017, 50.0647),
  Y = c(30.5234, 31.2893, 31.4454, 24.323, 23.7512, 19.945)
) # Define site coordinates
kyiv_coords <- d[d$Site == "Kyiv", c("X","Y")]
site_coords <- d[match(colnames(a_sel), d$Site), c("X","Y")]
dist_to_kyiv <- sqrt((site_coords$X - kyiv_coords$X)^2 + (site_coords$Y - kyiv_coords$Y)^2)

## --- Core plotting function ---

.plot_amphora_map_core <- function(df, groups_sel, scenario_name, group_colors = NULL) {
  if(is.null(group_colors)) {
    group_colors <- setNames(RColorBrewer::brewer.pal(length(groups_sel), "Set2"), groups_sel)
  }

  world_map <- map_data("world")

  ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "gray90", color = "gray70") +
    scatterpie::geom_scatterpie(data = df, aes(x = Y, y = X, group = Site, r = 1),
                                cols = groups_sel) +
    geom_text(data = df, aes(x = Y, y = X, label = Site), color = "black", size = 3, vjust = -1) +
    coord_fixed(xlim = c(15, 40), ylim = c(45, 55)) +
    theme_minimal() +
    labs(title = scenario_name, x = "Longitude", y = "Latitude") +
    scale_fill_manual(values = group_colors)
}

# ---

### --- STEP 1: EXPLORATORY ANALYSIS OF EMPIRICAL DATA WITH FISHER TEST ---

## --- Plot empirical amphora counts ---

# Plot function
plot_empirical_map <- function(a, groups_sel, d, scenario_name = "Empirical Data", group_colors = NULL) {
  df <- as.data.frame(t(a))
  df$Site <- rownames(df)
  df <- merge(df, d, by = "Site", all.x = TRUE)

  .plot_amphora_map_core(df, groups_sel, scenario_name, group_colors)
}

# Plot empirical data
group_colors <- c("Dardanelles" = "pink", "Chalcis" = "green", "S-Crimea" = "darkred")
plot_empirical_map(a, groups_sel, d, "Empirical Amphora Counts", group_colors)

## --- Pairwise Fisher Test on counts (upper triangle), with convertion to -log10(p) distance-like metrics ---

# Pairwise Fisher test with conversion function (stronger dissimilarities = larger values)
pairwise_fisher_logp <- function(count_data, simulate = FALSE, B = 10000, midp = FALSE) {
  count_mat <- as.matrix(count_data) # convert the data to a numeric matrix and extract site names and their count
  sites <- colnames(count_mat)
  n_sites <- length(sites)

  p_matrix <- matrix(NA, nrow = n_sites, ncol = n_sites,
                     dimnames = list(sites, sites)) # initialize an empty symmetric matrix to store pairwise p-values

  # Loop over all unique site pairs and extract their 2-column contingency subtable for testing
  for (i in 1:(n_sites - 1)) {
    for (j in (i + 1):n_sites) {
      subtable <- count_mat[, c(i, j)]
      # performs Fisher’s exact test for the two sites. If midp = TRUE, applies the mid-p correction. Otherwise, uses the standard p-value
      if (midp) {
        ft <- stats::fisher.test(subtable, simulate.p.value = simulate, B = B)
        obs_prob <- stats::dhyper(subtable[1, 1],
                                  sum(subtable[, 1]),
                                  sum(subtable) - sum(subtable[, 1]),
                                  sum(subtable[1, ]))
        p_val <- ft$p.value - 0.5 * obs_prob
        p_val <- max(min(p_val, 1), 0)
      } else {
        ft <- fisher.test(subtable, simulate.p.value = simulate, B = B)
        p_val <- ft$p.value
      }
      p_matrix[i, j] <- p_val
    }
  }

  # convert p-values to –log10(p) scale, making stronger dissimilarities appear as larger values
  logp_matrix <- -log10(p_matrix)
  logp_matrix[is.infinite(logp_matrix)] <-
    max(logp_matrix[is.finite(logp_matrix)], na.rm = TRUE) + 1
  logp_matrix <- round(logp_matrix, 3)

  return(list(
    p_values_upper = p_matrix,
    logp_upper = logp_matrix
  ))
}

# Execute (get similatiry/dissimilarity distance-like metrics)
results_logp <- pairwise_fisher_logp(a_sel, midp = TRUE)
round(results_logp$logp_upper, 3)

## --- Network graph from -log10(p) distance-like metrics ---

set.seed(100) # set seed

# Network graph function
pairwise_fisher_graph_inverse <- function(count_data, simulate = FALSE, B = 10000, midp = FALSE, alpha_edge = 0.6) {

  res <- pairwise_fisher_logp(count_data, simulate, B, midp)
  logp_matrix <- res$logp_upper
  p_matrix <- res$p_values_upper
  sites <- colnames(count_data)

  # Convert upper triangle to edge list
  edges <- as.data.frame(as.table(logp_matrix))
  colnames(edges) <- c("from", "to", "weight")
  edges <- edges %>%
    dplyr::filter(!is.na(weight) & from != to & as.numeric(from) < as.numeric(to)) %>%
    dplyr::mutate(color = ifelse(p_matrix[cbind(match(from, sites), match(to, sites))] < 0.05,
                                 "Non-significant (p<0.05)", "Significant (p≥0.05)"),
                  color = factor(color, levels = c("Non-significant (p<0.05)", "Significant (p≥0.05)")))

  # Invert edge width: higher weight = thinner edge
  max_w <- max(edges$weight, na.rm = TRUE)
  edges <- edges %>%
    dplyr::mutate(width_inv = max_w - weight + 0.1)

  # Build graph
  g <- igraph::graph_from_data_frame(d = edges, vertices = tibble(name = sites), directed = FALSE)

  # Plot
  plot_graph <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = width_inv, color = color), alpha = alpha_edge) +
    geom_node_point(color = "#F8766D", size = 5) +
    geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    scale_edge_width(range = c(0.5, 5)) +
    scale_edge_color_manual(values = c("Non-significant (p<0.05)" = "grey",
                                               "Significant (p≥0.05)" = "blue")) +
    labs(title = "Site Similarity Network (inverse -log10 Fisher p-values)",
                  edge_color = "Edge significance",
                  edge_width = "Inverse -log10(p)") +
    theme_void() +
    theme(legend.position = "right")

  print(plot_graph)

  return(list(
    p_values_upper = p_matrix,
    logp_upper = logp_matrix,
    graph = g,
    plot = plot_graph
  ))
}

# Plot network graph
results_graph_inv <- pairwise_fisher_graph_inverse(a_sel, midp = TRUE)

# ---

### -- STEP-2: SCENARIO FORMALISATION --

# Define model components
prop_sel <- sweep(a_sel, 2, colSums(a_sel), "/")  # convert counts to proportions
u_list <- lapply(1:ncol(a_sel), function(i) prop_sel[,i]) # local base
mu <- rowMeans(prop_sel)                                   # common pattern
kappa <- prop_sel[, "Kyiv"]                                # Kyiv profile

# Define parameters
S <- ncol(a_sel)       # number of sites
G <- nrow(a_sel)       # number of amphora groups
n_i <- colSums(a_sel)  # total counts per site

alpha <- 0.2           # global pattern weight
beta <- 0.5            # influence weight
phi <- 50              # Dirichlet concentration
lambda <- median(dist_to_kyiv)  # decay for spatial influence

# Define clusters
cluster_A <- 1:floor(S/2)
cluster_B <- (floor(S/2)+1):S

centroid_A <- colMeans(site_coords[cluster_A, ])
centroid_B <- colMeans(site_coords[cluster_B, ])

# ---

### --- STEP-3: SCENARIO SIMULATION ---

# set seed
set.seed(101)

## --- Function for simulating large datasets (x500) per scenario ---

simulate_scenario <- function(type, n_large = 500, phi_val = phi, N_sim_per_site = NULL){
  # N_sim_per_site: integer or vector of length S used for the "internal" total
  # if NULL, use a default large count per site: max( max(n_i), 200 )
  Y_large <- array(0, dim = c(G, S, n_large))
  dimnames(Y_large) <- list(rownames(a_sel), colnames(a_sel), NULL)

  # default internal counts if user didn't supply
  if(is.null(N_sim_per_site)){
    N_default <- max( max(n_i), 200 )   # choose something reasonably large
    N_sim_per_site <- rep(N_default, S)
  } else if(length(N_sim_per_site)==1){
    N_sim_per_site <- rep(as.integer(N_sim_per_site), S)
  }

  for(i in 1:S){
    site_name <- colnames(a_sel)[i]

    # Scenario-specific probability vector
    if(type == "independent"){
      mix <- stats::runif(G); mix <- mix / sum(mix)
    } else if(type == "common"){
      mix <- mu
    } else if(type == "kyiv_influence"){
      w_i <- exp(-dist_to_kyiv[i] / lambda)
      beta_i <- beta * w_i
      mix <- (1 - beta_i) * rep(1/G, G) + beta_i * kappa
    } else if(type == "two_cluster"){
      cluster_pattern_A <- c(0.6, 0.3, 0.1)
      cluster_pattern_B <- c(0.1, 0.4, 0.5)
      if(i %in% cluster_A){
        dist_i <- sqrt((site_coords$X[i] - centroid_A["X"])^2 +
                         (site_coords$Y[i] - centroid_A["Y"])^2)
        beta_i <- beta * exp(-dist_i / lambda)
        mix <- (1 - beta_i) * rep(1/G, G) + beta_i * cluster_pattern_A
      } else {
        dist_i <- sqrt((site_coords$X[i] - centroid_B["X"])^2 +
                         (site_coords$Y[i] - centroid_B["Y"])^2)
        beta_i <- beta * exp(-dist_i / lambda)
        mix <- (1 - beta_i) * rep(1/G, G) + beta_i * cluster_pattern_B
      }
    }

    # Dirichlet-multinomial sampling
    for(n_idx in 1:n_large){
      p_i <- as.numeric(MCMCpack::rdirichlet(1, phi_val * mix))
      Y_large[, i, n_idx] <- as.numeric(stats::rmultinom(1, N_sim_per_site[i], p_i))
    }
  }

  return(Y_large)
}

## --- Function for subsampling (downsample to empirical totals) ---

subsample_from_large <- function(Y_large, target_totals = n_i){
  # pick one draw (index) at random and return G x S matrix
  if(length(dim(Y_large)) != 3) stop("Y_large must be an array with 3 dimensions: G x S x n_large")
  n_large <- dim(Y_large)[3]
  idx_draw <- sample.int(n_large, 1)  # pick one draw
  Y_draw <- Y_large[ , , idx_draw, drop = FALSE]
  Y_draw <- matrix(Y_draw, nrow = dim(Y_draw)[1], ncol = dim(Y_draw)[2],
                   dimnames = list(rownames(Y_large), colnames(Y_large)))

  # Downsample each site's column to match empirical sample size (target_totals)
  Y_sub <- matrix(0, nrow = nrow(Y_draw), ncol = ncol(Y_draw),
                  dimnames = dimnames(Y_draw))
  for(i in 1:ncol(Y_draw)){
    col_counts <- Y_draw[, i]
    tot_sim <- sum(col_counts)
    target_n <- target_totals[i]

    if(tot_sim <= 0){
      # no simulated counts: produce zeros (or could fallback to a small pseudo-count)
      Y_sub[, i] <- rep(0, nrow(Y_draw))
    } else if(tot_sim == target_n){
      Y_sub[, i] <- col_counts
    } else {
      # probabilistic downsample: draw target_n items with probs = col_counts / tot_sim
      probs <- col_counts / tot_sim
      # if target_n is 0 -> result is zeros
      if(target_n > 0){
        Y_sub[, i] <- as.numeric(stats::rmultinom(1, target_n, probs))
      } else {
        Y_sub[, i] <- rep(0, nrow(Y_draw))
      }
    }
  }

  return(Y_sub)
}

# ---

### --- STEP-4: HYPERPARAMETER TUNING FOR SIMULATION MODEL ---

# Set seed
set.seed(102)

# Wrapper allowing simulate_scenario to see local param values
# IMPORTANT: simulate_scenario produces a large array; we subsample here and return a G x S matrix
simulate_with_params <- function(type, beta_val = beta, phi_val = phi, lambda_val = lambda, n_large = 500) {
  # save old globals
  old_beta <- beta; old_phi <- phi; old_lambda <- lambda
  # set globals used inside simulate_scenario
  beta <<- beta_val; phi <<- phi_val; lambda <<- lambda_val

  # simulate large dataset (G x S x n_large)
  Y_large <- simulate_scenario(type, n_large = n_large, phi_val = phi_val)

  # restore globals
  beta <<- old_beta; phi <<- old_phi; lambda <<- old_lambda

  # subsample one draw to produce a G x S counts matrix compatible with sim_to_features()
  Y_sub <- subsample_from_large(Y_large)

  return(Y_sub)
}

# Function for pairwise Fisher tests for simulations, with convertion to -log10(p) distance-like metrics
sim_to_features <- function(Y_sim){
  sites <- colnames(Y_sim)
  feats <- c()

  for(i in 1:(length(sites)-1)){
    for(j in (i+1):length(sites)){
      tab <- matrix(c(Y_sim[,i], Y_sim[,j]), ncol = 2)
      pval <- tryCatch({
        stats::fisher.test(tab, simulate.p.value = FALSE)$p.value
      }, error = function(e) 1)

      d <- -log10(pval)
      feats <- c(feats, d)
      names(feats)[length(feats)] <- paste0(sites[i], "_vs_", sites[j])
    }
  }
  return(feats)
}

# Function to evaluate parameter set with scenario-wise accuracy and confusion matrix
evaluate_params <- function(param, n_per_scenario = 100,
                            scenarios = c("independent","common","kyiv_influence","two_cluster"),
                            classifier_folds = 5, n_large = 500){
  beta_val   <- param$beta
  phi_val    <- param$phi
  lambda_val <- param$lambda

  feats_list <- list()
  labels <- character()
  for(scn in scenarios){
    for(i in 1:n_per_scenario){
      # simulate small-sample dataset by creating large dataset then subsampling
      Ysim <- simulate_with_params(scn, beta_val, phi_val, lambda_val, n_large = n_large)
      feat <- sim_to_features(Ysim)
      feats_list[[length(feats_list)+1]] <- feat
      labels <- c(labels, scn)
    }
  }

  feat_mat <- do.call(rbind, feats_list)
  feat_df <- as.data.frame(feat_mat)
  feat_df$label <- as.factor(labels)

  # Classification with cross-validation and saved predictions
  train_control <- caret::trainControl(method = "cv", number = classifier_folds,
                                       classProbs = FALSE, savePredictions = "final")

  rf_fit <- tryCatch({
    caret::train(label ~ ., data = feat_df, method = "rf",
                 trControl = train_control, ntree = 200, importance = FALSE)
  }, error = function(e){
    rf <- randomForest::randomForest(label ~ ., data = feat_df, ntree = 200)
    list(results = data.frame(Accuracy = NA), finalModel = rf, method="oob",
         pred = data.frame(pred = predict(rf, feat_df), obs = feat_df$label))
  })

  if(inherits(rf_fit, "train")){
    acc <- max(rf_fit$results$Accuracy, na.rm = TRUE)
    method_used <- "cv"
    preds <- rf_fit$pred
  } else {
    rf <- rf_fit$finalModel
    acc <- 1 - rf$err.rate[nrow(rf$err.rate), "OOB"]
    method_used <- "oob"
    preds <- rf_fit$pred
  }

  # Scenario-wise accuracy
  scenario_acc <- sapply(scenarios, function(scn){
    idx <- which(preds$obs == scn)
    mean(preds$pred[idx] == preds$obs[idx], na.rm = TRUE)
  })
  names(scenario_acc) <- scenarios

  # Confusion matrix
  conf_mat <- table(True = preds$obs, Predicted = preds$pred)

  return(list(accuracy = acc,
              method = method_used,
              scenario_accuracy = scenario_acc,
              confusion_matrix = conf_mat,
              n = n_per_scenario * length(scenarios)))
}

## --- Grid search ---

# Set candidate values for parameters
beta_vals   <- seq(0.1, 1.0, by = 0.2)
phi_vals    <- c(10, 30, 50, 100)
lambda_vals <- round(quantile(dist_to_kyiv, probs = c(0.25, 0.5, 0.75)))

# Build parameter grid
param_grid  <- expand.grid(beta = beta_vals, phi = phi_vals, lambda = lambda_vals)

# Storage for results
results_list <- vector("list", nrow(param_grid))

# Number of simulated datasets per scenario per parameter set
n_per_scenario <- 60

# Loop over parameter grid
for(i in 1:nrow(param_grid)){
  p <- as.list(param_grid[i,])
  cat(sprintf("Evaluating grid %d/%d: beta=%.2f phi=%d lambda=%d\n",
              i, nrow(param_grid), p$beta, p$phi, p$lambda))

  res <- evaluate_params(p, n_per_scenario = n_per_scenario,
                         classifier_folds = 5, n_large = 500)

  results_list[[i]] <- list(params = param_grid[i,],
                            accuracy = res$accuracy,
                            scenario_accuracy = res$scenario_accuracy,
                            confusion_matrix = res$confusion_matrix)
}

# Collect results into data.frame
results_df <- do.call(rbind, lapply(results_list, function(x){
  cbind(x$params, accuracy = x$accuracy)
}))
results_df <- as.data.frame(results_df)
results_df$accuracy <- as.numeric(as.character(results_df$accuracy))

# Find best by mean scenario-wise accuracy
scenario_acc_df <- do.call(rbind, lapply(results_list, function(x){
  cbind(x$params, as.data.frame(as.list(x$scenario_accuracy)))
}))
scenario_acc_df[ , 1:3] <- lapply(scenario_acc_df[ , 1:3], as.numeric)
scenario_acc_df[ , 4:7] <- lapply(scenario_acc_df[ , 4:7], as.numeric)
scenario_acc_df$mean_scenario_acc <- rowMeans(scenario_acc_df[, 4:7])

# Pick best row
best_params <- scenario_acc_df[which.max(scenario_acc_df$mean_scenario_acc),
                               c("beta", "phi", "lambda")]

# Assign values to "best" parameters
best_beta   <- best_params$beta
best_phi    <- best_params$phi
best_lambda <- best_params$lambda

cat("Best parameters:\n")
cat("beta =", best_beta, "phi =", best_phi, "lambda =", best_lambda, "\n")

# ---

### --- STEP 5: BAYESIAN MODEL (drop S-Crimea) ---

# Drop S-Crimea row to handle simplex constraint
a_sel_drop <- a_sel[!rownames(a_sel) %in% "S-Crimea", ]

stan_code_drop <- "
data {
  int<lower=1> G;              // number of groups (excluding S-Crimea)
  int<lower=1> S;              // number of sites
  int<lower=0> Y[G,S];         // observed counts (no S-Crimea)
}
parameters {
  simplex[G] theta[S];         // proportions for remaining groups
  real<lower=0> phi;           // concentration parameter
}
model {
  phi ~ exponential(1);
  for (s in 1:S) {
    theta[s] ~ dirichlet(rep_vector(phi, G));
    Y[, s] ~ multinomial(theta[s]);
  }
}
"

stan_model_drop <- rstan::stan_model(model_code = stan_code_drop)

stan_data_drop <- list(
  G = nrow(a_sel_drop),
  S = ncol(a_sel_drop),
  Y = as.matrix(a_sel_drop)
)

set.seed(103)
fit_drop <- rstan::sampling(
  stan_model_drop,
  data = stan_data_drop,
  iter = 2000,
  chains = 4,
  seed = 104,
  control = list(adapt_delta = 0.99)
)

post_drop <- rstan::extract(fit_drop)

## --- Plot phi convergence ---

posterior_drop <- as.array(fit_drop)
mcmc_trace(posterior_drop, pars = "phi") +
  ggtitle("Traceplot of the posterior samples of φ (without S-Crimea)")

posterior_predictive_drop <- function(post_samples, n_i, n_draws = 500) {
  S <- dim(post_samples$theta)[2]
  G <- dim(post_samples$theta)[3]
  draws <- sample(seq_len(dim(post_samples$theta)[1]), n_draws)
  sim_list <- vector("list", n_draws)

  for (d in seq_len(n_draws)) {
    theta_d <- post_samples$theta[draws[d], , ]  # [site × group]
    sim_mat <- matrix(0, G + 1, S)  # add row for reconstructed S-Crimea

    for (s in 1:S) {
      # reconstruct S-Crimea as remaining share
      theta_full <- c(theta_d[s, ], 1 - sum(theta_d[s, ]))
      theta_full[theta_full < 0] <- 0
      theta_full <- theta_full / sum(theta_full)
      sim_mat[, s] <- as.numeric(rmultinom(1, n_i[s], theta_full))
    }

    colnames(sim_mat) <- colnames(a_sel_drop)
    rownames(sim_mat) <- rownames(a_sel)
    sim_list[[d]] <- sim_mat
  }

  return(sim_list)
}

n_i <- colSums(a_sel)
set.seed(105)
pp_sims <- posterior_predictive_drop(post_drop, n_i = n_i, n_draws = 500)

# ---

### --- STEP 6: COMPARE BAYESIAN MODEL AND SCENARIO SIMULATIONS ---

# --- Simulate scenarios using best model parameters ---

set.seed(106)
n_sims <- 1000
scenarios <- c("common", "independent", "kyiv_influence", "two_cluster")

simulations <- lapply(scenarios, function(scn){
  replicate(n_sims,
            simulate_with_params(scn,
                                 beta_val = best_beta,
                                 phi_val = best_phi,
                                 lambda_val = best_lambda,
                                 n_large = 500),
            simplify = FALSE)
})
names(simulations) <- scenarios

## --- Compute -log10(p)-matrices for simulations ---

sim_to_logp_matrix <- function(sim_mat, midp = FALSE) {
  sim_mat <- as.matrix(sim_mat)
  sites <- colnames(sim_mat)
  n_sites <- length(sites)

  logp_matrix <- matrix(NA, nrow = n_sites, ncol = n_sites,
                        dimnames = list(sites, sites))

  for (i in 1:(n_sites - 1)) {
    for (j in (i + 1):n_sites) {
      subtable <- sim_mat[, c(i, j)]
      p_val <- tryCatch({
        if (midp) {
          ft <- fisher.test(subtable)
          obs_prob <- dhyper(subtable[1,1], sum(subtable[,1]), sum(subtable)-sum(subtable[,1]), sum(subtable[1,]))
          pval_adj <- max(min(ft$p.value - 0.5*obs_prob, 1), 0)
        } else {
          pval_adj <- fisher.test(subtable)$p.value
        }
        pval_adj
      }, error = function(e) 1)
      logp_matrix[i, j] <- -log10(p_val + 1e-300)
    }
  }

  # Replace infinities
  if (any(is.infinite(logp_matrix))) {
    logp_matrix[is.infinite(logp_matrix)] <- max(logp_matrix[is.finite(logp_matrix)], na.rm = TRUE) + 1
  }

  return(logp_matrix)
}

sim_logp_list <- lapply(simulations, function(scn_list) {
  lapply(scn_list, sim_to_logp_matrix)
})

## --- Plot the first simulation for each scenario ---

plot_simulation_map <- function(Y_sim, groups_sel, d, scenario_name, group_colors = NULL) {
  df <- as.data.frame(t(Y_sim))
  df$Site <- rownames(df)

  df <- df %>%
    tidyr::pivot_longer(-Site, names_to = "Group", values_to = "Count") %>%
    tidyr::pivot_wider(names_from = Group, values_from = Count)

  df <- merge(df, d, by = "Site", all.x = TRUE)
  df[groups_sel] <- lapply(df[groups_sel], function(x) ifelse(x < 0, 0, x))

  .plot_amphora_map_core(df, groups_sel, scenario_name, group_colors)
}

group_colors <- c("Dardanelles" = "pink", "Chalcis" = "green", "S-Crimea" = "darkred")

for (scn in scenarios) {
  plot_simulation_map(simulations[[scn]][[1]], groups_sel, d,
                      paste("Scenario -", scn), group_colors)
}

## --- Check the accuracy of scenario prediction using Random Forest ---

flatten_logp_matrix <- function(logp_mat) {
  as.vector(logp_mat[upper.tri(logp_mat)])
}

train_data <- do.call(rbind, lapply(names(sim_logp_list), function(scn) {
  scn_list <- sim_logp_list[[scn]]
  do.call(rbind, lapply(scn_list, flatten_logp_matrix))
}))

train_labels <- rep(names(sim_logp_list), each = n_sims)

train_df <- as.data.frame(train_data)
train_df$scenario <- factor(train_labels, levels = scenarios)

set.seed(107)
rf_model <- randomForest::randomForest(scenario ~ ., data = train_df,
                                       importance = TRUE, ntree = 500)

print(rf_model) # print confusion matrix

## --- Classification of posterior predictive draws against empirical data using Random Forest ---

# Flatten posterior predictive draws to match training features
pp_features <- do.call(rbind, lapply(pp_sims, function(sim_mat) {
  flatten_logp_matrix(sim_to_logp_matrix(sim_mat))
}))

pp_df <- as.data.frame(pp_features)
colnames(pp_df) <- colnames(train_df)[1:(ncol(train_df)-1)]  # match names exactly

# Predict probabilities
set.seed(108)
pp_probs <- predict(rf_model, newdata = pp_df, type = "prob")
pp_prob_summary <- colMeans(pp_probs)
print(pp_prob_summary) # print predicted probabilities

## --- Visualization of posterior predictive scenario probabilities ---

pp_probs_long <- reshape2::melt(
  data.frame(draw = 1:nrow(pp_probs), pp_probs),
  id.vars = "draw", variable.name = "scenario", value.name = "probability"
)

ggplot2::ggplot(pp_probs_long, aes(x = scenario, y = probability, fill = scenario)) +
  ggplot2::geom_boxplot(alpha = 0.7) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Posterior Predictive Scenario Probabilities",
    x = "Scenario", y = "Probability"
  )

### --- END ---
