# ============================================================
# Supplementary R script
# AIC-weighted stochastic character mapping and ASR visualization
# for Jones et al. (2024) and Hand et al. (2023) bat phylogenies
#
# Main analyses:
#   1. Load and clean phylogenies and tip-character data.
#   2. Use the already dated Bayesian Hand trees directly.
#   3. Generate stochastic tip-dated Jones trees with paleotree.
#   4. Fit ER, SYM and ARD transition models for each binary trait.
#   5. Weight stochastic character maps by AIC weights.
#   6. Produce ternary AIC-weight plots, all-tree node-pie PDFs,
#      reference ASR figures, and probability densitrees.
#
# Important:
#   - Hand trees are NOT recalibrated.
#   - Jones uses the original topology and random tip ages sampled
#     inside each taxon's FAD_low--FAD_up range.
#   - Tip values 0.5 are treated as prior probabilities:
#       P(present) = 0.5 and P(absent) = 0.5.
#
# Outputs are written to out_dir.
# ============================================================

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("ggtree", "treeio")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

packages_needed <- c(
  "ape",
  "phytools",
  "phangorn",
  "paleotree",
  "readxl",
  "dplyr",
  "tidyr",
  "Ternary",
  "ggplot2",
  "ggnewscale",
  "ggtree",
  "patchwork",
  "deeptime",
  "grid",
  "scales"
)

for (pkg in packages_needed) {
  install_if_missing(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

set.seed(123)

# ------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------

base_dir <- "C:/Users/qa19267/OneDrive/Documentos/Tesis/DOCENCIA/2025-2026/TFMs/Adriana/Data/R script"

data_file  <- file.path(base_dir, "Data.xlsx")
hand_file  <- file.path(base_dir, "Hand_et_al_2023.tre")
jones_file <- file.path(base_dir, "Jones_et_al_2024.nex")

out_dir <- file.path(base_dir, "ASR_stochastic_mapping_outputs_publication")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Reproducibility
RANDOM_SEED <- 123
set.seed(RANDOM_SEED)

# Choose "test" for a quick dry run or "final" for the full analysis.
# In final mode the intended analysis uses 1000 Jones trees and up to
# 1001 Hand trees, without replacement.
ANALYSIS_MODE <- "test"   # change to "test" for a fast 3-tree check

if (ANALYSIS_MODE == "test") {
  N_JONES_TREES_TO_USE <- 50
  N_HAND_TREES_TO_USE  <- 50
  N_SIMMAP_INTEGRATED_PER_TREE <- 20
} else if (ANALYSIS_MODE == "final") {
  N_JONES_TREES_TO_USE <- 100
  N_HAND_TREES_TO_USE  <- 100
  N_SIMMAP_INTEGRATED_PER_TREE <- 50
} else {
  stop("ANALYSIS_MODE must be either 'test' or 'final'.")
}

traits_to_analyse <- c(
  "Echolocation",
  "Cave_deposit",
  "Cave_karstic_deposit"
)

jones_taxon_col <- "taxa_jones"
hand_taxon_col  <- "taxa_hand"

fad_low_col <- "FAD_low"
fad_up_col  <- "FAD_up"

# Figure size controls for publication outputs
DENSITREE_WIDTH_MM  <- 560
DENSITREE_HEIGHT_MM <- 420
DENSITREE_DPI       <- 600

# ------------------------------------------------------------
# 2. Colours
# ------------------------------------------------------------

cols_echolocation <- c(
  "0" = "#E8D8B8",
  "1" = "#173F5F"
)

cols_cave <- c(
  "0" = "#BFE9F3",
  "1" = "#000000"
)

cols_cave_karstic <- cols_cave

get_trait_cols <- function(trait) {
  if (trait == "Echolocation") return(cols_echolocation)
  if (trait == "Cave_deposit") return(cols_cave)
  if (trait == "Cave_karstic_deposit") return(cols_cave_karstic)
  stop("Unknown trait: ", trait)
}

pal_echo <- c("#F2EFEA", "#7BA6A6", "#243B53")
pal_cave <- c("#BFE8FF", "#5B7C8F", "#000000")

# ------------------------------------------------------------
# 3. General helper functions
# ------------------------------------------------------------

clean_taxon_label <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "Na", "na", "N/A")] <- NA
  x <- trimws(x)
  x <- gsub("'", "", x)
  x <- gsub('"', "", x)
  x <- gsub("\\?", "", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[()]", "", x)
  x <- gsub("__+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

pretty_taxon_name <- function(x) {
  gsub("_", " ", as.character(x))
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

to_num <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
}

as_multiPhylo_safe <- function(trees) {
  if (inherits(trees, "phylo")) {
    trees <- list(trees)
    class(trees) <- "multiPhylo"
  }
  if (!inherits(trees, "multiPhylo")) {
    stop("Object is neither phylo nor multiPhylo.")
  }
  trees
}

read_trees_robust <- function(file) {
  
  trees <- tryCatch(
    ape::read.nexus(file),
    error = function(e) NULL
  )
  
  if (is.null(trees)) {
    trees <- tryCatch(
      ape::read.tree(file),
      error = function(e) {
        stop("Could not read tree file:\n", file, "\nError: ", e$message)
      }
    )
  }
  
  as_multiPhylo_safe(trees)
}

clean_tree_tip_labels <- function(tree) {
  tree$tip.label <- clean_taxon_label(tree$tip.label)
  
  if (anyDuplicated(tree$tip.label)) {
    stop(
      "Duplicated tip labels after cleaning:\n",
      paste(unique(tree$tip.label[duplicated(tree$tip.label)]), collapse = "\n")
    )
  }
  
  tree
}

fix_tree_lengths <- function(tree) {
  
  if (is.null(tree$edge.length)) {
    tree$edge.length <- rep(1, nrow(tree$edge))
  }
  
  tree$edge.length[is.na(tree$edge.length)] <- 1e-8
  tree$edge.length[tree$edge.length <= 0] <- 1e-8
  
  tree
}


fix_tree_for_analysis <- function(tree) {

  tree <- clean_tree_tip_labels(tree)
  tree <- fix_tree_lengths(tree)

  is_bin <- tryCatch(
    ape::is.binary(tree),
    error = function(e) ape::is.binary.phylo(tree)
  )

  if (!isTRUE(is_bin)) {
    tree <- ape::multi2di(tree, random = FALSE)
    tree <- fix_tree_lengths(tree)
  }

  tree
}



write_single_tree_nexus <- function(tree, file) {
  
  if (inherits(tree, "multiPhylo")) {
    tree <- tree[[1]]
  }
  
  if (!inherits(tree, "phylo")) {
    stop("Object to write is not a phylo object: ", file)
  }
  
  one_tree <- list(tree)
  class(one_tree) <- "multiPhylo"
  
  ape::write.nexus(one_tree, file = file)
}

write_multi_tree_nexus <- function(trees, file) {
  trees <- as_multiPhylo_safe(trees)
  ape::write.nexus(trees, file = file)
}

get_mcc_tree <- function(trees) {
  
  trees <- as_multiPhylo_safe(trees)
  
  mcc <- tryCatch(
    phangorn::maxCladeCred(trees, tree = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(mcc)) {
    mcc <- phangorn::maxCladeCred(trees)
  }
  
  if (inherits(mcc, "phylo")) {
    return(fix_tree_for_analysis(mcc))
  }
  
  if (inherits(mcc, "multiPhylo")) {
    return(fix_tree_for_analysis(mcc[[1]]))
  }
  
  if (is.numeric(mcc) && length(mcc) == 1) {
    return(fix_tree_for_analysis(trees[[as.integer(mcc)]]))
  }
  
  stop(
    "phangorn::maxCladeCred() did not return a phylo tree or a usable tree index."
  )
}

sample_tree_set <- function(trees, n, dataset_label) {
  
  trees <- as_multiPhylo_safe(trees)
  
  n_available <- length(trees)
  n_use <- min(n, n_available)
  
  idx <- sample(seq_len(n_available), size = n_use, replace = FALSE)
  
  sampled <- trees[idx]
  class(sampled) <- "multiPhylo"
  
  attr(sampled, "sampled_indices") <- idx
  
  write.csv(
    data.frame(
      dataset = dataset_label,
      sampled_tree_index = idx
    ),
    file.path(out_dir, paste0(dataset_label, "_sampled_tree_indices.csv")),
    row.names = FALSE
  )
  
  sampled
}

get_allowed_taxa <- function(dat, taxon_col, required_value_cols = NULL) {
  
  dd <- dat
  
  dd[[taxon_col]] <- clean_taxon_label(dd[[taxon_col]])
  
  keep <- !is.na(dd[[taxon_col]]) & dd[[taxon_col]] != ""
  
  if (!is.null(required_value_cols)) {
    for (cc in required_value_cols) {
      if (!cc %in% names(dd)) {
        stop("Column not found in Data.xlsx: ", cc)
      }
      keep <- keep & !is.na(to_num(dd[[cc]]))
    }
  }
  
  unique(dd[[taxon_col]][keep])
}

prune_tree_to_allowed_taxa <- function(tree, allowed_taxa, dataset_label, tree_id = NA) {
  
  tree <- clean_tree_tip_labels(tree)
  
  drop_tips <- setdiff(tree$tip.label, allowed_taxa)
  
  if (length(drop_tips) > 0) {
    
    drop_file <- file.path(
      out_dir,
      paste0(dataset_label, "_dropped_tips_not_in_Data_tree_", tree_id, ".csv")
    )
    
    write.csv(
      data.frame(dropped_tip = drop_tips),
      drop_file,
      row.names = FALSE
    )
    
    tree <- ape::drop.tip(tree, drop_tips)
  }
  
  tree <- ape::collapse.singles(tree)
  fix_tree_for_analysis(tree)
}

prune_tree_set_to_allowed_taxa <- function(trees, allowed_taxa, dataset_label) {
  
  trees <- as_multiPhylo_safe(trees)
  
  out <- vector("list", length(trees))
  
  for (i in seq_along(trees)) {
    out[[i]] <- prune_tree_to_allowed_taxa(
      tree = trees[[i]],
      allowed_taxa = allowed_taxa,
      dataset_label = dataset_label,
      tree_id = i
    )
  }
  
  class(out) <- "multiPhylo"
  out
}

get_root_age <- function(tree) {
  max(ape::node.depth.edgelength(tree), na.rm = TRUE)
}

# ------------------------------------------------------------
# 4. Load Data.xlsx
# ------------------------------------------------------------

dat <- readxl::read_excel(data_file)
dat <- as.data.frame(dat)
names(dat) <- trimws(names(dat))

required_cols <- c(
  jones_taxon_col,
  hand_taxon_col,
  fad_low_col,
  fad_up_col,
  traits_to_analyse
)

missing_cols <- setdiff(required_cols, names(dat))

if (length(missing_cols) > 0) {
  stop(
    "These columns are missing from Data.xlsx:\n",
    paste(missing_cols, collapse = ", ")
  )
}

dat[[jones_taxon_col]] <- clean_taxon_label(dat[[jones_taxon_col]])
dat[[hand_taxon_col]]  <- clean_taxon_label(dat[[hand_taxon_col]])

dat[[fad_low_col]] <- to_num(dat[[fad_low_col]])
dat[[fad_up_col]]  <- to_num(dat[[fad_up_col]])

for (trait in traits_to_analyse) {
  dat[[trait]] <- to_num(dat[[trait]])
}

# ------------------------------------------------------------
# 5. Load and prepare Hand trees
#    Hand is already dated. Do not recalibrate.
# ------------------------------------------------------------

hand_allowed_taxa <- get_allowed_taxa(
  dat = dat,
  taxon_col = hand_taxon_col,
  required_value_cols = traits_to_analyse
)

hand_trees_all <- read_trees_robust(hand_file)
hand_trees_all <- prune_tree_set_to_allowed_taxa(
  trees = hand_trees_all,
  allowed_taxa = hand_allowed_taxa,
  dataset_label = "Hand"
)

hand_trees_all <- lapply(hand_trees_all, fix_tree_for_analysis)
class(hand_trees_all) <- "multiPhylo"

# ------------------------------------------------------------
# 6. Load Jones original topology and date tips with paleotree
# ------------------------------------------------------------

jones_allowed_taxa <- get_allowed_taxa(
  dat = dat,
  taxon_col = jones_taxon_col,
  required_value_cols = c(fad_low_col, fad_up_col, traits_to_analyse)
)

jones_source_trees <- read_trees_robust(jones_file)
jones_source_tree <- jones_source_trees[[1]]

jones_source_tree <- prune_tree_to_allowed_taxa(
  tree = jones_source_tree,
  allowed_taxa = jones_allowed_taxa,
  dataset_label = "Jones_original",
  tree_id = 1
)

get_tip_ages_for_jones <- function(dat, taxon_col, low_col, up_col, mode = c("random", "mean")) {
  
  mode <- match.arg(mode)
  
  dd <- dat
  dd[[taxon_col]] <- clean_taxon_label(dd[[taxon_col]])
  
  dd <- dd[!is.na(dd[[taxon_col]]) & dd[[taxon_col]] != "", ]
  
  low <- to_num(dd[[low_col]])
  up  <- to_num(dd[[up_col]])
  
  younger <- pmin(low, up, na.rm = TRUE)
  older   <- pmax(low, up, na.rm = TRUE)
  
  if (mode == "random") {
    ages <- stats::runif(length(younger), min = younger, max = older)
  } else {
    ages <- (younger + older) / 2
  }
  
  names(ages) <- dd[[taxon_col]]
  ages
}

date_one_jones_tree_paleotree <- function(tree, tip_ages) {
  
  tree <- fix_tree_for_analysis(tree)
  
  ages <- tip_ages[tree$tip.label]
  
  if (any(!is.finite(ages))) {
    missing <- tree$tip.label[!is.finite(ages)]
    stop(
      "Missing FAD ages for Jones tips:\n",
      paste(missing, collapse = "\n")
    )
  }
  
  time_data <- cbind(
    first = ages,
    last = ages
  )
  
  rownames(time_data) <- names(ages)
  
  dated <- paleotree::timePaleoPhy(
    tree = tree,
    timeData = time_data,
    type = "mbl",
    vartime = 1,
    ntrees = 1,
    dateTreatment = "firstLast"
  )
  
  if (inherits(dated, "multiPhylo")) {
    dated <- dated[[1]]
  }
  
  fix_tree_for_analysis(dated)
}

generate_jones_dated_trees <- function(tree, dat, n_trees) {
  
  out <- vector("list", n_trees)
  age_records <- vector("list", n_trees)
  
  for (i in seq_len(n_trees)) {
    
    message("Dating Jones tree ", i, " / ", n_trees)
    
    ages_i <- get_tip_ages_for_jones(
      dat = dat,
      taxon_col = jones_taxon_col,
      low_col = fad_low_col,
      up_col = fad_up_col,
      mode = "random"
    )
    
    out[[i]] <- date_one_jones_tree_paleotree(
      tree = tree,
      tip_ages = ages_i
    )
    
    age_records[[i]] <- data.frame(
      tree_number = i,
      taxon = names(ages_i),
      sampled_tip_age_Ma = as.numeric(ages_i),
      stringsAsFactors = FALSE
    )
  }
  
  class(out) <- "multiPhylo"
  
  write.csv(
    dplyr::bind_rows(age_records),
    file.path(out_dir, "Jones_random_tip_ages_used_for_paleotree.csv"),
    row.names = FALSE
  )
  
  out
}

jones_trees_all <- generate_jones_dated_trees(
  tree = jones_source_tree,
  dat = dat,
  n_trees = N_JONES_TREES_TO_USE
)

write_multi_tree_nexus(
  jones_trees_all,
  file.path(out_dir, "Jones_randomly_dated_trees_paleotree_TEST_OR_FINAL.nex")
)

saveRDS(
  jones_trees_all,
  file.path(out_dir, "Jones_randomly_dated_trees_paleotree_TEST_OR_FINAL.rds")
)

# Jones summary tree:
# original topology dated using mean age inside each FAD interval.
jones_mean_ages <- get_tip_ages_for_jones(
  dat = dat,
  taxon_col = jones_taxon_col,
  low_col = fad_low_col,
  up_col = fad_up_col,
  mode = "mean"
)

jones_mean_tree <- date_one_jones_tree_paleotree(
  tree = jones_source_tree,
  tip_ages = jones_mean_ages
)

write_single_tree_nexus(
  jones_mean_tree,
  file.path(out_dir, "Jones_original_topology_mean_tip_age_tree.nex")
)

# ------------------------------------------------------------
# 7. Select Hand trees and get Hand MCC tree
# ------------------------------------------------------------

hand_trees <- sample_tree_set(
  trees = hand_trees_all,
  n = N_HAND_TREES_TO_USE,
  dataset_label = "Hand"
)

jones_trees <- sample_tree_set(
  trees = jones_trees_all,
  n = N_JONES_TREES_TO_USE,
  dataset_label = "Jones"
)

hand_mcc_tree <- get_mcc_tree(hand_trees)

write_single_tree_nexus(
  hand_mcc_tree,
  file.path(out_dir, "Hand_MCC_tree.nex")
)

saveRDS(hand_mcc_tree, file.path(out_dir, "Hand_MCC_tree.rds"))
saveRDS(jones_mean_tree, file.path(out_dir, "Jones_original_topology_mean_tip_age_tree.rds"))

# ------------------------------------------------------------
# 8. Trait prior matrices
# ------------------------------------------------------------

make_tip_prior_matrix <- function(tree, dat, taxon_col, trait_col) {
  
  if (!taxon_col %in% names(dat)) {
    stop("Column not found in Data.xlsx: ", taxon_col)
  }
  
  if (!trait_col %in% names(dat)) {
    stop("Column not found in Data.xlsx: ", trait_col)
  }
  
  dd <- dat
  dd[[taxon_col]] <- clean_taxon_label(dd[[taxon_col]])
  dd[[trait_col]] <- to_num(dd[[trait_col]])
  
  dd <- dd[!is.na(dd[[taxon_col]]) & dd[[taxon_col]] != "", ]
  dd <- dd[!is.na(dd[[trait_col]]), ]
  
  p1 <- dd[[trait_col]]
  names(p1) <- dd[[taxon_col]]
  
  missing_taxa <- setdiff(tree$tip.label, names(p1))
  
  if (length(missing_taxa) > 0) {
    
    missing_file <- file.path(
      out_dir,
      paste0("missing_taxa_", taxon_col, "_", trait_col, ".csv")
    )
    
    write.csv(
      data.frame(
        missing_taxon = missing_taxa,
        dataset_taxon_column = taxon_col,
        trait = trait_col
      ),
      missing_file,
      row.names = FALSE
    )
    
    stop(
      "Some tree taxa are missing from Data.xlsx for trait ",
      trait_col, ".\nMissing taxa saved in:\n",
      missing_file
    )
  }
  
  p1 <- p1[tree$tip.label]
  
  if (any(!p1 %in% c(0, 0.5, 1))) {
    stop(
      "Trait ", trait_col,
      " contains values other than 0, 0.5, 1 for selected tree taxa."
    )
  }
  
  prior <- cbind(
    "0" = 1 - p1,
    "1" = p1
  )
  
  rownames(prior) <- tree$tip.label
  
  prior
}

# ------------------------------------------------------------
# 9. AIC-weighted stochastic mapping
# ------------------------------------------------------------

aic_fun <- function(logL, k) {
  2 * k - 2 * logL
}

aic_weights_fun <- function(aic_values) {
  delta <- aic_values - min(aic_values, na.rm = TRUE)
  w <- exp(-0.5 * delta)
  w / sum(w, na.rm = TRUE)
}

n_rate_parameters <- function(model, n_states = 2) {
  
  if (model == "ER") {
    return(1)
  }
  
  if (model == "SYM") {
    return(n_states * (n_states - 1) / 2)
  }
  
  if (model == "ARD") {
    return(n_states * (n_states - 1))
  }
  
  stop("Unknown model: ", model)
}

allocate_nsim_by_weights <- function(weights, total_nsim) {
  
  weights <- weights / sum(weights)
  
  raw <- total_nsim * weights
  nsim <- floor(raw)
  
  remainder <- total_nsim - sum(nsim)
  
  if (remainder > 0) {
    add_to <- order(raw - nsim, decreasing = TRUE)[seq_len(remainder)]
    nsim[add_to] <- nsim[add_to] + 1
  }
  
  names(nsim) <- names(weights)
  nsim
}

as_multi_simmap <- function(x) {
  
  if (inherits(x, "simmap")) {
    xx <- list(x)
    class(xx) <- c("multiSimmap", "multiPhylo")
    return(xx)
  }
  
  if (inherits(x, "multiPhylo")) {
    class(x) <- c("multiSimmap", "multiPhylo")
  }
  
  x
}

extract_logLik_from_make_simmap_fit <- function(fit) {
  
  if (inherits(fit, "multiPhylo") && length(fit) >= 1) {
    fit1 <- fit[[1]]
  } else {
    fit1 <- fit
  }
  
  ll <- fit1$logL
  
  if (is.null(ll)) {
    ll <- attr(fit1, "logL")
  }
  
  if (is.null(ll) || !is.finite(as.numeric(ll))) {
    stop("Could not extract logLik from make.simmap result.")
  }
  
  as.numeric(ll)
}

fit_models_and_make_integrated_simmaps <- function(tree,
                                                   tip_prior,
                                                   total_nsim = 20,
                                                   root_prior = "estimated") {
  
  tree <- fix_tree_for_analysis(tree)
  
  models <- c("ER", "SYM", "ARD")
  
  logL <- sapply(
    models,
    function(m) {
      fit <- phytools::make.simmap(
        tree = tree,
        x = tip_prior,
        model = m,
        nsim = 1,
        pi = root_prior,
        message = FALSE
      )
      extract_logLik_from_make_simmap_fit(fit)
    }
  )
  
  k <- sapply(models, n_rate_parameters, n_states = 2)
  
  aic <- mapply(aic_fun, logL, k)
  names(aic) <- models
  
  weights <- aic_weights_fun(aic)
  names(weights) <- models
  
  nsim_by_model <- allocate_nsim_by_weights(
    weights = weights,
    total_nsim = total_nsim
  )
  
  simmap_list <- list()
  
  for (m in models) {
    
    if (nsim_by_model[m] > 0) {
      
      sm <- phytools::make.simmap(
        tree = tree,
        x = tip_prior,
        model = m,
        nsim = nsim_by_model[m],
        pi = root_prior,
        message = FALSE
      )
      
      sm <- as_multi_simmap(sm)
      simmap_list[[m]] <- sm
    }
  }
  
  integrated <- unlist(simmap_list, recursive = FALSE)
  class(integrated) <- c("multiSimmap", "multiPhylo")
  
  aic_table <- data.frame(
    model = models,
    logLik = as.numeric(logL[models]),
    k = as.numeric(k[models]),
    AIC = as.numeric(aic[models]),
    AIC_weight = as.numeric(weights[models]),
    nsim = as.numeric(nsim_by_model[models]),
    stringsAsFactors = FALSE
  )
  
  list(
    simmaps = integrated,
    aic_table = aic_table
  )
}

run_dataset_asr <- function(dataset_label,
                            trees,
                            dat,
                            taxon_col,
                            traits,
                            total_nsim_per_tree) {
  
  message("\nRunning ASR for ", dataset_label, "...")
  
  dataset_dir <- file.path(out_dir, dataset_label)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- list()
  
  for (trait in traits) {
    
    message("  Trait: ", trait)
    
    trait_dir <- file.path(dataset_dir, trait)
    dir.create(trait_dir, recursive = TRUE, showWarnings = FALSE)
    
    maps_by_tree <- vector("list", length(trees))
    aic_all <- list()
    
    for (i in seq_along(trees)) {
      
      message("    Tree ", i, " / ", length(trees))
      
      tree_i <- fix_tree_for_analysis(trees[[i]])
      
      tip_prior <- make_tip_prior_matrix(
        tree = tree_i,
        dat = dat,
        taxon_col = taxon_col,
        trait_col = trait
      )
      
      fit_i <- fit_models_and_make_integrated_simmaps(
        tree = tree_i,
        tip_prior = tip_prior,
        total_nsim = total_nsim_per_tree,
        root_prior = "estimated"
      )
      
      maps_by_tree[[i]] <- fit_i$simmaps
      
      aic_i <- fit_i$aic_table
      aic_i$tree_number <- i
      aic_i$dataset <- dataset_label
      aic_i$trait <- trait
      
      aic_all[[i]] <- aic_i
    }
    
    aic_all_df <- dplyr::bind_rows(aic_all)
    
    write.csv(
      aic_all_df,
      file.path(trait_dir, paste0(dataset_label, "_", trait, "_AIC_by_tree_long.csv")),
      row.names = FALSE
    )
    
    aic_wide <- aic_all_df %>%
      dplyr::select(tree_number, model, AIC_weight) %>%
      tidyr::pivot_wider(
        names_from = model,
        values_from = AIC_weight,
        names_prefix = "w"
      )
    
    write.csv(
      aic_wide,
      file.path(trait_dir, paste0(dataset_label, "_", trait, "_AIC_weights_by_tree_wide.csv")),
      row.names = FALSE
    )
    
    saveRDS(
      maps_by_tree,
      file.path(trait_dir, paste0(dataset_label, "_", trait, "_integrated_simmaps_by_tree.rds"))
    )
    
    results[[trait]] <- list(
      maps_by_tree = maps_by_tree,
      aic_long = aic_all_df,
      aic_wide = aic_wide
    )
  }
  
  saveRDS(
    results,
    file.path(dataset_dir, paste0(dataset_label, "_ASR_integrated_results.rds"))
  )
  
  results
}

hand_asr <- run_dataset_asr(
  dataset_label = "Hand",
  trees = hand_trees,
  dat = dat,
  taxon_col = hand_taxon_col,
  traits = traits_to_analyse,
  total_nsim_per_tree = N_SIMMAP_INTEGRATED_PER_TREE
)

jones_asr <- run_dataset_asr(
  dataset_label = "Jones",
  trees = jones_trees,
  dat = dat,
  taxon_col = jones_taxon_col,
  traits = traits_to_analyse,
  total_nsim_per_tree = N_SIMMAP_INTEGRATED_PER_TREE
)

# ------------------------------------------------------------
# 10. Ternary plots: one point = one tree
# ------------------------------------------------------------
# 
# plot_ternary_trait <- function(aic_wide,
#                                main_title,
#                                output_file) {
#   
#   required_cols <- c("wER", "wSYM", "wARD")
#   missing_cols <- setdiff(required_cols, names(aic_wide))
#   
#   if (length(missing_cols) > 0) {
#     stop("Missing AIC-weight columns: ", paste(missing_cols, collapse = ", "))
#   }
#   
#   pdf(output_file, width = 6, height = 6)
#   
#   Ternary::TernaryPlot(
#     alab = "SYM",
#     blab = "ARD",
#     clab = "ER",
#     axis.labels = seq(0, 1, by = 0.1),
#     main = main_title
#   )
#   
#   # Same order as your reference script:
#   # ER, SYM, ARD.
#   Ternary::TernaryPoints(
#     aic_wide[, c("wER", "wSYM", "wARD")],
#     cex = 0.85,
#     col = "black",
#     pch = 16
#   )
#   
#   mean_point <- matrix(
#     colMeans(aic_wide[, c("wER", "wSYM", "wARD")], na.rm = TRUE),
#     nrow = 1
#   )
#   
#   Ternary::TernaryPoints(
#     mean_point,
#     cex = 1.5,
#     col = "red",
#     pch = 4,
#     lwd = 2
#   )
#   
#   dev.off()
# }
# 
# plot_ternary_dataset <- function(dataset_label, asr_results) {
#   
#   dataset_dir <- file.path(out_dir, dataset_label)
#   
#   for (trait in names(asr_results)) {
#     
#     trait_clean <- gsub("_", " ", trait)
#     
#     plot_ternary_trait(
#       aic_wide = asr_results[[trait]]$aic_wide,
#       main_title = paste0(dataset_label, " - ", trait_clean),
#       output_file = file.path(
#         dataset_dir,
#         paste0(dataset_label, "_", trait, "_ternary_AIC_weights.pdf")
#       )
#     )
#   }
#   
#   pdf(
#     file.path(dataset_dir, paste0(dataset_label, "_ALL_TRAITS_ternary_AIC_weights.pdf")),
#     width = 18,
#     height = 6
#   )
#   
#   par(mfrow = c(1, 3), mar = c(1, 1, 4, 1))
#   
#   for (trait in names(asr_results)) {
#     
#     aic_wide <- asr_results[[trait]]$aic_wide
#     
#     Ternary::TernaryPlot(
#       alab = "SYM",
#       blab = "ARD",
#       clab = "ER",
#       axis.labels = seq(0, 1, by = 0.1),
#       main = gsub("_", " ", trait)
#     )
#     
#     Ternary::TernaryPoints(
#       aic_wide[, c("wER", "wSYM", "wARD")],
#       cex = 0.75,
#       col = "black",
#       pch = 16
#     )
#     
#     mean_point <- matrix(
#       colMeans(aic_wide[, c("wER", "wSYM", "wARD")], na.rm = TRUE),
#       nrow = 1
#     )
#     
#     Ternary::TernaryPoints(
#       mean_point,
#       cex = 1.4,
#       col = "red",
#       pch = 4,
#       lwd = 2
#     )
#   }
#   
#   dev.off()
# }
# 
# plot_ternary_dataset("Hand", hand_asr)
# plot_ternary_dataset("Jones", jones_asr)



# ------------------------------------------------------------
# 10. Ternary plots: one point = one tree
# ------------------------------------------------------------

plot_ternary_trait <- function(aic_wide,
                               main_title,
                               output_file) {
  
  required_cols <- c("wER", "wSYM", "wARD")
  missing_cols <- setdiff(required_cols, names(aic_wide))
  
  if (length(missing_cols) > 0) {
    stop("Missing AIC-weight columns: ", paste(missing_cols, collapse = ", "))
  }
  
  pdf(output_file, width = 7, height = 7)
  
  Ternary::TernaryPlot(
    alab = "SYM",
    blab = "ER",         # <--- CORREGIDO: ER a la derecha
    clab = "ARD",        # <--- CORREGIDO: ARD abajo
    axis.labels = seq(0, 1, by = 0.1),
    main = main_title,
    cex.main = 1.8,
    lab.cex = 1.5,
    axis.cex = 1.1
  )
  
  # CORREGIDO: Los datos ahora se pasan en el mismo orden que las etiquetas (SYM, ER, ARD)
  Ternary::TernaryPoints(
    aic_wide[, c("wSYM", "wER", "wARD")], 
    cex = 0.85,
    col = "black",
    pch = 16
  )
  
  mean_point <- matrix(
    colMeans(aic_wide[, c("wSYM", "wER", "wARD")], na.rm = TRUE),
    nrow = 1
  )
  
  Ternary::TernaryPoints(
    mean_point,
    cex = 1.8,
    col = "red",
    pch = 4,
    lwd = 3
  )
  
  dev.off()
}

plot_ternary_dataset <- function(dataset_label, asr_results) {
  
  dataset_dir <- file.path(out_dir, dataset_label)
  
  for (trait in names(asr_results)) {
    
    trait_clean <- gsub("_", " ", trait)
    
    plot_ternary_trait(
      aic_wide = asr_results[[trait]]$aic_wide,
      main_title = paste0(dataset_label, " - ", trait_clean),
      output_file = file.path(
        dataset_dir,
        paste0(dataset_label, "_", trait, "_ternary_AIC_weights.pdf")
      )
    )
  }
  
  pdf(
    file.path(dataset_dir, paste0(dataset_label, "_ALL_TRAITS_ternary_AIC_weights.pdf")),
    width = 18,
    height = 6.5
  )
  
  par(mfrow = c(1, 3), mar = c(2, 2, 4, 2), oma = c(0, 0, 4, 0))
  
  for (trait in names(asr_results)) {
    
    aic_wide <- asr_results[[trait]]$aic_wide
    
    Ternary::TernaryPlot(
      alab = "SYM",
      blab = "ER",         # <--- CORREGIDO
      clab = "ARD",        # <--- CORREGIDO
      axis.labels = seq(0, 1, by = 0.1),
      main = gsub("_", " ", trait),
      cex.main = 2.5,
      lab.cex = 2.0,
      axis.cex = 1.3
    )
    
    # CORREGIDO: Orden sincronizado
    Ternary::TernaryPoints(
      aic_wide[, c("wSYM", "wER", "wARD")],
      cex = 0.75,
      col = "black",
      pch = 16
    )
    
    mean_point <- matrix(
      colMeans(aic_wide[, c("wSYM", "wER", "wARD")], na.rm = TRUE),
      nrow = 1
    )
    
    Ternary::TernaryPoints(
      mean_point,
      cex = 2.0,
      col = "red",
      pch = 4,
      lwd = 3
    )
  }
  
  mtext(paste0("Model AIC Weights - ", dataset_label, " Dataset"), outer = TRUE, cex = 2.5, font = 2)
  
  dev.off()
}

plot_ternary_dataset("Hand", hand_asr)
plot_ternary_dataset("Jones", jones_asr)
# ------------------------------------------------------------
# 11. Chronostratigraphic chart for base R plots
# ------------------------------------------------------------
#
add_geo_scale_base <- function(tree,
                               y_offset_tips = 1.0,   # Ajustado: más pegado al árbol
                               height_tips = 0.8,     # Ajustado: cajas más finas y proporcionadas
                               cex_labels = 0.9) {    # Letra de las épocas
  
  root_age <- get_root_age(tree)
  usr <- par("usr")
  
  y_top_epoch <- usr[3] - y_offset_tips
  h <- height_tips
  
  y_bottom_epoch  <- y_top_epoch - h
  y_bottom_period <- y_bottom_epoch - h
  
  epoch_df <- data.frame(
    name = c(
      "Pleist.", "Pliocene", "Miocene", "Oligocene",
      "Eocene", "Paleocene", "Late Cret.", "Early Cret.", "Jurassic"
    ),
    young = c(0, 2.58, 5.333, 23.03, 33.9, 56.0, 66.0, 100.5, 145.0),
    old   = c(2.58, 5.333, 23.03, 33.9, 56.0, 66.0, 100.5, 145.0, 201.4),
    col   = c(
      "#FFF7BC", "#FFFFB2", "#FFFF00", "#FDD49E",
      "#FDAE6B", "#FDBF6F", "#B3CDE3", "#6497B1", "#CAB2D6"
    )
  )
  
  period_df <- data.frame(
    name = c("Quat.", "Neogene", "Paleogene", "Cretaceous", "Jurassic"),
    young = c(0, 2.58, 23.03, 66.0, 145.0),
    old   = c(2.58, 23.03, 66.0, 145.0, 201.4),
    col   = c("#FFF7BC", "#FFFF00", "#FD8D3C", "#80B1D3", "#BC80BD")
  )
  
  par(xpd = NA)
  
  for (i in seq_len(nrow(epoch_df))) {
    x_left  <- root_age - min(epoch_df$old[i], root_age)
    x_right <- root_age - epoch_df$young[i]
    
    if (x_right >= 0 && x_left <= usr[2]) {
      rect(
        xleft = max(x_left, 0), ybottom = y_bottom_epoch,
        xright = min(x_right, usr[2]), ytop = y_top_epoch,
        col = epoch_df$col[i], border = "black", lwd = 0.45
      )
      text(
        x = (max(x_left, 0) + min(x_right, usr[2])) / 2,
        y = (y_bottom_epoch + y_top_epoch) / 2,
        labels = epoch_df$name[i], srt = 90, cex = cex_labels
      )
    }
  }
  
  for (i in seq_len(nrow(period_df))) {
    x_left  <- root_age - min(period_df$old[i], root_age)
    x_right <- root_age - period_df$young[i]
    
    if (x_right >= 0 && x_left <= usr[2]) {
      rect(
        xleft = max(x_left, 0), ybottom = y_bottom_period,
        xright = min(x_right, usr[2]), ytop = y_bottom_epoch,
        col = period_df$col[i], border = "black", lwd = 0.7
      )
      text(
        x = (max(x_left, 0) + min(x_right, usr[2])) / 2,
        y = (y_bottom_period + y_bottom_epoch) / 2,
        labels = period_df$name[i], cex = cex_labels + 0.05
      )
    }
  }
  
  axis_ages <- seq(0, ceiling(root_age / 10) * 10, by = 10)
  axis_at <- root_age - axis_ages
  
  axis(
    side = 1, at = axis_at, labels = axis_ages,
    pos = y_bottom_period, las = 2, 
    cex.axis = 0.8  # Hacemos los números del eje más grandes
  )
  
  text(
    x = root_age / 2, y = y_bottom_period - (h * 1.5),
    labels = "Age Ma", 
    cex = 1.3       # Letra de "Age Ma" más grande
  )
  
  par(xpd = FALSE)
}

plot_tree_with_chronoscale <- function(tree, main_title, fsize = NULL) {
  
  tree <- fix_tree_for_analysis(tree)
  n_tips <- length(tree$tip.label)
  
  if (is.null(fsize)) {
    fsize <- max(0.22, min(0.75, 18 / n_tips))
  }
  
  root_age <- get_root_age(tree)
  
  # Aumentamos el margen de arriba (de 4 a 8) para el título grande
  par(mar = c(12, 1, 8, 18)) 
  
  ape::plot.phylo(
    tree,
    direction = "rightwards", show.tip.label = TRUE,
    cex = fsize, label.offset = root_age * 0.012,
    no.margin = FALSE, x.lim = c(0, root_age * 1.65),
    y.lim = c(1, n_tips + (n_tips * 0.06)),
    main = main_title,
    cex.main = 2.5,   # NUEVO: Título principal más grande
    font.main = 2     # NUEVO: Título en negrita
  )
  
  add_geo_scale_base(tree)
}



# add_geo_scale_base <- function(tree,
#                                y_offset_tips = 2.5,   # NUEVO: Distancia absoluta (equivale a 2.5 especies)
#                                height_tips = 2.0,     # NUEVO: Altura de las cajas (equivale a 2 especies)
#                                cex_labels = 0.85) {
#   
#   root_age <- get_root_age(tree)
#   usr <- par("usr")
#   
#   # AHORA EL TAMAÑO ES ABSOLUTO Y NO DEPENDE DEL TAMAÑO DEL ÁRBOL
#   y_top_epoch <- usr[3] - y_offset_tips
#   h <- height_tips
#   
#   y_bottom_epoch  <- y_top_epoch - h
#   y_bottom_period <- y_bottom_epoch - h
#   
#   epoch_df <- data.frame(
#     name = c(
#       "Pleist.", "Pliocene", "Miocene", "Oligocene",
#       "Eocene", "Paleocene", "Late Cret.", "Early Cret.", "Jurassic"
#     ),
#     young = c(0, 2.58, 5.333, 23.03, 33.9, 56.0, 66.0, 100.5, 145.0),
#     old   = c(2.58, 5.333, 23.03, 33.9, 56.0, 66.0, 100.5, 145.0, 201.4),
#     col   = c(
#       "#FFF7BC", "#FFFFB2", "#FFFF00", "#FDD49E",
#       "#FDAE6B", "#FDBF6F", "#B3CDE3", "#6497B1", "#CAB2D6"
#     )
#   )
#   
#   period_df <- data.frame(
#     name = c("Quat.", "Neogene", "Paleogene", "Cretaceous", "Jurassic"),
#     young = c(0, 2.58, 23.03, 66.0, 145.0),
#     old   = c(2.58, 23.03, 66.0, 145.0, 201.4),
#     col   = c("#FFF7BC", "#FFFF00", "#FD8D3C", "#80B1D3", "#BC80BD")
#   )
#   
#   par(xpd = NA)
#   
#   for (i in seq_len(nrow(epoch_df))) {
#     
#     x_left  <- root_age - min(epoch_df$old[i], root_age)
#     x_right <- root_age - epoch_df$young[i]
#     
#     if (x_right >= 0 && x_left <= usr[2]) {
#       rect(
#         xleft = max(x_left, 0),
#         ybottom = y_bottom_epoch,
#         xright = min(x_right, usr[2]),
#         ytop = y_top_epoch,
#         col = epoch_df$col[i],
#         border = "black",
#         lwd = 0.45
#       )
#       
#       text(
#         x = (max(x_left, 0) + min(x_right, usr[2])) / 2,
#         y = (y_bottom_epoch + y_top_epoch) / 2,
#         labels = epoch_df$name[i],
#         srt = 90,
#         cex = cex_labels
#       )
#     }
#   }
#   
#   for (i in seq_len(nrow(period_df))) {
#     
#     x_left  <- root_age - min(period_df$old[i], root_age)
#     x_right <- root_age - period_df$young[i]
#     
#     if (x_right >= 0 && x_left <= usr[2]) {
#       rect(
#         xleft = max(x_left, 0),
#         ybottom = y_bottom_period,
#         xright = min(x_right, usr[2]),
#         ytop = y_bottom_epoch,
#         col = period_df$col[i],
#         border = "black",
#         lwd = 0.7
#       )
#       
#       text(
#         x = (max(x_left, 0) + min(x_right, usr[2])) / 2,
#         y = (y_bottom_period + y_bottom_epoch) / 2,
#         labels = period_df$name[i],
#         cex = cex_labels + 0.05
#       )
#     }
#   }
#   
#   axis_ages <- seq(0, ceiling(root_age / 10) * 10, by = 10)
#   axis_at <- root_age - axis_ages
#   
#   # EJE TEMPORAL ANCLADO EXACTAMENTE DEBAJO DE LAS CAJAS
#   axis(
#     side = 1,
#     at = axis_at,
#     labels = axis_ages,
#     pos = y_bottom_period,  # <--- NUEVO: Se engancha al borde de los rectángulos
#     las = 2,
#     cex.axis = 0.65
#   )
#   
#   # TEXTO "Age Ma" 
#   text(
#     x = root_age / 2, 
#     y = y_bottom_period - (h * 1.5), # <--- NUEVO: Baja el texto proporcionalmente
#     labels = "Age Ma",
#     cex = 0.85
#   )
#   
#   par(xpd = FALSE)
# }
# 
# plot_tree_with_chronoscale <- function(tree,
#                                        main_title,
#                                        fsize = NULL) {
#   
#   tree <- fix_tree_for_analysis(tree)
#   n_tips <- length(tree$tip.label)
#   
#   if (is.null(fsize)) {
#     fsize <- max(0.22, min(0.75, 18 / n_tips))
#   }
#   
#   root_age <- get_root_age(tree)
#   
#   par(mar = c(16, 1, 4, 18)) # <--- CAMBIO: Aumentamos el margen inferior de 10 a 16
#   
#   ape::plot.phylo(
#     tree,
#     direction = "rightwards",
#     show.tip.label = TRUE,
#     cex = fsize,
#     label.offset = root_age * 0.012,
#     no.margin = FALSE,
#     x.lim = c(0, root_age * 1.65),
#     y.lim = c(1, n_tips + (n_tips * 0.06)), # <--- NUEVO: Añade un 6% de espacio en blanco arriba
#     main = main_title
#   )
#   
#   # HEMOS BORRADO EL BLOQUE ape::axisPhylo(...) QUE ESTABA AQUÍ PARA EVITAR DUPLICADOS
#   
#   add_geo_scale_base(tree)
# }

# ------------------------------------------------------------
# 12. Node pies: Echolocation + Cave trait
#     One PDF with one tree per page
# ------------------------------------------------------------

# ============================================================
# PATCH: robust node pies from simmap objects
# Fixes:
#   Error: Could not align ACE matrix rows with internal node numbers
# Also removes deprecated ape::is.binary.tree() warning.
# fix_tree_for_analysis() is defined in Section 3 and is repeated here
# harmlessly for compatibility if this section is sourced independently.
# ============================================================

fix_tree_for_analysis <- function(tree) {
  
  tree <- clean_tree_tip_labels(tree)
  tree <- fix_tree_lengths(tree)
  
  is_bin <- tryCatch(
    ape::is.binary(tree),
    error = function(e) ape::is.binary.phylo(tree)
  )
  
  if (!isTRUE(is_bin)) {
    tree <- ape::multi2di(tree, random = FALSE)
    tree <- fix_tree_lengths(tree)
  }
  
  tree
}

strip_simmap_to_phylo <- function(smap) {
  
  tr <- smap
  
  tr$maps <- NULL
  tr$mapped.edge <- NULL
  
  class(tr) <- "phylo"
  
  tr <- fix_tree_for_analysis(tr)
  
  tr
}

get_tree_from_simmaps <- function(simmaps) {
  
  if (inherits(simmaps, "simmap")) {
    return(strip_simmap_to_phylo(simmaps))
  }
  
  if (inherits(simmaps, "multiSimmap") || inherits(simmaps, "multiPhylo") || is.list(simmaps)) {
    return(strip_simmap_to_phylo(simmaps[[1]]))
  }
  
  stop("The object supplied is not a simmap or multiSimmap object.")
}

get_node_states_from_one_simmap <- function(smap) {
  
  ntip <- ape::Ntip(smap)
  nnode <- ape::Nnode(smap)
  internal_nodes <- (ntip + 1):(ntip + nnode)
  
  all_states <- rep(NA_character_, ntip + nnode)
  names(all_states) <- as.character(seq_len(ntip + nnode))
  
  root_node <- setdiff(smap$edge[, 1], smap$edge[, 2])
  root_node <- root_node[1]
  
  for (ii in seq_len(nrow(smap$edge))) {
    
    parent <- smap$edge[ii, 1]
    child  <- smap$edge[ii, 2]
    
    branch_map <- smap$maps[[ii]]
    branch_states <- names(branch_map)
    
    if (length(branch_states) == 0) next
    
    start_state <- branch_states[1]
    end_state   <- branch_states[length(branch_states)]
    
    # State at the child node = final state along the incoming branch
    all_states[as.character(child)] <- end_state
    
    # State at the root = initial state along one of its outgoing branches
    if (parent == root_node && is.na(all_states[as.character(parent)])) {
      all_states[as.character(parent)] <- start_state
    }
  }
  
  all_states[as.character(internal_nodes)]
}

get_ace_internal_nodes <- function(simmaps) {
  
  if (inherits(simmaps, "simmap")) {
    simmaps <- list(simmaps)
    class(simmaps) <- c("multiSimmap", "multiPhylo")
  }
  
  base_tree <- get_tree_from_simmaps(simmaps)
  
  ntip <- ape::Ntip(base_tree)
  nnode <- ape::Nnode(base_tree)
  internal_nodes <- (ntip + 1):(ntip + nnode)
  
  state_mat <- do.call(
    cbind,
    lapply(simmaps, get_node_states_from_one_simmap)
  )
  
  if (is.null(dim(state_mat))) {
    state_mat <- matrix(state_mat, ncol = 1)
  }
  
  rownames(state_mat) <- as.character(internal_nodes)
  
  p0 <- rowMeans(state_mat == "0", na.rm = TRUE)
  p1 <- rowMeans(state_mat == "1", na.rm = TRUE)
  
  # If something is still unresolved at a node, set it to 0.5 / 0.5
  bad <- !is.finite(p0) | !is.finite(p1) | ((p0 + p1) == 0)
  
  p0[bad] <- 0.5
  p1[bad] <- 0.5
  
  ace <- cbind(
    "0" = p0,
    "1" = p1
  )
  
  rownames(ace) <- as.character(internal_nodes)
  
  ace
}

make_dual_trait_pie_matrix <- function(simmaps_trait_1,
                                       simmaps_trait_2) {
  
  tree_1 <- get_tree_from_simmaps(simmaps_trait_1)
  tree_2 <- get_tree_from_simmaps(simmaps_trait_2)
  
  if (!setequal(tree_1$tip.label, tree_2$tip.label)) {
    stop("The two simmap objects do not contain the same tip labels.")
  }
  
  ace1 <- get_ace_internal_nodes(simmaps_trait_1)
  ace2 <- get_ace_internal_nodes(simmaps_trait_2)
  
  common_nodes <- intersect(rownames(ace1), rownames(ace2))
  
  ace1 <- ace1[common_nodes, , drop = FALSE]
  ace2 <- ace2[common_nodes, , drop = FALSE]
  
  pie_mat <- cbind(
    Trait1_absent  = 0.5 * ace1[, "0"],
    Trait1_present = 0.5 * ace1[, "1"],
    Trait2_absent  = 0.5 * ace2[, "0"],
    Trait2_present = 0.5 * ace2[, "1"]
  )
  
  list(
    pie = pie_mat,
    nodes = as.integer(common_nodes),
    tree = tree_1
  )
}

plot_all_tree_pages_with_dual_pies <- function(dataset_label,
                                               trees,
                                               asr_results,
                                               trait_1,
                                               trait_2,
                                               output_pdf) {
  
  # Use the tree stored inside the simmap object, not the original tree.
  # This prevents node-number mismatches if multi2di() changed the topology.
  first_plot_tree <- get_tree_from_simmaps(asr_results[[trait_1]]$maps_by_tree[[1]])
  n_tips <- length(first_plot_tree$tip.label)
  
  pdf_height <- max(12, min(34, n_tips * 0.32))
  
  pdf(output_pdf, width = 26, height = pdf_height, onefile = TRUE)
  
  cols_1 <- get_trait_cols(trait_1)
  cols_2 <- get_trait_cols(trait_2)
  
  add_dual_pie_legend_page(dataset_label, trait_1, trait_2, cols_1, cols_2)
  
  for (i in seq_along(asr_results[[trait_1]]$maps_by_tree)) {
    
    pie_obj <- make_dual_trait_pie_matrix(
      simmaps_trait_1 = asr_results[[trait_1]]$maps_by_tree[[i]],
      simmaps_trait_2 = asr_results[[trait_2]]$maps_by_tree[[i]]
    )
    
    tree_i <- pie_obj$tree
    
    pie_cols <- c(
      cols_1["0"],
      cols_1["1"],
      cols_2["0"],
      cols_2["1"]
    )
    
    plot_tree_with_chronoscale(
      tree_i,
      main_title = paste0(
        dataset_label,
        " tree ",
        i,
        " - node pies: ",
        trait_1,
        " + ",
        trait_2
      )
    )
    
    ape::nodelabels(
      node = pie_obj$nodes,
      pie = pie_obj$pie,
      piecol = pie_cols,
      cex = 0.42,
      frame = "none"
    )
  }
  
  dev.off()
}

add_dual_pie_legend_page <- function(dataset_label,
                                     trait_1,
                                     trait_2,
                                     cols_1,
                                     cols_2) {
  
  plot.new()
  
  title(
    main = paste0(dataset_label, ": ", trait_1, " + ", trait_2),
    cex.main = 1.3
  )
  
  legend(
    "center",
    legend = c(
      paste0(trait_1, " absent"),
      paste0(trait_1, " present"),
      paste0(trait_2, " absent"),
      paste0(trait_2, " present")
    ),
    fill = c(cols_1["0"], cols_1["1"], cols_2["0"], cols_2["1"]),
    bty = "n",
    cex = 1.15,
    ncol = 2
  )
}

plot_all_tree_pages_with_dual_pies(
  dataset_label = "Hand",
  trees = hand_trees,
  asr_results = hand_asr,
  trait_1 = "Echolocation",
  trait_2 = "Cave_deposit",
  output_pdf = file.path(out_dir, "Hand_Echolocation_vs_CaveDeposit_node_pies_ALL_TREES.pdf")
)

plot_all_tree_pages_with_dual_pies(
  dataset_label = "Hand",
  trees = hand_trees,
  asr_results = hand_asr,
  trait_1 = "Echolocation",
  trait_2 = "Cave_karstic_deposit",
  output_pdf = file.path(out_dir, "Hand_Echolocation_vs_CaveKarsticDeposit_node_pies_ALL_TREES.pdf")
)

plot_all_tree_pages_with_dual_pies(
  dataset_label = "Jones",
  trees = jones_trees,
  asr_results = jones_asr,
  trait_1 = "Echolocation",
  trait_2 = "Cave_deposit",
  output_pdf = file.path(out_dir, "Jones_Echolocation_vs_CaveDeposit_node_pies_ALL_TREES.pdf")
)

plot_all_tree_pages_with_dual_pies(
  dataset_label = "Jones",
  trees = jones_trees,
  asr_results = jones_asr,
  trait_1 = "Echolocation",
  trait_2 = "Cave_karstic_deposit",
  output_pdf = file.path(out_dir, "Jones_Echolocation_vs_CaveKarsticDeposit_node_pies_ALL_TREES.pdf")
)

# ------------------------------------------------------------
# 13. Reference tree plots
#     Plain reference phylogenies and ASR node pies on:
#       - Hand maximum clade credibility tree
#       - Jones original topology with mean tip ages
# ------------------------------------------------------------

reference_dir <- file.path(out_dir, "Reference_ASR_trees")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)

plot_summary_tree_plain <- function(tree,
                                    dataset_label,
                                    output_pdf) {

  tree <- fix_tree_for_analysis(tree)

  pdf(
    output_pdf,
    width = 24,
    height = max(12, min(34, length(tree$tip.label) * 0.32))
  )

  plot_tree_with_chronoscale(
    tree,
    main_title = paste0(dataset_label, " reference tree")
  )

  ape::nodelabels(cex = 0.35, frame = "circle", bg = "white")

  dev.off()
}

descendant_tip_signature <- function(tree, node) {

  tree <- fix_tree_for_analysis(tree)
  tips <- phangorn::Descendants(tree, node, type = "tips")[[1]]
  paste(sort(tree$tip.label[tips]), collapse = "||")
}

get_clade_p1_from_simmaps <- function(simmaps) {

  tree <- get_tree_from_simmaps(simmaps)
  ace  <- get_ace_internal_nodes(simmaps)

  nodes <- as.integer(rownames(ace))
  sigs <- vapply(nodes, function(nd) {
    descendant_tip_signature(tree, nd)
  }, character(1))

  p1 <- ace[, "1"]
  names(p1) <- sigs

  p1
}

aggregate_trait_p1_to_reference <- function(reference_tree,
                                            asr_results,
                                            trait_name) {

  reference_tree <- fix_tree_for_analysis(reference_tree)

  ref_nodes <- (ape::Ntip(reference_tree) + 1):(ape::Ntip(reference_tree) + ape::Nnode(reference_tree))
  ref_sigs <- vapply(ref_nodes, function(nd) {
    descendant_tip_signature(reference_tree, nd)
  }, character(1))

  p1_by_tree <- lapply(
    asr_results[[trait_name]]$maps_by_tree,
    get_clade_p1_from_simmaps
  )

  out <- lapply(seq_along(ref_nodes), function(i) {

    values_i <- vapply(p1_by_tree, function(v) {
      if (ref_sigs[i] %in% names(v)) {
        as.numeric(v[ref_sigs[i]])
      } else {
        NA_real_
      }
    }, numeric(1))

    data.frame(
      node = ref_nodes[i],
      clade_signature = ref_sigs[i],
      p1_mean = mean(values_i, na.rm = TRUE),
      p1_sd = stats::sd(values_i, na.rm = TRUE),
      n_matching_trees = sum(is.finite(values_i)),
      n_total_trees = length(values_i),
      stringsAsFactors = FALSE
    )
  })

  out <- dplyr::bind_rows(out)

  out$p1_mean[!is.finite(out$p1_mean)] <- NA_real_
  out$p1_sd[!is.finite(out$p1_sd)] <- NA_real_

  out
}

make_reference_dual_trait_pies <- function(reference_tree,
                                           asr_results,
                                           trait_1,
                                           trait_2,
                                           output_csv) {

  p1_trait_1 <- aggregate_trait_p1_to_reference(
    reference_tree = reference_tree,
    asr_results = asr_results,
    trait_name = trait_1
  )

  p1_trait_2 <- aggregate_trait_p1_to_reference(
    reference_tree = reference_tree,
    asr_results = asr_results,
    trait_name = trait_2
  )

  names(p1_trait_1)[names(p1_trait_1) %in% c("p1_mean", "p1_sd", "n_matching_trees")] <-
    c("trait_1_p1_mean", "trait_1_p1_sd", "trait_1_n_matching_trees")

  names(p1_trait_2)[names(p1_trait_2) %in% c("p1_mean", "p1_sd", "n_matching_trees")] <-
    c("trait_2_p1_mean", "trait_2_p1_sd", "trait_2_n_matching_trees")

  merged <- dplyr::left_join(
    p1_trait_1,
    p1_trait_2[, c(
      "node",
      "trait_2_p1_mean",
      "trait_2_p1_sd",
      "trait_2_n_matching_trees"
    )],
    by = "node"
  )

  merged$trait_1 <- trait_1
  merged$trait_2 <- trait_2

  write.csv(merged, output_csv, row.names = FALSE)

  keep <- is.finite(merged$trait_1_p1_mean) &
    is.finite(merged$trait_2_p1_mean) &
    merged$trait_1_n_matching_trees > 0 &
    merged$trait_2_n_matching_trees > 0

  merged_keep <- merged[keep, , drop = FALSE]

  pie_mat <- cbind(
    Trait1_absent  = 0.5 * (1 - merged_keep$trait_1_p1_mean),
    Trait1_present = 0.5 * merged_keep$trait_1_p1_mean,
    Trait2_absent  = 0.5 * (1 - merged_keep$trait_2_p1_mean),
    Trait2_present = 0.5 * merged_keep$trait_2_p1_mean
  )

  rownames(pie_mat) <- as.character(merged_keep$node)

  list(
    nodes = merged_keep$node,
    pie = pie_mat,
    table = merged
  )
}

plot_reference_tree_with_dual_asr_pies <- function(reference_tree,
                                                   dataset_label,
                                                   asr_results,
                                                   trait_1,
                                                   trait_2,
                                                   output_prefix,
                                                   fsize_tips = 2.0,       # NUEVO
                                                   espacio_tips = 0.8) {   # NUEVO
  
  reference_tree <- fix_tree_for_analysis(reference_tree)
  
  cols_1 <- get_trait_cols(trait_1)
  cols_2 <- get_trait_cols(trait_2)
  
  pie_cols <- c(
    cols_1["0"],
    cols_1["1"],
    cols_2["0"],
    cols_2["1"]
  )
  
  pie_obj <- make_reference_dual_trait_pies(
    reference_tree = reference_tree,
    asr_results = asr_results,
    trait_1 = trait_1,
    trait_2 = trait_2,
    output_csv = paste0(output_prefix, "_node_probabilities.csv")
  )
  
  # BLOQUE PDF
  pdf(
    paste0(output_prefix, ".pdf"),
    width = 24,
    height = max(12, length(reference_tree$tip.label) * espacio_tips) # USAMOS LA VARIABLE
  )
  
  plot_tree_with_chronoscale(
    reference_tree,
    main_title = paste0(
      dataset_label,
      " reference ASR: ",
      trait_1,
      " + ",
      trait_2
    ),
    fsize = fsize_tips   # USAMOS LA VARIABLE
  )
  
  ape::nodelabels(
    node = pie_obj$nodes,
    pie = pie_obj$pie,
    piecol = pie_cols,
    cex = 0.65,          # Aumentamos un poco el tamaño de los "pies" para que acompañen a la letra
    frame = "none"
  )
  
  legend(
    "topleft",
    inset = c(0, 0),
    xpd = TRUE,
    legend = c(
      paste0(trait_1, " absent"),
      paste0(trait_1, " present"),
      paste0(trait_2, " absent"),
      paste0(trait_2, " present")
    ),
    fill = pie_cols,
    bty = "n",
    cex = 1.6,
    ncol = 1
  )
  
  dev.off()
  
  # BLOQUE PNG
  png(
    paste0(output_prefix, ".png"),
    width = 7200,
    height = max(3600, length(reference_tree$tip.label) * (espacio_tips * 300)), # USAMOS LA VARIABLE
    res = 300
  )
  
  # (Deja el resto de la función png igual, pero recuerda poner fsize = fsize_tips en el plot_tree_with_chronoscale de aquí también y el cex=0.70 en los nodelabels)

  # ... (código anterior: pdf(...) o png(...) )
  
  plot_tree_with_chronoscale(
    reference_tree,
    main_title = paste0(
      dataset_label,
      " reference ASR: ",
      trait_1,
      " + ",
      trait_2
    ),
    fsize = 1.5   # <--- NUEVO: Aumenta el tamaño de la letra de las especies (el original era ~0.75)
  )
  
  ape::nodelabels(
    node = pie_obj$nodes,
    pie = pie_obj$pie,
    piecol = pie_cols,
    cex = 0.65,
    frame = "none"
  )
  
  legend(
    "topleft",
    inset = c(0, -0.08),  # <--- NUEVO: Desplaza la leyenda hacia el margen superior en blanco
    xpd = TRUE,           # <--- NUEVO: Permite que la leyenda se dibuje fuera del límite estricto del gráfico
    legend = c(
      paste0(trait_1, " absent"),
      paste0(trait_1, " present"),
      paste0(trait_2, " absent"),
      paste0(trait_2, " present")
    ),
    fill = pie_cols,
    bty = "n",
    cex = 1.6,            # <--- MODIFICADO: He subido un poco el tamaño de la leyenda (antes 0.85) para que sea más legible
    ncol = 1
  )
  
  dev.off()

  invisible(pie_obj)
}

plot_summary_tree_plain(
  hand_mcc_tree,
  "Hand MCC",
  file.path(reference_dir, "Hand_MCC_tree_chronostratigraphic.pdf")
)

plot_summary_tree_plain(
  jones_mean_tree,
  "Jones original topology with mean tip ages",
  file.path(reference_dir, "Jones_original_topology_mean_tip_age_tree_chronostratigraphic.pdf")
)

hand_reference_cave <- plot_reference_tree_with_dual_asr_pies(
  reference_tree = hand_mcc_tree,
  dataset_label = "Hand MCC",
  asr_results = hand_asr,
  trait_1 = "Echolocation",
  trait_2 = "Cave_deposit",
  output_prefix = file.path(reference_dir, "Hand_MCC_reference_ASR_Echolocation_vs_Cave_deposit"),
  fsize_tips = 1.5,      # Letra gigante para Hand
  espacio_tips = 1.0     # Mucho espacio vertical
)

hand_reference_karst <- plot_reference_tree_with_dual_asr_pies(
  reference_tree = hand_mcc_tree,
  dataset_label = "Hand MCC",
  asr_results = hand_asr,
  trait_1 = "Echolocation",
  trait_2 = "Cave_karstic_deposit",
  output_prefix = file.path(reference_dir, "Hand_MCC_reference_ASR_Echolocation_vs_Cave_karstic_deposit"),
  fsize_tips = 3.5,
  espacio_tips = 1.0
)

jones_reference_cave <- plot_reference_tree_with_dual_asr_pies(
  reference_tree = jones_mean_tree,
  dataset_label = "Jones mean-age tree",
  asr_results = jones_asr,
  trait_1 = "Echolocation",
  trait_2 = "Cave_deposit",
  output_prefix = file.path(reference_dir, "Jones_mean_age_reference_ASR_Echolocation_vs_Cave_deposit"),
  fsize_tips = 1.5,      # Letra muy grande para Jones
  espacio_tips = 0.8     # Documento muy alto para que no se pisen las ramas
)

jones_reference_karst <- plot_reference_tree_with_dual_asr_pies(
  reference_tree = jones_mean_tree,
  dataset_label = "Jones mean-age tree",
  asr_results = jones_asr,
  trait_1 = "Echolocation",
  trait_2 = "Cave_karstic_deposit",
  output_prefix = file.path(reference_dir, "Jones_mean_age_reference_ASR_Echolocation_vs_Cave_karstic_deposit"),
  fsize_tips = 2.2,
  espacio_tips = 0.8
)

# ------------------------------------------------------------
# 14. FAST PROBABILITY DENSITREES, OSTRACODERM-STYLE
#     Fast approach:
#       1. Extract P(state = 1) for tips and internal nodes.
#       2. Join those probabilities to each tree object.
#       3. Fortify trees.
#       4. Plot with ggtree::ggdensitree().
#
#     This replaces the heavier manual-segment densitree.
# ------------------------------------------------------------

if (!exists("out_dir")) {
  out_dir <- file.path(
    getwd(),
    "outputs",
    paste0("ASR_stochastic_mapping_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  )
}

if (!exists("figures_dir")) {
  figures_dir <- file.path(out_dir, "figures")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

densitree_dir <- file.path(figures_dir, "08_probability_densitrees")
dir.create(densitree_dir, recursive = TRUE, showWarnings = FALSE)


DENSITREE_WIDTH_MM  <- 520
DENSITREE_HEIGHT_MM <- 360
DENSITREE_DPI       <- 300

# If central labels appear vertically inverted relative to the tips,
# change this to TRUE and re-run only the densitree section.
REVERSE_CENTRAL_LABEL_ORDER <- FALSE

# ------------------------------------------------------------
# 14.1 Helpers: extract tree and node probabilities from simmaps
# ------------------------------------------------------------

strip_simmap_to_phylo <- function(smap) {
  tr <- smap
  tr$maps <- NULL
  tr$mapped.edge <- NULL
  class(tr) <- "phylo"
  fix_tree_for_analysis(tr)
}

get_tree_from_simmaps <- function(simmaps) {
  if (inherits(simmaps, "simmap")) {
    return(strip_simmap_to_phylo(simmaps))
  }
  if (inherits(simmaps, "multiSimmap") || inherits(simmaps, "multiPhylo") || is.list(simmaps)) {
    return(strip_simmap_to_phylo(simmaps[[1]]))
  }
  stop("Object is not a simmap or multiSimmap.")
}

get_node_states_from_one_simmap <- function(smap) {
  
  ntip <- ape::Ntip(smap)
  nnode <- ape::Nnode(smap)
  internal_nodes <- (ntip + 1):(ntip + nnode)
  
  all_states <- rep(NA_character_, ntip + nnode)
  names(all_states) <- as.character(seq_len(ntip + nnode))
  
  root_node <- setdiff(smap$edge[, 1], smap$edge[, 2])[1]
  
  for (ii in seq_len(nrow(smap$edge))) {
    
    parent <- smap$edge[ii, 1]
    child  <- smap$edge[ii, 2]
    
    branch_map <- smap$maps[[ii]]
    branch_states <- names(branch_map)
    
    if (length(branch_states) == 0) next
    
    start_state <- branch_states[1]
    end_state   <- branch_states[length(branch_states)]
    
    all_states[as.character(child)] <- end_state
    
    if (parent == root_node && is.na(all_states[as.character(parent)])) {
      all_states[as.character(parent)] <- start_state
    }
  }
  
  all_states[as.character(internal_nodes)]
}

get_internal_node_probabilities_from_simmaps <- function(simmaps) {
  
  if (inherits(simmaps, "simmap")) {
    simmaps <- list(simmaps)
    class(simmaps) <- c("multiSimmap", "multiPhylo")
  }
  
  base_tree <- get_tree_from_simmaps(simmaps)
  
  ntip <- ape::Ntip(base_tree)
  nnode <- ape::Nnode(base_tree)
  internal_nodes <- (ntip + 1):(ntip + nnode)
  
  state_mat <- do.call(
    cbind,
    lapply(simmaps, get_node_states_from_one_simmap)
  )
  
  if (is.null(dim(state_mat))) {
    state_mat <- matrix(state_mat, ncol = 1)
  }
  
  rownames(state_mat) <- as.character(internal_nodes)
  
  p1 <- rowMeans(state_mat == "1", na.rm = TRUE)
  p0 <- rowMeans(state_mat == "0", na.rm = TRUE)
  
  bad <- !is.finite(p1) | !is.finite(p0) | ((p0 + p1) == 0)
  p1[bad] <- 0.5
  
  data.frame(
    node = internal_nodes,
    state = as.numeric(p1),
    stringsAsFactors = FALSE
  )
}

get_tip_probabilities_for_trait <- function(tree, dat, taxon_col, trait_col) {
  
  prior <- make_tip_prior_matrix(
    tree = tree,
    dat = dat,
    taxon_col = taxon_col,
    trait_col = trait_col
  )
  
  data.frame(
    node = seq_len(ape::Ntip(tree)),
    state = as.numeric(prior[tree$tip.label, "1"]),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 14.2 Build fortified trees for ggdensitree
# ------------------------------------------------------------

build_fast_probability_fortified_trees <- function(trees,
                                                   asr_results,
                                                   dat,
                                                   taxon_col,
                                                   trait_name,
                                                   tip_order) {
  
  fortified_trees <- vector("list", length(trees))
  
  for (i in seq_along(trees)) {
    
    message("Preparing densitree data: ", trait_name, " | tree ", i, " / ", length(trees))
    
    simmaps_i <- asr_results[[trait_name]]$maps_by_tree[[i]]
    tree_i <- get_tree_from_simmaps(simmaps_i)
    
    tip_df <- get_tip_probabilities_for_trait(
      tree = tree_i,
      dat = dat,
      taxon_col = taxon_col,
      trait_col = trait_name
    )
    
    node_df <- get_internal_node_probabilities_from_simmaps(simmaps_i)
    
    state_df <- dplyr::bind_rows(tip_df, node_df)
    state_df$node <- as.numeric(state_df$node)
    
    tree_with_states <- dplyr::full_join(tree_i, state_df, by = "node")
    
    fort_i <- ggplot2::fortify(tree_with_states)
    fort_i$tree <- i
    
    # Convert root-to-tip distance into geological age.
    # In dated trees with extant or youngest tips near 0 Ma:
    #   x = max depth - distance from root.
    fort_i$x <- max(fort_i$x, na.rm = TRUE) - fort_i$x
    
    # ggtree uses 'branch' internally for branch segments.
    # It must be transformed consistently with x.
    if ("branch" %in% names(fort_i)) {
      fort_i$branch <- max(fort_i$branch, na.rm = TRUE) - fort_i$branch
    }
    
    fortified_trees[[i]] <- fort_i
  }
  
  tips_present <- fortified_trees[[1]]$label[fortified_trees[[1]]$isTip]
  tip_order_for_plot <- tip_order[tip_order %in% tips_present]
  
  list(
    fortified_trees = fortified_trees,
    tip_order = tip_order_for_plot
  )
}

mirror_fast_fortified_trees <- function(fortified_trees) {
  
  lapply(fortified_trees, function(df) {
    
    if ("x" %in% names(df)) {
      df$x <- -df$x
    }
    
    if ("branch" %in% names(df)) {
      df$branch <- -df$branch
    }
    
    if ("xend" %in% names(df)) {
      df$xend <- -df$xend
    }
    
    df
  })
}

# ------------------------------------------------------------
# 14.3 Single fast densitree plot
# ------------------------------------------------------------

make_fast_ggdensitree_plot <- function(fortified_trees,
                                       tip_order,
                                       title,
                                       palette,
                                       side = c("right", "left"),
                                       show_legend = TRUE) {
  
  side <- match.arg(side)
  
  plot_trees <- fortified_trees
  
  if (side == "left") {
    plot_trees <- mirror_fast_fortified_trees(plot_trees)
  }
  
  x_vals <- unlist(lapply(plot_trees, function(z) z$x))
  x_vals <- x_vals[is.finite(x_vals)]
  
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  
  y_vals <- unlist(lapply(plot_trees, function(z) z$y))
  y_vals <- y_vals[is.finite(y_vals)]
  
  y_min <- min(y_vals, na.rm = TRUE)
  y_max <- max(y_vals, na.rm = TRUE)
  
  if (side == "left") {
    x_limits <- c(x_min * 1.03, 0)
    use_negative_geo <- TRUE
  } else {
    x_limits <- c(0, x_max * 1.03)
    use_negative_geo <- FALSE
  }
  
  p <- suppressWarnings(
    suppressMessages(
      ggtree::ggdensitree(
        plot_trees,
        ggplot2::aes(colour = state),
        continuous = "colour",
        tip.order = tip_order,
        align.tips = FALSE,
        layout = "rectangular",
        jitter = 0.1
      )
    )
  ) +
    ggtree::geom_tiplab(
      colour = "transparent",
      size = 1.5
    ) +
    ggplot2::scale_colour_gradientn(
      colours = palette,
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("0", "0.5", "1"),
      oob = scales::squish,
      name = "P(presence)"
    ) +
    ggplot2::scale_x_continuous(
      limits = x_limits,
      labels = function(x) abs(x),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(y_min - 2, y_max + 2),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = title,
      x = "Age Ma",
      y = NULL
    ) +
    deeptime::coord_geo(
      dat = list("periods", "epochs"),
      pos = list("bottom", "bottom"),
      # NUEVO: Aumentamos el grosor de las barras de la escala (antes 0.30 y 0.34)
      height = list(grid::unit(0.60, "cm"), grid::unit(0.70, "cm")),
      # NUEVO: Aumentamos el tamaño de la letra de la escala (antes 2.2 y 1.8)
      size = list(4.0, 3.5),
      neg = use_negative_geo,
      abbrv = TRUE,
      expand = FALSE
    ) +
    ggtree::theme_tree2(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10, colour = "black"),
      axis.title.x = ggplot2::element_text(size = 12, colour = "black"),
      legend.position = ifelse(show_legend, "bottom", "none"),
      legend.title = ggplot2::element_text(face = "bold", size = 11, colour = "black"),
      legend.text = ggplot2::element_text(size = 10, colour = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 15, colour = "black"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      # NUEVO: Aumentamos el margen inferior de 38 a 70 para hacer sitio a la escala más grande
      plot.margin = ggplot2::margin(8, 8, 70, 8) 
    )
  
  p
}

# ------------------------------------------------------------
# 14.4 Central labels panel
# ------------------------------------------------------------

make_central_tip_label_panel <- function(tip_order,
                                         title = NULL,
                                         reverse_order = FALSE) {
  
  label_order <- tip_order
  
  if (isTRUE(reverse_order)) {
    label_order <- rev(label_order)
  }
  
  label_df <- data.frame(
    label = label_order,
    y = seq_along(label_order),
    stringsAsFactors = FALSE
  )
  
  label_df$label_plot <- pretty_taxon_name(label_df$label)
  
  ggplot2::ggplot(label_df, ggplot2::aes(x = 0, y = .data$y)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = -0.45, xend = -0.10, yend = .data$y),
      linetype = "dashed",
      linewidth = 0.25,
      colour = "grey55"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0.10, xend = 0.45, yend = .data$y),
      linetype = "dashed",
      linewidth = 0.25,
      colour = "grey55"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label_plot),
      size = 3.1,
      fontface = "italic",
      colour = "black"
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.5, length(label_order) + 0.5),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-0.55, 0.55),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12),
      # NUEVO: Aumentamos el margen inferior a 70 para que se alinee con los árboles
      plot.margin = ggplot2::margin(8, 2, 70, 2),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# ------------------------------------------------------------
# 14.5 Facing densitree with fast ggdensitree panels
# ------------------------------------------------------------

# ============================================================
# PATCH: get fixed tip order from reference tree
# Needed by make_fast_facing_probability_pair()
# ============================================================

get_tip_order_from_tree <- function(tree, ladderize_tree = TRUE, reverse = FALSE) {
  
  if (inherits(tree, "multiPhylo")) {
    tree <- tree[[1]]
  }
  
  if (!inherits(tree, "phylo")) {
    stop("base_order_tree must be a phylo or multiPhylo object.")
  }
  
  # Clean labels if the cleaning function exists
  if (exists("clean_tree_tip_labels")) {
    tree <- clean_tree_tip_labels(tree)
  }
  
  # Fix branch lengths/topology if the helper exists
  if (exists("fix_tree_for_analysis")) {
    tree <- fix_tree_for_analysis(tree)
  }
  
  if (isTRUE(ladderize_tree)) {
    tree <- ape::ladderize(tree, right = FALSE)
  }
  
  tip_order <- tree$tip.label
  
  if (isTRUE(reverse)) {
    tip_order <- rev(tip_order)
  }
  
  tip_order
}

make_fast_facing_probability_pair <- function(dataset_label,
                                              trees,
                                              asr_results,
                                              dat,
                                              taxon_col,
                                              base_order_tree,
                                              right_trait) {
  
  tip_order <- get_tip_order_from_tree(base_order_tree)
  
  left_data <- build_fast_probability_fortified_trees(
    trees = trees,
    asr_results = asr_results,
    dat = dat,
    taxon_col = taxon_col,
    trait_name = "Echolocation",
    tip_order = tip_order
  )
  
  right_data <- build_fast_probability_fortified_trees(
    trees = trees,
    asr_results = asr_results,
    dat = dat,
    taxon_col = taxon_col,
    trait_name = right_trait,
    tip_order = tip_order
  )
  
  tip_order_plot <- left_data$tip_order
  tip_order_plot <- tip_order_plot[tip_order_plot %in% right_data$tip_order]
  
  p_left <- make_fast_ggdensitree_plot(
    fortified_trees = left_data$fortified_trees,
    tip_order = tip_order_plot,
    title = "Echolocation",
    palette = pal_echo,
    side = "left",
    show_legend = TRUE
  )
  
  p_mid <- make_central_tip_label_panel(
    tip_order = tip_order_plot,
    title = dataset_label,
    reverse_order = REVERSE_CENTRAL_LABEL_ORDER
  )
  
  p_right <- make_fast_ggdensitree_plot(
    fortified_trees = right_data$fortified_trees,
    tip_order = tip_order_plot,
    title = gsub("_", " ", right_trait),
    palette = pal_cave,
    side = "right",
    show_legend = TRUE
  )
  
  p_left + p_mid + p_right +
    patchwork::plot_layout(widths = c(1.05, 0.58, 1.05), guides = "keep") +
    patchwork::plot_annotation(
      title = paste0(dataset_label, " - probability density trees"),
      subtitle = paste0(
        "Left: echolocation | Centre: species labels | Right: ",
        gsub("_", " ", right_trait)
      )
    ) &
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# ------------------------------------------------------------
# 14.6 Generate and save fast densitrees
# ------------------------------------------------------------

message("Building fast ggdensitree figures...")

jones_pair_cave <- make_fast_facing_probability_pair(
  dataset_label = "Jones",
  trees = jones_trees,
  asr_results = jones_asr,
  dat = dat,
  taxon_col = jones_taxon_col,
  base_order_tree = jones_mean_tree,
  right_trait = "Cave_deposit"
)

jones_pair_karst <- make_fast_facing_probability_pair(
  dataset_label = "Jones",
  trees = jones_trees,
  asr_results = jones_asr,
  dat = dat,
  taxon_col = jones_taxon_col,
  base_order_tree = jones_mean_tree,
  right_trait = "Cave_karstic_deposit"
)

hand_pair_cave <- make_fast_facing_probability_pair(
  dataset_label = "Hand",
  trees = hand_trees,
  asr_results = hand_asr,
  dat = dat,
  taxon_col = hand_taxon_col,
  base_order_tree = hand_mcc_tree,
  right_trait = "Cave_deposit"
)

hand_pair_karst <- make_fast_facing_probability_pair(
  dataset_label = "Hand",
  trees = hand_trees,
  asr_results = hand_asr,
  dat = dat,
  taxon_col = hand_taxon_col,
  base_order_tree = hand_mcc_tree,
  right_trait = "Cave_karstic_deposit"
)

# ============================================================
# PATCH: save plots safely
# Needed by save_plot_all_formats()
# ============================================================

if (!exists("DENSITREE_WIDTH_MM")) {
  DENSITREE_WIDTH_MM <- 520
}

if (!exists("DENSITREE_HEIGHT_MM")) {
  DENSITREE_HEIGHT_MM <- 360
}

if (!exists("DENSITREE_DPI")) {
  DENSITREE_DPI <- 300
}

save_plot_all_formats <- function(p,
                                  file_base,
                                  width_mm = 520,
                                  height_mm = 360,
                                  dpi = 300,
                                  save_pdf = TRUE,
                                  save_png = TRUE,
                                  save_tiff = FALSE) {
  
  out_folder <- dirname(file_base)
  dir.create(out_folder, recursive = TRUE, showWarnings = FALSE)
  
  if (save_pdf) {
    message("Saving PDF: ", paste0(file_base, ".pdf"))
    
    ggplot2::ggsave(
      filename = paste0(file_base, ".pdf"),
      plot = p,
      width = width_mm,
      height = height_mm,
      units = "mm",
      bg = "white",
      limitsize = FALSE
    )
  }
  
  if (save_png) {
    message("Saving PNG: ", paste0(file_base, ".png"))
    
    ggplot2::ggsave(
      filename = paste0(file_base, ".png"),
      plot = p,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
  
  if (save_tiff) {
    message("Saving TIFF: ", paste0(file_base, ".tiff"))
    
    ggplot2::ggsave(
      filename = paste0(file_base, ".tiff"),
      plot = p,
      width = width_mm,
      height = height_mm,
      units = "mm",
      dpi = dpi,
      compression = "lzw",
      bg = "white",
      limitsize = FALSE
    )
  }
  
  invisible(file_base)
}

save_plot_all_formats(
  jones_pair_cave,
  file.path(densitree_dir, "Jones_FAST_ggdensitree_Echolocation_vs_Cave_deposit"),
  width_mm = DENSITREE_WIDTH_MM,
  height_mm = DENSITREE_HEIGHT_MM,
  dpi = DENSITREE_DPI
)


save_plot_all_formats(
  jones_pair_karst,
  file.path(densitree_dir, "Jones_FAST_ggdensitree_Echolocation_vs_Cave_karstic_deposit"),
  width_mm = DENSITREE_WIDTH_MM,
  height_mm = DENSITREE_HEIGHT_MM,
  dpi = DENSITREE_DPI
)

save_plot_all_formats(
  hand_pair_cave,
  file.path(densitree_dir, "Hand_FAST_ggdensitree_Echolocation_vs_Cave_deposit"),
  width_mm = DENSITREE_WIDTH_MM,
  height_mm = DENSITREE_HEIGHT_MM,
  dpi = DENSITREE_DPI
)

save_plot_all_formats(
  hand_pair_karst,
  file.path(densitree_dir, "Hand_FAST_ggdensitree_Echolocation_vs_Cave_karstic_deposit"),
  width_mm = DENSITREE_WIDTH_MM,
  height_mm = DENSITREE_HEIGHT_MM,
  dpi = DENSITREE_DPI
)

composite_cave <- jones_pair_cave / hand_pair_cave +
  patchwork::plot_annotation(
    title = "Echolocation vs cave deposit probability density trees",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 28, face = "bold")) # <--- TAMAÑO SEGURO: 28
  )

composite_karst <- jones_pair_karst / hand_pair_karst +
  patchwork::plot_annotation(
    title = "Echolocation vs cave/karstic deposit probability density trees",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 28, face = "bold")) # <--- TAMAÑO SEGURO: 28
  )

save_plot_all_formats(
  composite_cave,
  file.path(densitree_dir, "COMPOSITE_FAST_ggdensitree_Jones_Hand_Echolocation_vs_Cave_deposit"),
  width_mm = DENSITREE_WIDTH_MM,
  height_mm = DENSITREE_HEIGHT_MM * 2,
  dpi = DENSITREE_DPI
)

save_plot_all_formats(
  composite_karst,
  file.path(densitree_dir, "COMPOSITE_FAST_ggdensitree_Jones_Hand_Echolocation_vs_Cave_karstic_deposit"),
  width_mm = DENSITREE_WIDTH_MM,
  height_mm = DENSITREE_HEIGHT_MM * 2,
  dpi = DENSITREE_DPI
)

# ------------------------------------------------------------
# 15. Save final objects
# ------------------------------------------------------------

saveRDS(hand_asr, file.path(out_dir, "Hand_ASR_integrated_results.rds"))
saveRDS(jones_asr, file.path(out_dir, "Jones_ASR_integrated_results.rds"))
saveRDS(hand_trees, file.path(out_dir, "Hand_selected_trees.rds"))
saveRDS(jones_trees, file.path(out_dir, "Jones_selected_dated_trees.rds"))
saveRDS(hand_reference_cave, file.path(out_dir, "Hand_MCC_reference_ASR_Echolocation_vs_Cave_deposit.rds"))
saveRDS(hand_reference_karst, file.path(out_dir, "Hand_MCC_reference_ASR_Echolocation_vs_Cave_karstic_deposit.rds"))
saveRDS(jones_reference_cave, file.path(out_dir, "Jones_mean_age_reference_ASR_Echolocation_vs_Cave_deposit.rds"))
saveRDS(jones_reference_karst, file.path(out_dir, "Jones_mean_age_reference_ASR_Echolocation_vs_Cave_karstic_deposit.rds"))

message("\nDONE.")
message("Outputs saved in:")
message(normalizePath(out_dir, winslash = "/"))
