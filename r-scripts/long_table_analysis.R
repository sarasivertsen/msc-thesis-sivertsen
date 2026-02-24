# README
# How to run:
#   Rscript: long_table_analysis.R, dataset: fd_long.csv
# Workflow:
# 1. PERMANOVA (Bray-Curtis distances, avg. linkage) 
# 2. Targeted GLMM (logistic, random participant)
# Notes:
#   - The script installs relevant packages if needed 
#   - Warnings/errors are written to long_analysis_log.txt (not printed)
#   - Outputs go to an output folder

# ----- Installing packages, setting up directories, defining datasets, storing results ----- #

# Installing packages 
required_pkgs <- c("vegan", "lme4", "stats", "methods")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

# Creating directories and output folder
data_dir <- "flavor_data"
results_dir <- "flavor_results"

if (!dir.exists(data_dir)) {
  stop("Data directory does not exist: ", normalizePath(data_dir, mustWork = FALSE))
}

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

out_path <- function(filename) file.path(results_dir, filename)

# Loading input data 
long_file <- file.path(data_dir, "fd_long.csv")
fisher_file <- file.path(results_dir, "fisher_flavor_additive_pairwise.csv")

if (!file.exists(long_file)) {
  stop("Input file not found: ", long_file)
}

# Reading datasets
long_df <- read.csv(long_file, check.names = FALSE, stringsAsFactors = FALSE)
required_cols <- c("participant_id", "additive", "flavors", "gender", "covid_status")
missing_cols <- setdiff(required_cols, names(long_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

long_df <- long_df[!(is.na(long_df$additive) | is.na(long_df$participant_id)), ]
long_df$gender[is.na(long_df$gender)] <- "Unknown"
long_df$covid_status[is.na(long_df$covid_status)] <- "Unknown"
long_df$flavors[is.na(long_df$flavors)] <- ""

# Writing logs 
log_msgs <- c()
capture_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    log_msgs <<- c(log_msgs, paste("Warning:", conditionMessage(w)))
    invokeRestart("muffleWarning")
  })
}

# ---------- Preparing the datasets for PERMANOVA and GLMM analyses ---------- #

# PERMANOVA
split_flavors <- strsplit(long_df$flavors, ";")
split_flavors <- lapply(split_flavors, function(x) trimws(x[x != ""]))
all_flavs <- sort(unique(unlist(split_flavors)))

rows <- data.frame(
  participant_id = trimws(long_df$participant_id),
  additive = trimws(long_df$additive),
  gender = trimws(long_df$gender),
  covid_status = trimws(long_df$covid_status),
  stringsAsFactors = FALSE
)

mat <- matrix(
  0,
  nrow = nrow(long_df),
  ncol = length(all_flavs),
  dimnames = list(NULL, all_flavs))

for (i in seq_len(nrow(long_df))) {
  if (length(split_flavors[[i]]) > 0) {
    mat[i, split_flavors[[i]]] <- 1
  }
}

row_sums <- rowSums(mat)
if (any(row_sums == 0)) {
  dropped <- sum(row_sums == 0)
  log_msgs <- c(log_msgs, paste0("Dropped ", dropped, " rows with no reported flavors."))
  keep <- row_sums > 0
  mat <- mat[keep, , drop = FALSE]
  rows <- rows[keep, , drop = FALSE]
}

# GLMM 
glmm_rows <- data.frame(
  participant_id = rows$participant_id,
  additive = rows$additive,
  gender = rows$gender,
  covid_status = rows$covid_status,
  stringsAsFactors = FALSE
)
glmm_results <- list()
glmm_df <- data.frame()

# Determining target flavors from Fisher's exact output file (odds_ratio > 1 and BH < 0.05)
target_flavs <- all_flavs
if (file.exists(fisher_file)) {
  fisher_df <- try(read.csv(fisher_file, stringsAsFactors = FALSE), silent = TRUE)
  if (!inherits(fisher_df, "try-error")) {
    target_flavs <- unique(fisher_df$flavor[(fisher_df$odds_ratio > 1) & (fisher_df$p_adj_BH < 0.05)])
    log_msgs <- c(log_msgs, paste0("Using ", length(target_flavs), " flavors from Fisher hits in ", fisher_file, "."))
  }
}

# ----- Running PERMANOVA and GLMM ----- #

# PERMANOVA (Bray-Curtis distances, avg. linkage) 
rows[] <- lapply(rows, factor)
dist_bc <- vegan::vegdist(mat, method = "bray", na.rm = TRUE)

permanova_res <- vegan::adonis2(
  dist_bc ~ additive + gender + covid_status,
  data = rows,
  permutations = 999,
  strata = rows$participant_id,
  by = "margin"
)

write.csv(as.data.frame(permanova_res),
          out_path("permanova_results.csv")
)

# GLMM (logistic with random participant) for targeted flavors 
if (length(target_flavs) > 0) {
  capture_warnings({
    for (flav in target_flavs) {
      y <- mat[, flav]
      if (length(unique(y)) < 2) next 
      dat <- cbind(glmm_rows, y = y)
      dat$additive <- factor(dat$additive)
      dat$gender <- factor(dat$gender)
      dat$covid_status <- factor(dat$covid_status)
      dat$participant_id <- factor(dat$participant_id)
      ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5),
                                 check.conv.singular = "ignore",
                                 check.conv.hess = "ignore",
                                 check.rankX = "ignore")
      m <- try(
        suppressWarnings(
          glmer(y ~ additive + gender + covid_status + (1 | participant_id),
                data = dat, family = binomial(), control = ctrl)
        ),
        silent = TRUE
      )
      fit_mode <- "glmer"
      if (inherits(m, "try-error")) {
        log_msgs <- c(log_msgs, paste0("GLMM failed for flavor ", flav, ": ", conditionMessage(attr(m, "condition")),
                                       ". Falling back to plain glm without random effect."))
        m <- try(suppressWarnings(glm(y ~ additive + gender + covid_status,
                                      data = dat, family = binomial())), silent = TRUE)
        fit_mode <- "glm_fallback"
      }
      if (inherits(m, "try-error")) {
        log_msgs <- c(log_msgs, paste0("GLM fallback failed for flavor ", flav, ": ", conditionMessage(attr(m, "condition"))))
        next
      }
      is_sing <- FALSE
      if (fit_mode == "glmer") {
        is_sing <- lme4::isSingular(m, tol = 1e-4)
      }
      sm <- try(suppressWarnings(summary(m)), silent = TRUE)
      if (inherits(sm, "try-error")) {
        log_msgs <- c(log_msgs, paste0("Summary failed for flavor ", flav, ": ", conditionMessage(attr(sm, "condition"))))
        next
      }
      coefs <- sm$coefficients
      glmm_results[[length(glmm_results) + 1]] <- data.frame(
        flavor = flav,
        term = rownames(coefs),
        estimate = coefs[, "Estimate"],
        std_error = coefs[, "Std. Error"],
        z_value = coefs[, "z value"],
        p_value = coefs[, "Pr(>|z|)"],
        singular = is_sing,
        fit_mode = fit_mode,
        stringsAsFactors = FALSE
      )
    }
  })
  if (length(glmm_results) > 0) {
    glmm_df <- do.call(rbind, glmm_results)
    write.csv(glmm_df, out_path("glmm_flavor_additive_results.csv"), row.names = FALSE)
  }
}

cat("Outputs written:\n",
    "- ", out_path("permanova_results.csv"), " (additive/gender/covid marginal effects)\n",
    if (nrow(glmm_df) > 0) paste0("- ", out_path("glmm_flavor_additive_results.csv"),
                                  " (target flavors from Fisher hits)\n") else "",
    sep = "")

# ---------- Writing log messages (if any) ---------- #
if (length(log_msgs) > 0) {
  writeLines(log_msgs, out_path("long_analysis_log.txt"))
}

# ----- Session information ----- #
writeLines(
  capture.output(sessionInfo()),
  out_path("sessionInfo_long.txt")
)