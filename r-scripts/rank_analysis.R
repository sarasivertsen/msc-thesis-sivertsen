# README
# How to run:
#   R-script: rank_analysis.R, datasets: rank_v.csv, rank_aa.csv, rank_rk.csv, rank_ia.csv 
# Workflow:
#   1. Kendall's W test to test panel agreement 
#   2. Friedman tests (overall and subgroups) 
#   3. Nemenyi post-hoc analysis (if significant Friedman) 
#   4. Exact binomial tests 

# ----- Installing packages, setting up directories, defining datasets, storing results ----- #

# Installing packages 
required_pkgs <- c("tidyverse", "DescTools", "PMCMRplus")

for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}

# Creating directories and output folder
data_dir <- "rank_data"
results_dir <- "rank_results"

if (!dir.exists(data_dir)) {
  stop("Data directory does not exist: ", normalizePath(data_dir, mustWork = FALSE))
}

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
  
out_path <- function(filename) file.path(results_dir, filename)

# Loading datasets
get_data_file <- function(filename) {
  f <- file.path(data_dir, filename)
  if (!file.exists(f)) stop("File not found: ", f)
  f
}

datasets <- list(
  vanillin = get_data_file("rank_v.csv"),           
  acetic_acid = get_data_file("rank_aa.csv"),        
  raspberry_ketone = get_data_file("rank_rk.csv"),   
  isoamyl_acetate = get_data_file("rank_ia.csv")     
)

# Defining subgroups 
subgroups <- list(
  all = function(df) df,
  male = function(df) subset(df, Gender == "male"),
  female = function(df) subset(df, Gender == "female"),
  COVID_yes = function(df) subset(df, COVID == "yes"),
  COVID_no = function(df) subset(df, COVID == "no")
)

# Storing results 
kendall_results  <- list()
friedman_results <- list()
posthoc_results  <- list()
threshold_results <- list()
log_msgs <- c()


# ----- Kendall's W test ----- #   

for (dataset_name in names(datasets)) {
  df <- read.csv(datasets[[dataset_name]], stringsAsFactors = FALSE)
  df$rank <- as.numeric(df$Rank)
  df$Gender <- trimws(tolower(df$Gender))
  df$COVID <- trimws(tolower(df$COVID))
  
  kendall_w <- tryCatch({
    rank_mat <- df %>%
      select(Participant, Sample, Rank) %>%
      pivot_wider(names_from = Sample, values_from = Rank) %>%
      drop_na()
    
    rank_matrix <- as.matrix(rank_mat[,-1])
    res <- DescTools::KendallW(rank_matrix, correct = TRUE)
    
    list(
      W = unname(if (is.list(res)) res$W else res),
      Chi.sq = if (is.list(res)) res$Chi.sq else NA_real_,
      df = if (is.list(res)) res$df else NA_real_,
      p.value = if (is.list(res)) res$p.value else NA_real_
    )
  }, error = function(e) {
    log_msgs <<- c(log_msgs, paste("KendallW failed:", dataset_name, e$message))
    return(NULL)
  })
  
  kendall_results[[dataset_name]] <- kendall_w
}

# ----- Friedman test and Nemenyi post-hoc analysis ----- #

for (dataset_name in names(datasets)) {
  df <- read.csv(datasets[[dataset_name]], stringsAsFactors = FALSE)
  df$Gender <- trimws(tolower(df$Gender))
  df$COVID <- trimws(tolower(df$COVID))
  
  for (subgroup_name in names(subgroups)) {
    sub_df <- tryCatch(subgroups[[subgroup_name]](df), error = function(e) NULL)
      
    if (!is.data.frame(sub_df) || nrow(sub_df) < 2) {
      log_msgs <- c(log_msgs, paste("Skipping", dataset_name, subgroup_name, " - not enough data."))
      next
  }
    
    key <- paste(dataset_name, subgroup_name, sep = "_")
    
    # Friedman test
    fried_test <- tryCatch(
      friedman.test(Rank ~ Sample | Participant, data = sub_df),
      error = function(e) {
        log_msgs <<- c(log_msgs, paste("Friedman failed:", key, e$message))
        return(NULL)
      }
    )
    
    friedman_results[[key]] <- fried_test
    
    # Nemenyi post-hoc (if significant Friedman)
    posthoc <- NULL
    if (inherits(fried_test, "htest") &&
        fried_test$p.value < 0.05 &&
        length(unique(sub_df$Sample)) >= 3) {
      
      posthoc <- tryCatch(
        frdAllPairsNemenyiTest(Rank ~ Sample | Participant, data = sub_df),
        error = function(e) {
          log_msgs <<- c(log_msgs, paste("Nemenyi failed:", key, e$message))
          return(NULL)
        }
      )
    }
    posthoc_results[[key]] <- posthoc
  }
}

# ----- Threshold discrimination analysis ----- #
threshold_result <- list()

for (dataset_name in names(datasets)) {
  df <- read.csv(datasets[[dataset_name]], stringsAsFactors = FALSE)
  df$Gender <- trimws(tolower(df$Gender))
  df$COVID <- trimws(tolower(df$COVID))
  
  for (subgroup_name in names(subgroups)) {
    sub_df <- tryCatch(subgroups[[subgroup_name]](df), error = function(e) NULL)
    if (!is.data.frame(sub_df) || nrow(sub_df) == 0) {
      log_msgs <- c(log_msgs, paste("Skipping", dataset_name, subgroup_name, "- no data"))
      next 
    }

    sub_df <- sub_df %>%
      mutate(threshold_group = case_when(
        Sample %in% c("REF", "LOW") ~ "below",
        Sample %in% c("MED", "HIGH") ~ "above"
      ))
    
    participant_scores <- sub_df %>%
      group_by(Participant) %>%
      summarise(
        correct = as.integer(
          all(Rank[threshold_group == "below"] <= 2) &
          all(Rank[threshold_group == "above"] >= 3)
        ),
        .groups = "drop"
      )
    
   bt <- binom.test(
     x = sum(participant_scores$correct),
     n = nrow(participant_scores),
     p = 0.5,
     alternative = "greater"
    )
  
   threshold_result[[paste(dataset_name, subgroup_name, sep = "_")]] <- data.frame(
     dataset = dataset_name,
     subgroup = subgroup_name,
     n_participants = nrow(participant_scores),
     n_correct = sum(participant_scores$correct), 
     proportion_correct = mean(participant_scores$correct),
     p_value = bt$p.value
    )
  }
}


# ----- Saving the results ----- #

# Kendall's W
kendall_summary <- bind_rows(
  lapply(names(kendall_results), function(nm) {
    x <- kendall_results[[nm]]
    if (is.null(x)) return(NULL)
    data.frame(
      dataset = nm,
      W = x$W,
      chi_sq = x$Chi.sq,
      df = x$df,
      p_value = x$p.value
    )
  })
)

write.csv(kendall_summary, out_path("kendall_W_results.csv"), row.names = FALSE)

# Friedman test
friedman_summary <- bind_rows(lapply(names(friedman_results), function(nm) {
  x <- friedman_results[[nm]]
  if (!inherits(x, "htest")) return(NULL)
  data.frame(
    dataset_subgroup = nm,
    statistic = unname(x$statistic),
    df = unname(x$parameter),
    p_value = x$p.value
  )
}))

write.csv(friedman_summary, out_path("friedman_results_summary.csv"), row.names = FALSE)

# Nemenyi post-hoc analysis
posthoc_summary <- bind_rows(lapply(names(posthoc_results), function(nm) {
  x <- posthoc_results[[nm]]
  if (is.null(x)) return(NULL)
  as.data.frame(as.table(x$p.value)) %>%
    rename(Sample_1 = Var1, Sample_2 = Var2, p_value = Freq) %>%
    mutate(dataset_subgroup = nm)
}))

write.csv(posthoc_summary, out_path("posthoc_nemenyi_results.csv"), row.names = FALSE)

# Threshold discrimination analysis 
threshold_summary <- bind_rows(threshold_result)
write.csv(threshold_summary, out_path("threshold_discrimination_results.csv"), row.names = FALSE)

# ----- Writing logs ------ #
if (length(log_msgs) > 0) {
  writeLines(log_msgs, file.path(results_dir, "friedman_nemenyi_log.txt"))
}

cat("Analysis complete: Kendall's W, Friedman, Nemenyi, and Threshold discrimination results saved.\n")

# ----- Session information ----- #
writeLines(
  capture.output(sessionInfo()),
  out_path("sessionInfo_rank.txt")
)
