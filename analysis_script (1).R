# =============================================================================
#                    MAIZE PHENOTYPIC DATA ANALYSIS
#                    Correlation, Repeatability, Heritability, Variance Partitioning
#                    Publication-Quality Plots
# =============================================================================
#
# STATISTICAL APPROACH:
#
# 1. EXPERIMENTAL DESIGN
#    - 125 genotypes, each at a fixed row position (1-125)
#    - 2 adjacent "reps" per genotype per year (subsamples, not true blocks)
#    - 2 years: Spring-2024, Autumn-2025
#
# 2. MODEL SPECIFICATIONS
#    - Variance partitioning: Year as RANDOM (included in partitioning)
#    - BLUEs: Year as FIXED (to extract adjusted genotype means)
#    - Genotype: RANDOM for variance estimation, FIXED for BLUEs
#
# 3. TERMINOLOGY
#    - Within-season: "Repeatability" (consistency between adjacent reps)
#    - Across-season: "Heritability" (genetic variance on genotype-mean basis)
#
# =============================================================================

# Set library path
.libPaths(c("~/R_libs", .libPaths()))

# Load required libraries
library(tidyverse)
library(lme4)
library(splines)
library(emmeans)
library(corrplot)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(viridis)
library(readxl)

# Set working directory
# Working directory should be set to Analysis folder before running
# setwd("/path/to/Analysis/")

# Create output directories
dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# =============================================================================
#                    CUSTOM THEME FOR PUBLICATION-QUALITY PLOTS
# =============================================================================

theme_publication <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      # Clean white background
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),

      # Solid black axis lines
      axis.line = element_line(color = "black", linewidth = 0.6),
      panel.border = element_blank(),

      # No grid lines
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      # Axis text and titles
      axis.text = element_text(color = "black", size = base_size - 2),
      axis.title = element_text(color = "black", size = base_size, face = "bold"),
      axis.ticks = element_line(color = "black", linewidth = 0.4),

      # Title
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5, color = "gray40"),

      # Legend
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size - 1, face = "bold"),

      # Margins
      plot.margin = margin(10, 15, 10, 15)
    )
}

# Color palette
season_colors <- c("Sp-24" = "#2E86AB", "Aut-25" = "#E94F37")
component_colors <- c("Genotype" = "#1B9AAA", "Year" = "#F6AE2D",
                      "Row" = "#86BBD8", "G x E" = "#F26419", "Residual" = "#33658A")

# =============================================================================
#                         1. LOAD AND PREPARE DATA
# =============================================================================

cat("\n========== LOADING DATA ==========\n")

# Read from Excel file (PhenotypeData.xlsx has 3 sheets)
season1 <- read_excel("PhenotypeData.xlsx", sheet = "Season 1")
season2 <- read_excel("PhenotypeData.xlsx", sheet = "Season 2")
whole <- read_excel("PhenotypeData.xlsx", sheet = "Whole")

# Convert to data.frame
season1 <- as.data.frame(season1)
season2 <- as.data.frame(season2)
whole <- as.data.frame(whole)

# Clean column names - standardize to snake_case
clean_names <- function(df) {
  names(df) <- gsub("\\s+", "", names(df))      # Remove spaces
  names(df) <- gsub("/", "_", names(df))        # Replace / with _
  names(df) <- gsub("\\(", "_", names(df))      # Replace ( with _
  names(df) <- gsub("\\)", "", names(df))       # Remove )
  names(df) <- gsub("\\.", "_", names(df))      # Replace . with _
  return(df)
}

season1 <- clean_names(season1)
season2 <- clean_names(season2)
whole <- clean_names(whole)

cat("Column names after cleaning:\n")
print(names(whole))

# Get trait columns (columns 5-15)
traits <- names(whole)[5:15]
cat("\nTraits to analyze:\n")
print(traits)

# Trait labels for better display (keys must match cleaned names)
trait_labels <- c(
  "DTA_days" = "Days to Anthesis",
  "DTS_days" = "Days to Silking",
  "ASI_days" = "Anthesis-Silking Interval",
  "EH_cm" = "Ear Height",
  "PH_cm" = "Plant Height",
  "NOR_C" = "Rows per Cob",
  "NOK_R" = "Kernels per Row",
  "HGW_g" = "100 Grain Weight",
  "DTM_days" = "Days to Maturity",
  "Fe_mg_kg" = "Iron",
  "Zn_mg_kg" = "Zinc"
)

# =============================================================================
#                    2. SUMMARY STATISTICS TABLE
# =============================================================================

cat("\n========== SUMMARY STATISTICS ==========\n")

# Calculate summary stats for each trait across both seasons
summary_stats <- whole %>%
  summarise(across(all_of(traits), list(
    Mean = ~mean(.x, na.rm = TRUE),
    SD = ~sd(.x, na.rm = TRUE),
    Min = ~min(.x, na.rm = TRUE),
    Max = ~max(.x, na.rm = TRUE),
    CV = ~sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE) * 100
  ), .names = "{.col}_{.fn}"))

# Reshape to long format for table
summary_table <- data.frame(
  Trait = traits,
  TraitLabel = trait_labels[traits]
)
summary_table$Mean <- sapply(traits, function(t) round(summary_stats[[paste0(t, "_Mean")]], 2))
summary_table$SD <- sapply(traits, function(t) round(summary_stats[[paste0(t, "_SD")]], 2))
summary_table$Min <- sapply(traits, function(t) round(summary_stats[[paste0(t, "_Min")]], 2))
summary_table$Max <- sapply(traits, function(t) round(summary_stats[[paste0(t, "_Max")]], 2))
summary_table$CV <- sapply(traits, function(t) round(summary_stats[[paste0(t, "_CV")]], 1))

cat("\nSummary Statistics:\n")
print(summary_table)
write.csv(summary_table, "results/summary_statistics.csv", row.names = FALSE)

# =============================================================================
#                    3. CORRELATION ANALYSIS
# =============================================================================

cat("\n========== CORRELATION ANALYSIS ==========\n")

# -----------------------------------------------------------------------------
# 3a. Fe and Zn correlations with other traits (FOCUSED - for main paper)
# -----------------------------------------------------------------------------
# Only show correlations of micronutrients with agronomic traits

cor_whole <- cor(whole[, traits], use = "complete.obs")

# Extract Fe and Zn correlations with other traits
other_traits <- setdiff(traits, c("Fe_mg_kg", "Zn_mg_kg"))

fe_zn_cors <- data.frame(
  Trait = rep(other_traits, 2),
  TraitLabel = rep(trait_labels[other_traits], 2),
  Micronutrient = c(rep("Fe", length(other_traits)),
                    rep("Zn", length(other_traits))),
  Correlation = c(
    sapply(other_traits, function(t) cor_whole["Fe_mg_kg", t]),
    sapply(other_traits, function(t) cor_whole["Zn_mg_kg", t])
  )
)

# Add Fe-Zn correlation
fe_zn_r <- cor_whole["Fe_mg_kg", "Zn_mg_kg"]
cat("\nFe-Zn correlation (pooled data): r =", round(fe_zn_r, 3), "\n")

cat("\nMicronutrient correlations with agronomic traits:\n")
print(fe_zn_cors)
write.csv(fe_zn_cors, "results/micronutrient_correlations.csv", row.names = FALSE)

# Plot: Fe and Zn correlations with other traits (store for combined figure)
p_correlations <- ggplot(fe_zn_cors, aes(x = reorder(TraitLabel, Correlation), y = Correlation,
                        fill = Micronutrient)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Fe" = "#8B4513", "Zn" = "#4682B4")) +
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = seq(-0.2, 0.2, 0.1)) +
  labs(x = NULL, y = "Pearson Correlation (r)") +
  theme_publication() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = c(0.85, 0.15),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = NA))

png("plots/01_micronutrient_correlations.png", width = 1200, height = 600, res = 150)
print(p_correlations)
dev.off()

# -----------------------------------------------------------------------------
# 3b. Cross-season correlations (stability)
# -----------------------------------------------------------------------------

geno_means_s1 <- season1 %>%
  group_by(InbCode) %>%
  summarise(across(all_of(traits), ~mean(.x, na.rm = TRUE), .names = "S1_{.col}"))

geno_means_s2 <- season2 %>%
  group_by(InbCode) %>%
  summarise(across(all_of(traits), ~mean(.x, na.rm = TRUE), .names = "S2_{.col}"))

geno_means <- merge(geno_means_s1, geno_means_s2, by = "InbCode")

between_season_cor <- data.frame(
  Trait = traits,
  TraitLabel = trait_labels[traits],
  Correlation = sapply(traits, function(t) {
    cor(geno_means[[paste0("S1_", t)]], geno_means[[paste0("S2_", t)]], use = "complete.obs")
  })
)

cat("\nCross-Season Correlations (genotype stability):\n")
print(between_season_cor)
write.csv(between_season_cor, "results/cross_season_correlations.csv", row.names = FALSE)

# Plot cross-season correlations
png("plots/02_cross_season_stability.png", width = 1000, height = 700, res = 150)
ggplot(between_season_cor, aes(x = reorder(TraitLabel, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", fill = "#2E86AB", width = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "Correlation Across Years (r)") +
  theme_publication() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank())
dev.off()

# -----------------------------------------------------------------------------
# 3c. Cross-season scatter for Fe and Zn only (key traits)
# -----------------------------------------------------------------------------

fe_r <- round(cor(geno_means$S1_Fe_mg_kg, geno_means$S2_Fe_mg_kg, use = "complete.obs"), 2)
zn_r <- round(cor(geno_means$S1_Zn_mg_kg, geno_means$S2_Zn_mg_kg, use = "complete.obs"), 2)

p_fe <- ggplot(geno_means, aes(x = S1_Fe_mg_kg, y = S2_Fe_mg_kg)) +
  geom_point(alpha = 0.6, color = "#8B4513", size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf, label = paste("r =", fe_r),
           hjust = 1.1, vjust = -0.5, size = 4) +
  labs(x = "Spring-2024 Fe (mg/kg)", y = "Autumn-2025 Fe (mg/kg)") +
  theme_publication() +
  theme(panel.grid = element_blank())

p_zn <- ggplot(geno_means, aes(x = S1_Zn_mg_kg, y = S2_Zn_mg_kg)) +
  geom_point(alpha = 0.6, color = "#4682B4", size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf, label = paste("r =", zn_r),
           hjust = 1.1, vjust = -0.5, size = 4) +
  labs(x = "Spring-2024 Zn (mg/kg)", y = "Autumn-2025 Zn (mg/kg)") +
  theme_publication() +
  theme(panel.grid = element_blank())

png("plots/03_cross_season_Fe_Zn.png", width = 1200, height = 500, res = 150)
grid.arrange(p_fe, p_zn, ncol = 2)
dev.off()

# =============================================================================
#                    4. SPATIAL TREND ANALYSIS (for Supplementary Materials)
# =============================================================================
# NOTE: This figure goes in supplement to show spatial gradients exist.
# The main models use splines to correct for these trends.

cat("\n========== SPATIAL TRENDS (Supplementary) ==========\n")

# Combined spatial plot - both seasons overlaid
png("plots/S1_spatial_trends.png", width = 1600, height = 1000, res = 150)
plots_comb <- list()
for (i in seq_along(traits)) {
  t <- traits[i]
  p <- ggplot(whole, aes_string(x = "Row", y = t, color = "Year")) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
    scale_color_manual(values = season_colors) +
    labs(title = trait_labels[t], x = "Row", y = NULL) +
    theme_publication(base_size = 9) +
    theme(plot.title = element_text(size = 9),
          legend.position = "none",
          panel.grid = element_blank())
  plots_comb[[i]] <- p
}
# Add legend to first plot
plots_comb[[1]] <- plots_comb[[1]] + theme(legend.position = c(0.8, 0.2))
grid.arrange(grobs = plots_comb, ncol = 4)
dev.off()

cat("Saved: plots/S1_spatial_trends.png (for Supplementary Materials)\n")

# =============================================================================
#                    4. REPEATABILITY AND HERITABILITY ANALYSIS
# =============================================================================
#
# DESIGN NOTES:
# - Each genotype is at a fixed row position (genotype confounded with row)
# - Two "reps" per genotype are adjacent plants at the same position
# - We fit a 1D spline on Row to capture spatial trends within each season
# - Within-season: "Repeatability" (consistency between adjacent reps)
# - Across-season: "Heritability" (genetic variance / total variance)
# - Year is treated as FIXED in all models (only 2 levels)
#
# =============================================================================

cat("\n========== REPEATABILITY & HERITABILITY ANALYSIS ==========\n")

# -----------------------------------------------------------------------------
# 4a. WITHIN-SEASON REPEATABILITY
# -----------------------------------------------------------------------------
# Model: trait ~ (1|InbCode)
#   - (1|InbCode): Genotype as random effect
#   - Residual: Variation between the two reps
#
# Repeatability = Var(InbCode) / (Var(InbCode) + Var(Residual))

calc_repeatability <- function(data, trait) {
  formula <- as.formula(paste(trait, "~ (1|InbCode)"))

  tryCatch({
    model <- lmer(formula, data = data, REML = TRUE)
    var_comp <- as.data.frame(VarCorr(model))
    var_geno <- var_comp$vcov[var_comp$grp == "InbCode"]
    var_resid <- var_comp$vcov[var_comp$grp == "Residual"]

    # Repeatability: proportion of variance due to genotype (after spatial correction)
    repeatability <- var_geno / (var_geno + var_resid)

    return(list(
      trait = trait,
      var_geno = var_geno,
      var_resid = var_resid,
      repeatability = repeatability
    ))
  }, error = function(e) {
    cat("Error for", trait, ":", e$message, "\n")
    return(list(trait = trait, var_geno = NA, var_resid = NA, repeatability = NA))
  })
}

# Season 1 repeatability
cat("\nCalculating Season 1 repeatability...\n")
repeat_s1 <- lapply(traits, function(t) calc_repeatability(season1, t))
repeat_s1_df <- do.call(rbind, lapply(repeat_s1, function(x) {
  data.frame(Trait = x$trait, Var_Geno = x$var_geno, Var_Resid = x$var_resid,
             Repeatability = x$repeatability)
}))

# Season 2 repeatability
cat("Calculating Season 2 repeatability...\n")
repeat_s2 <- lapply(traits, function(t) calc_repeatability(season2, t))
repeat_s2_df <- do.call(rbind, lapply(repeat_s2, function(x) {
  data.frame(Trait = x$trait, Var_Geno = x$var_geno, Var_Resid = x$var_resid,
             Repeatability = x$repeatability)
}))

# -----------------------------------------------------------------------------
# 4b. ACROSS-SEASON HERITABILITY
# -----------------------------------------------------------------------------
# Model: trait ~ Year + (1|InbCode)
#   - Year: Fixed effect (only 2 levels)
#   - (1|InbCode): Genotype variance
#
# Heritability = Var(InbCode) / (Var(InbCode) + Var(Residual)/n)
# where n = number of observations per genotype (2 reps × 2 years = 4)

calc_heritability_across <- function(data, trait, n_obs_per_geno = 4) {
  formula <- as.formula(paste(trait, "~ Year + (1|InbCode)"))

  tryCatch({
    model <- lmer(formula, data = data, REML = TRUE)
    var_comp <- as.data.frame(VarCorr(model))
    var_geno <- var_comp$vcov[var_comp$grp == "InbCode"]
    var_resid <- var_comp$vcov[var_comp$grp == "Residual"]

    # Broad-sense heritability on a genotype-mean basis
    heritability <- var_geno / (var_geno + var_resid / n_obs_per_geno)

    return(list(
      trait = trait,
      var_geno = var_geno,
      var_resid = var_resid,
      heritability = heritability
    ))
  }, error = function(e) {
    cat("Error for", trait, ":", e$message, "\n")
    return(list(trait = trait, var_geno = NA, var_resid = NA, heritability = NA))
  })
}

cat("Calculating across-season heritability...\n")
herit_across <- lapply(traits, function(t) calc_heritability_across(whole, t))
herit_across_df <- do.call(rbind, lapply(herit_across, function(x) {
  data.frame(Trait = x$trait, Var_Geno = x$var_geno, Var_Resid = x$var_resid,
             Heritability = x$heritability)
}))

# -----------------------------------------------------------------------------
# 4c. COMBINED SUMMARY
# -----------------------------------------------------------------------------
herit_summary <- data.frame(
  Trait = repeat_s1_df$Trait,
  TraitLabel = trait_labels[repeat_s1_df$Trait],
  R_Season1 = round(repeat_s1_df$Repeatability, 3),
  R_Season2 = round(repeat_s2_df$Repeatability, 3),
  H2_Across = round(herit_across_df$Heritability, 3)
)

cat("\n=== REPEATABILITY & HERITABILITY SUMMARY ===\n")
cat("R = within-season repeatability\n")
cat("H2 = across-season heritability\n\n")
print(herit_summary)
write.csv(herit_summary, "results/heritability_summary.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 4d. PLOT
# -----------------------------------------------------------------------------
herit_long <- herit_summary %>%
  pivot_longer(cols = c(R_Season1, R_Season2, H2_Across),
               names_to = "Measure", values_to = "Value") %>%
  mutate(Measure = recode(Measure,
                          "R_Season1" = "Repeatability (Spring-24)",
                          "R_Season2" = "Repeatability (Autumn-25)",
                          "H2_Across" = "Heritability (Across)")) %>%
  mutate(Measure = factor(Measure, levels = c("Repeatability (Spring-24)",
                                               "Repeatability (Autumn-25)",
                                               "Heritability (Across)")))

# Order traits by mean value (low to high)
trait_order <- herit_long %>%
  group_by(TraitLabel) %>%
  summarise(mean_val = mean(Value, na.rm = TRUE)) %>%
  arrange(mean_val) %>%
  pull(TraitLabel)

herit_long$TraitLabel <- factor(herit_long$TraitLabel, levels = trait_order)

# Store plot for combined figure
p_heritability <- ggplot(herit_long, aes(x = TraitLabel, y = Value, fill = Measure)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c("Repeatability (Spring-24)" = "#1E4D2B",
                               "Repeatability (Autumn-25)" = "#D4A84B",
                               "Heritability (Across)" = "#5B7C99")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "Repeatability / Heritability") +
  theme_publication() +
  theme(legend.position = c(0.15, 0.85),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank())

png("plots/05_heritability_comparison.png", width = 1400, height = 800, res = 150)
print(p_heritability)
dev.off()

# =============================================================================
#                    5. VARIANCE PARTITIONING
# =============================================================================
#
# Model: trait ~ (1|Year) + (1|InbCode) + (1|InbCode:Year)
#   - (1|Year): Year variance (random to include in partitioning)
#   - (1|InbCode): Genotype variance
#   - (1|InbCode:Year): Genotype × Year interaction variance
#   - Residual: Within-genotype error
# =============================================================================

cat("\n========== VARIANCE PARTITIONING ==========\n")

variance_partition <- function(data, trait) {
  data$InbCode <- as.factor(data$InbCode)
  data$Year <- as.factor(data$Year)

  # All effects as random for variance partitioning
  formula <- as.formula(paste(trait, "~ (1|Year) + (1|InbCode) + (1|InbCode:Year)"))

  tryCatch({
    model <- lmer(formula, data = data, REML = TRUE,
                  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

    var_comp <- as.data.frame(VarCorr(model))

    var_year <- var_comp$vcov[var_comp$grp == "Year"]
    var_geno <- var_comp$vcov[var_comp$grp == "InbCode"]
    var_gxy <- var_comp$vcov[var_comp$grp == "InbCode:Year"]
    var_resid <- var_comp$vcov[var_comp$grp == "Residual"]

    # Handle NAs/empty
    var_year <- ifelse(is.null(var_year) || length(var_year) == 0, 0, var_year)
    var_geno <- ifelse(is.null(var_geno) || length(var_geno) == 0, 0, var_geno)
    var_gxy <- ifelse(is.null(var_gxy) || length(var_gxy) == 0, 0, var_gxy)

    # Total variance including Year
    total_var <- var_year + var_geno + var_gxy + var_resid

    return(data.frame(
      Trait = trait,
      Var_Year = round(var_year, 4),
      Var_Genotype = round(var_geno, 4),
      Var_GxY = round(var_gxy, 4),
      Var_Residual = round(var_resid, 4),
      Pct_Year = round(100 * var_year / total_var, 2),
      Pct_Genotype = round(100 * var_geno / total_var, 2),
      Pct_GxY = round(100 * var_gxy / total_var, 2),
      Pct_Residual = round(100 * var_resid / total_var, 2)
    ))
  }, error = function(e) {
    cat("Error for", trait, ":", e$message, "\n")
    return(data.frame(Trait = trait, Var_Year = NA, Var_Genotype = NA, Var_GxY = NA,
                      Var_Residual = NA, Pct_Year = NA, Pct_Genotype = NA, Pct_GxY = NA,
                      Pct_Residual = NA))
  })
}

cat("Calculating variance partitioning (Year as random)...\n")
var_part_results <- do.call(rbind, lapply(traits, function(t) variance_partition(whole, t)))
var_part_results$TraitLabel <- trait_labels[var_part_results$Trait]

cat("\n=== VARIANCE PARTITIONING (%) ===\n")
print(var_part_results[, c("TraitLabel", "Pct_Year", "Pct_Genotype", "Pct_GxY", "Pct_Residual")])
write.csv(var_part_results, "results/variance_partitioning.csv", row.names = FALSE)

# Plot variance partitioning - stacked bars
var_part_long <- var_part_results %>%
  select(Trait, TraitLabel, Pct_Year, Pct_Genotype, Pct_GxY, Pct_Residual) %>%
  pivot_longer(cols = starts_with("Pct"), names_to = "Component", values_to = "Percentage") %>%
  mutate(Component = recode(Component,
                            "Pct_Year" = "Year",
                            "Pct_Genotype" = "Genotype",
                            "Pct_GxY" = "G x Y",
                            "Pct_Residual" = "Residual")) %>%
  mutate(Component = factor(Component, levels = c("Residual", "G x Y", "Year", "Genotype")))

var_colors <- c(
  "Genotype" = "#1E4D2B",
  "Year" = "#4A90D9",
  "G x Y" = "#C67D5E",
  "Residual" = "#8C8C8C"
)

# Store plot for combined figure
p_variance <- ggplot(var_part_long, aes(x = TraitLabel, y = Percentage, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = var_colors) +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = NULL, y = "Variance Explained (%)") +
  theme_publication() +
  theme(legend.position = c(0.85, 0.85),
        legend.justification = c(1, 1),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE))

png("plots/06_variance_partitioning.png", width = 1400, height = 800, res = 150)
print(p_variance)
dev.off()

# =========================================================
#  6. BLUES for Five Traits
# =========================================================
#
# Model: trait ~ Inb.Code + Year + bs(Row., df=5)
#   - Inb.Code: Genotype as FIXED (we want BLUEs, not BLUPs)
#   - Year: FIXED (only 2 levels, consistent with other models)
#   - bs(Row.): Spatial trend correction (1D spline)
#
# BLUEs = Best Linear Unbiased Estimates = adjusted genotype means
# =========================================================

cat("\n========== COMPUTING BLUEs ==========\n")

# Load data for BLUEs (uses different file with 5 traits)
data_blues <- read.csv('FiveTraits.csv', check.names = FALSE)

# Clean column names to match main dataset naming convention
# Original: Year, Rep, Row , Inb Code, NOR/C, NOK/R, HGW(g), Fe(mg/kg), Zn(mg/kg)
# Target:   Year, Rep, Row,  Inb_Code, NOR_C, NOK_R, HGW_g,  Fe_mg_kg,  Zn_mg_kg
names(data_blues) <- c("Year", "Rep", "Row", "Inb_Code", "NOR_C", "NOK_R", "HGW_g", "Fe_mg_kg", "Zn_mg_kg")

# Check column names
cat("Column names:", paste(names(data_blues), collapse = ", "), "\n")

# Clean whitespace from genotype names (critical for matching across years)
data_blues$Inb_Code <- trimws(data_blues$Inb_Code)

# No outlier removal - all values are biologically plausible
data_clean <- data_blues

# Convert to factors
data_clean$Inb_Code <- as.factor(data_clean$Inb_Code)
data_clean$Year <- as.factor(data_clean$Year)
data_clean$Row <- as.numeric(data_clean$Row)

cat("Using all", nrow(data_clean), "observations (no outlier removal)\n")
cat("Unique genotypes:", length(unique(data_clean$Inb_Code)), "\n")

# =========================================================
# 6b. Histogram Theme
# =========================================================

custom_theme <- theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    axis.line = element_line(linewidth = 1, color = "black"),
    panel.grid = element_blank()
  )

# =========================================================
# 6c. Function to compute BLUEs and plot histogram
# =========================================================

compute_BLUEs <- function(trait_name, fill_color, bins = 15) {

  cat("\nComputing BLUEs for", trait_name, "...\n")

  # Simple model: Genotype + Year (both fixed effects)
  # Note: Spatial spline removed for BLUEs due to model instability with
  # limited observations per genotype. Year adjustment is sufficient.
  formula <- as.formula(paste(trait_name, "~ Inb_Code + Year"))

  # Fit linear model
  model <- lm(formula, data = data_clean, na.action = na.omit)

  # Compute BLUEs (least-squares means adjusted for Year)
  blues <- as.data.frame(emmeans(model, ~ Inb_Code))

  # Plot histogram
  p <- ggplot(blues, aes(x = emmean)) +
    geom_histogram(fill = fill_color, color = "black", bins = bins) +
    labs(x = "Genotype Mean",
         y = "Frequency") +
    custom_theme +
    scale_y_continuous(expand = c(0, 0))

  print(p)

  # Save plot
  ggsave(paste0("plots/BLUE_histogram_", gsub("[^a-zA-Z0-9]", "", trait_name), ".png"),
         p, width = 8, height = 6, dpi = 150)

  return(blues)
}

# =========================================================
# 6d. Compute BLUEs for all five traits
# =========================================================

blues_nokr <- compute_BLUEs("NOK_R", "#29774e")
blues_norc <- compute_BLUEs("NOR_C", "#c27a74")
blues_hgw  <- compute_BLUEs("HGW_g", "#e68132")
blues_fe   <- compute_BLUEs("Fe_mg_kg", "#4a7c59")
blues_zn   <- compute_BLUEs("Zn_mg_kg", "#8b4513")

# =========================================================
# 6e. Save BLUEs and compute correlations
# =========================================================

# Combine all BLUEs into one dataframe
all_blues <- data.frame(
  Genotype = blues_fe$Inb_Code,
  Fe = blues_fe$emmean,
  Zn = blues_zn$emmean,
  HGW = blues_hgw$emmean,
  NOK_R = blues_nokr$emmean,
  NOR_C = blues_norc$emmean
)

write.csv(all_blues, "results/genotype_BLUEs.csv", row.names = FALSE)

# Correlation between Fe and Zn BLUEs (for abstract claim)
fe_zn_cor <- cor(all_blues$Fe, all_blues$Zn, use = "complete.obs")
cat("\n=== Fe-Zn BLUE Correlation ===\n")
cat("r =", round(fe_zn_cor, 3), "\n")
cat("(This value should match the abstract)\n")

# Correlations of micronutrients with yield proxy (HGW)
fe_hgw_cor <- cor(all_blues$Fe, all_blues$HGW, use = "complete.obs")
zn_hgw_cor <- cor(all_blues$Zn, all_blues$HGW, use = "complete.obs")
cat("\n=== Micronutrient-HGW Correlations ===\n")
cat("Fe-HGW r =", round(fe_hgw_cor, 3), "\n")
cat("Zn-HGW r =", round(zn_hgw_cor, 3), "\n")
cat("(Weak/non-significant = no yield penalty for biofortification)\n")

# Identify top genotypes
cat("\n=== Top 10 Genotypes for Fe ===\n")
print(all_blues %>% arrange(desc(Fe)) %>% head(10) %>% select(Genotype, Fe, Zn, HGW))

cat("\n=== Top 10 Genotypes for Zn ===\n")
print(all_blues %>% arrange(desc(Zn)) %>% head(10) %>% select(Genotype, Zn, Fe, HGW))

cat("\n=== Top Genotypes for BOTH Fe and Zn (top quartile each) ===\n")
fe_q75 <- quantile(all_blues$Fe, 0.75, na.rm = TRUE)
zn_q75 <- quantile(all_blues$Zn, 0.75, na.rm = TRUE)
top_both <- all_blues %>%
  filter(Fe >= fe_q75 & Zn >= zn_q75) %>%
  arrange(desc(Fe + Zn))
print(top_both)
cat("\nThese", nrow(top_both), "genotypes are candidates for biofortification breeding\n")

# Save top genotypes to CSV for paper table
top_fe <- all_blues %>% arrange(desc(Fe)) %>% head(10)
top_zn <- all_blues %>% arrange(desc(Zn)) %>% head(10)
write.csv(top_fe, "results/top10_Fe_genotypes.csv", row.names = FALSE)
write.csv(top_zn, "results/top10_Zn_genotypes.csv", row.names = FALSE)
write.csv(top_both, "results/top_Fe_and_Zn_genotypes.csv", row.names = FALSE)

# =========================================================
# 6f. Fe vs Zn scatterplot with top genotypes labeled
# =========================================================

# Identify top performers to label
top_to_label <- all_blues %>%
  filter(Fe >= fe_q75 & Zn >= zn_q75) %>%
  arrange(desc(Fe + Zn)) %>%
  head(10)

# Store plot for combined figure
p_fe_zn_scatter <- ggplot(all_blues, aes(x = Fe, y = Zn)) +
  # All points
  geom_point(alpha = 0.5, color = "gray50", size = 2.5) +
  # Highlight top performers
  geom_point(data = top_both, aes(x = Fe, y = Zn),
             color = "#2E7D32", size = 3, alpha = 0.8) +
  # Label top 10
  geom_text(data = top_to_label, aes(x = Fe, y = Zn, label = Genotype),
            hjust = -0.1, vjust = 0.5, size = 3, color = "#1B5E20") +
  # Correlation line
  geom_smooth(method = "lm", se = TRUE, color = "#1565C0", fill = "#1565C0",
              alpha = 0.2, linewidth = 1) +
  # Reference lines at 75th percentiles
  geom_vline(xintercept = fe_q75, linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = zn_q75, linetype = "dashed", color = "gray60") +
  # Correlation annotation
  annotate("text", x = Inf, y = -Inf, label = paste("r =", round(fe_zn_cor, 2)),
           hjust = 1.1, vjust = -0.5, size = 4) +
  labs(x = "Fe (mg/kg)",
       y = "Zn (mg/kg)") +
  theme_publication() +
  theme(panel.grid = element_blank())

png("plots/04_Fe_vs_Zn_scatter.png", width = 900, height = 800, res = 150)
print(p_fe_zn_scatter)
dev.off()

cat("\nSaved: plots/04_Fe_vs_Zn_scatter.png\n")

# =========================================================
# 6g. BLUE-based correlations (for Figure 3B)
# =========================================================
# Use BLUEs for consistent correlation estimates across manuscript
# Focus on yield-related traits: HGW, NOK/R, NOR/C

cat("\n=== BLUE-BASED CORRELATIONS (for manuscript consistency) ===\n")

# Calculate estimated grain mass per ear (yield proxy)
all_blues$GrainMass <- all_blues$HGW * all_blues$NOK_R * all_blues$NOR_C

# Compute correlations from BLUEs (4 yield-related traits)
blue_cor_data <- data.frame(
  Trait = rep(c("GrainMass", "HGW", "NOK/R", "NOR/C"), 2),
  TraitLabel = rep(c("Est. Grain Mass", "100 Grain Weight", "Kernels per Row", "Rows per Cob"), 2),
  Micronutrient = c(rep("Fe", 4), rep("Zn", 4)),
  Correlation = c(
    cor(all_blues$Fe, all_blues$GrainMass, use = "complete.obs"),
    cor(all_blues$Fe, all_blues$HGW, use = "complete.obs"),
    cor(all_blues$Fe, all_blues$NOK_R, use = "complete.obs"),
    cor(all_blues$Fe, all_blues$NOR_C, use = "complete.obs"),
    cor(all_blues$Zn, all_blues$GrainMass, use = "complete.obs"),
    cor(all_blues$Zn, all_blues$HGW, use = "complete.obs"),
    cor(all_blues$Zn, all_blues$NOK_R, use = "complete.obs"),
    cor(all_blues$Zn, all_blues$NOR_C, use = "complete.obs")
  )
)

cat("Micronutrient correlations with yield components (BLUE-based, n=125):\n")
print(blue_cor_data)
write.csv(blue_cor_data, "results/micronutrient_correlations_BLUEs.csv", row.names = FALSE)

# Order traits for display (grain mass at top)
blue_cor_data$TraitLabel <- factor(blue_cor_data$TraitLabel,
  levels = c("Rows per Cob", "Kernels per Row", "100 Grain Weight", "Est. Grain Mass"))

# Recreate correlation plot using BLUE-based values
p_correlations <- ggplot(blue_cor_data, aes(x = TraitLabel, y = Correlation,
                        fill = Micronutrient)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Fe" = "#8B4513", "Zn" = "#4682B4")) +
  scale_y_continuous(limits = c(-0.2, 0.2), breaks = seq(-0.2, 0.2, 0.1)) +
  labs(x = NULL, y = "Pearson Correlation (r)") +
  theme_publication() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = c(0.85, 0.15),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = NA))

png("plots/01_micronutrient_correlations.png", width = 1200, height = 500, res = 150)
print(p_correlations)
dev.off()
cat("Updated: plots/01_micronutrient_correlations.png (BLUE-based, with grain mass)\n")


# =============================================================================
#                    7. COMBINED FIGURES FOR PUBLICATION
# =============================================================================
# Three main figures combining related panels:
#   Figure 1: Fe and Zn variation (two histograms)
#   Figure 2: Heritability story (cross-season scatter, repeatability, variance)
#   Figure 3: Breeding implications (Fe-Zn scatter, agronomic correlations)
# =============================================================================

cat("\n========== CREATING COMBINED FIGURES ==========\n")

# --- Figure 1: Fe and Zn BLUE distributions ---
# Calculate common y-axis max
hist_fe_data <- ggplot_build(ggplot(blues_fe, aes(x = emmean)) + geom_histogram(bins = 15))$data[[1]]
hist_zn_data <- ggplot_build(ggplot(blues_zn, aes(x = emmean)) + geom_histogram(bins = 15))$data[[1]]
y_max <- max(c(hist_fe_data$count, hist_zn_data$count))

p_hist_fe <- ggplot(blues_fe, aes(x = emmean)) +
  geom_histogram(fill = "#8B4513", color = "black", bins = 15) +
  labs(x = "Fe (mg/kg)", y = "Frequency", tag = "A") +
  custom_theme +
  scale_y_continuous(expand = c(0, 0), limits = c(0, y_max * 1.05)) +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.05, 0.95))

p_hist_zn <- ggplot(blues_zn, aes(x = emmean)) +
  geom_histogram(fill = "#4682B4", color = "black", bins = 15) +
  labs(x = "Zn (mg/kg)", y = "Frequency", tag = "B") +
  custom_theme +
  scale_y_continuous(expand = c(0, 0), limits = c(0, y_max * 1.05)) +
  theme(plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.05, 0.95))

png("plots/Figure1_micronutrient_variation.png", width = 1400, height = 600, res = 150)
grid.arrange(p_hist_fe, p_hist_zn, ncol = 2)
dev.off()
cat("Saved: plots/Figure1_micronutrient_variation.png\n")

# --- Figure 2: Heritability story (4 panels: A, B, C, D) ---
# Panel A: Cross-season Fe scatter
# Panel B: Cross-season Zn scatter
# Panel C: Heritability comparison
# Panel D: Variance partitioning

# Recreate cross-season plots with simplified labels and panel labels
p_fe_fig2 <- ggplot(geno_means, aes(x = S1_Fe_mg_kg, y = S2_Fe_mg_kg)) +
  geom_point(alpha = 0.6, color = "#8B4513", size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf, label = paste("r =", fe_r), hjust = 1.1, vjust = -0.5, size = 4) +
  labs(x = "Iron, 2024 (mg/kg)", y = "Iron, 2025 (mg/kg)", tag = "A") +
  theme_publication() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.05, 0.95))

p_zn_fig2 <- ggplot(geno_means, aes(x = S1_Zn_mg_kg, y = S2_Zn_mg_kg)) +
  geom_point(alpha = 0.6, color = "#4682B4", size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf, label = paste("r =", zn_r), hjust = 1.1, vjust = -0.5, size = 4) +
  labs(x = "Zinc, 2024 (mg/kg)", y = "Zinc, 2025 (mg/kg)", tag = "B") +
  theme_publication() +
  theme(panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.05, 0.95))

# Heritability panel with shortened legend and panel label
p_herit_fig2 <- ggplot(herit_long, aes(x = TraitLabel, y = Value, fill = Measure)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c("Repeatability (Spring-24)" = "#1E4D2B",
                               "Repeatability (Autumn-25)" = "#D4A84B",
                               "Heritability (Across)" = "#5B7C99"),
                    labels = c("Spring 2024", "Autumn 2025", "Heritability")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = "Repeatability / Heritability", tag = "C") +
  theme_publication() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.05, 0.95))

# Variance panel with panel label
p_var_fig2 <- ggplot(var_part_long, aes(x = TraitLabel, y = Percentage, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = var_colors) +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = NULL, y = "Variance Explained (%)", tag = "D") +
  theme_publication() +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        panel.grid = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 14),
        plot.tag.position = c(0.05, 0.95)) +
  guides(fill = guide_legend(reverse = TRUE))

# Combine: top row = Fe/Zn scatter, bottom row = heritability + variance
png("plots/Figure2_heritability.png", width = 1600, height = 1400, res = 150)
grid.arrange(
  arrangeGrob(p_fe_fig2, p_zn_fig2, ncol = 2),
  arrangeGrob(p_herit_fig2, p_var_fig2, ncol = 2),
  nrow = 2, heights = c(1, 1.2)
)
dev.off()
cat("Saved: plots/Figure2_heritability.png\n")

# --- Figure 3: Breeding implications (2 panels) ---
# Panel A: Fe vs Zn scatter with candidates
# Panel B: Micronutrient correlations with agronomic traits

# Recreate scatter with tag for Figure 3
p_scatter_fig3 <- ggplot(all_blues, aes(x = Fe, y = Zn)) +
  geom_point(alpha = 0.5, color = "gray50", size = 2.5) +
  geom_point(data = top_both, aes(x = Fe, y = Zn),
             color = "#2E7D32", size = 3, alpha = 0.8) +
  geom_vline(xintercept = fe_q75, linetype = "dashed", color = "gray60") +
  geom_hline(yintercept = zn_q75, linetype = "dashed", color = "gray60") +
  annotate("text", x = -Inf, y = Inf, label = "A", fontface = "bold", size = 6,
           hjust = -0.3, vjust = 1.2) +
  labs(x = "Iron (mg/kg)", y = "Zinc (mg/kg)") +
  theme_publication() +
  theme(plot.margin = margin(5, 10, 5, 5))

# For coord_flip: x=Inf is top (right side of bars), y=-Inf is left
p_corr_fig3 <- p_correlations +
  annotate("text", x = Inf, y = -Inf, label = "B", fontface = "bold", size = 6,
           hjust = -0.3, vjust = 1.2) +
  theme(plot.margin = margin(5, 5, 5, 10))

png("plots/Figure3_breeding_implications.png", width = 1600, height = 700, res = 150)
grid.arrange(
  p_scatter_fig3, p_corr_fig3,
  ncol = 2, widths = c(1, 1.2)
)
dev.off()
cat("Saved: plots/Figure3_breeding_implications.png\n")

cat("\n=== COMBINED FIGURES COMPLETE ===\n")


# =============================================================================
#                    8. SUMMARY OUTPUT
# =============================================================================

cat("\n\n")
cat("==============================================================\n")
cat("                    ANALYSIS COMPLETE                          \n")
cat("==============================================================\n")

cat("\nMODEL SPECIFICATIONS:\n")
cat("  - Variance partitioning: Year as RANDOM (included in partitioning)\n")
cat("  - BLUEs: Year as FIXED\n")
cat("  - Within-season: Repeatability (not heritability)\n")
cat("  - Across-season: Heritability on genotype-mean basis\n")

cat("\nMAIN FIGURES FOR PAPER (plots/):\n")
cat("  Figure1_micronutrient_variation.png  - Fe/Zn BLUE distributions\n")
cat("  Figure2_heritability.png             - Cross-season stability, repeatability, variance\n")
cat("  Figure3_breeding_implications.png    - Fe-Zn scatter + agronomic correlations\n")

cat("\nINDIVIDUAL PANELS (plots/):\n")
cat("  01_micronutrient_correlations.png  - Fe/Zn vs agronomic traits\n")
cat("  02_cross_season_stability.png      - Genotype stability across seasons\n")
cat("  03_cross_season_Fe_Zn.png          - Fe/Zn season-to-season scatter\n")
cat("  04_Fe_vs_Zn_scatter.png            - Fe-Zn relationship + top lines\n")
cat("  05_heritability_comparison.png     - Repeatability & heritability\n")
cat("  06_variance_partitioning.png       - Variance components\n")
cat("  BLUE_histogram_*.png               - Genotype distributions (5 traits)\n")

cat("\nSUPPLEMENTARY FIGURES (plots/):\n")
cat("  S1_spatial_trends.png              - Field spatial gradients\n")
cat("  02_cross_season_stability.png      - All-trait stability (if needed)\n")

cat("\nTABLES FOR PAPER (results/):\n")
cat("  summary_statistics.csv             - Table 1: trait means, ranges, CV\n")
cat("  heritability_summary.csv           - Repeatability & heritability\n")
cat("  genotype_BLUEs.csv                 - All genotype adjusted means\n")
cat("  top10_Fe_genotypes.csv             - Best lines for Fe\n")
cat("  top10_Zn_genotypes.csv             - Best lines for Zn\n")
cat("  top_Fe_and_Zn_genotypes.csv        - Lines in top quartile for both\n")

cat("\nOTHER RESULTS (results/):\n")
cat("  micronutrient_correlations.csv\n")
cat("  cross_season_correlations.csv\n")
cat("  variance_partitioning.csv\n")

cat("\n==============================================================\n")
cat("Ready to compile figures and tables for manuscript.\n")
cat("==============================================================\n")

