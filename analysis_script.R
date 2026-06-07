# =============================================================================
#                    MAIZE PHENOTYPIC DATA ANALYSIS
#                    Correlation, Spatial Trends, Heritability
#                    Publication-Quality Plots
# =============================================================================

# Set library path
.libPaths(c("~/R_libs", .libPaths()))

# Load required libraries
library(tidyverse)
library(lme4)
library(corrplot)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(viridis)

# Set working directory
setwd("/work/schnablelab/waqarali/Hira_Aslam")

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

      # Grid lines - subtle
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
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

season1 <- read.csv("Season_1.csv")
season2 <- read.csv("Season_2.csv")
whole <- read.csv("Whole.csv")

# Clean column names
clean_names <- function(df) {
  names(df) <- gsub("\\.", "", names(df))
  names(df) <- gsub("\\s+", "", names(df))
  names(df) <- gsub("/", "_", names(df))
  names(df) <- gsub("\\(", "_", names(df))
  names(df) <- gsub("\\)", "", names(df))
  return(df)
}

season1 <- clean_names(season1)
season2 <- clean_names(season2)
whole <- clean_names(whole)

# Get trait columns
traits <- names(whole)[5:15]
cat("Traits to analyze:\n")
print(traits)

# Trait labels for better display
trait_labels <- c(
  "DTAdays" = "Days to Anthesis",
  "DTSdays" = "Days to Silking",
  "ASIdays" = "Anthesis-Silking Interval",
  "EHcm" = "Ear Height (cm)",
  "PHcm" = "Plant Height (cm)",
  "NORC" = "Rows per Cob",
  "NOKR" = "Kernels per Row",
  "HGWg" = "100 Grain Weight (g)",
  "DTMdays" = "Days to Maturity",
  "Femgkg" = "Iron (mg/kg)",
  "Znmgkg" = "Zinc (mg/kg)"
)

# =============================================================================
#                    2. CORRELATION ANALYSIS
# =============================================================================

cat("\n========== CORRELATION ANALYSIS ==========\n")

# --- 2a. Season 1 Correlation ---
cor_season1 <- cor(season1[, traits], use = "complete.obs")
colnames(cor_season1) <- trait_labels[colnames(cor_season1)]
rownames(cor_season1) <- trait_labels[rownames(cor_season1)]

png("plots/01_correlation_season1.png", width = 1200, height = 1000, res = 150)
par(mar = c(1, 1, 3, 1))
corrplot(cor_season1, method = "color", type = "lower",
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.6,
         col = colorRampPalette(c("#2E86AB", "white", "#E94F37"))(100),
         mar = c(0, 0, 1, 0),
         cl.cex = 0.8)
dev.off()

# --- 2b. Season 2 Correlation ---
cor_season2 <- cor(season2[, traits], use = "complete.obs")
colnames(cor_season2) <- trait_labels[colnames(cor_season2)]
rownames(cor_season2) <- trait_labels[rownames(cor_season2)]

png("plots/02_correlation_season2.png", width = 1200, height = 1000, res = 150)
par(mar = c(1, 1, 3, 1))
corrplot(cor_season2, method = "color", type = "lower",
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.6,
         col = colorRampPalette(c("#2E86AB", "white", "#E94F37"))(100),
         mar = c(0, 0, 1, 0),
         cl.cex = 0.8)
dev.off()

# --- 2c. Combined Correlation ---
cor_whole <- cor(whole[, traits], use = "complete.obs")
colnames(cor_whole) <- trait_labels[colnames(cor_whole)]
rownames(cor_whole) <- trait_labels[rownames(cor_whole)]

png("plots/03_correlation_combined.png", width = 1200, height = 1000, res = 150)
par(mar = c(1, 1, 3, 1))
corrplot(cor_whole, method = "color", type = "lower",
         tl.col = "black", tl.srt = 45, tl.cex = 0.8,
         addCoef.col = "black", number.cex = 0.6,
         col = colorRampPalette(c("#2E86AB", "white", "#E94F37"))(100),
         mar = c(0, 0, 1, 0),
         cl.cex = 0.8)
dev.off()

# --- 2d. Between Season Correlation ---
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

cat("\nBetween-Season Correlations:\n")
print(between_season_cor)
write.csv(between_season_cor, "results/between_season_correlations.csv", row.names = FALSE)

# Plot between-season correlations
png("plots/04_between_season_correlations.png", width = 1000, height = 700, res = 150)
ggplot(between_season_cor, aes(x = reorder(TraitLabel, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", fill = "#2E86AB", width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  coord_flip() +
  scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Correlation Coefficient (r)") +
  theme_publication() +
  theme(panel.grid.major.y = element_blank())
dev.off()

# --- 2e. Scatter plots between seasons ---
png("plots/05_between_season_scatter.png", width = 1600, height = 1400, res = 150)
plots_list <- list()
for (i in seq_along(traits)) {
  t <- traits[i]
  s1_col <- paste0("S1_", t)
  s2_col <- paste0("S2_", t)
  r_val <- round(cor(geno_means[[s1_col]], geno_means[[s2_col]], use = "complete.obs"), 2)

  p <- ggplot(geno_means, aes_string(x = s1_col, y = s2_col)) +
    geom_point(alpha = 0.6, color = "#2E86AB", size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "#E94F37", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, color = "gray50", linetype = "dashed") +
    labs(title = trait_labels[t],
         subtitle = paste("r =", r_val),
         x = "Season 1", y = "Season 2") +
    theme_publication(base_size = 10) +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 9))
  plots_list[[i]] <- p
}
do.call(grid.arrange, c(plots_list, ncol = 4))
dev.off()

# =============================================================================
#                    3. SPATIAL TREND ANALYSIS
# =============================================================================

cat("\n========== SPATIAL TREND ANALYSIS ==========\n")

# --- 3a. Spatial trends - Season 1 ---
png("plots/06_spatial_trends_season1.png", width = 1600, height = 1400, res = 150)
plots_s1 <- list()
for (i in seq_along(traits)) {
  t <- traits[i]
  p <- ggplot(season1, aes_string(x = "Row", y = t)) +
    geom_point(alpha = 0.4, color = "#2E86AB", size = 1.5) +
    geom_smooth(method = "loess", se = TRUE, color = "#E94F37", fill = "#E94F37", alpha = 0.2, linewidth = 1) +
    labs(title = trait_labels[t], x = "Row Number", y = NULL) +
    theme_publication(base_size = 10) +
    theme(plot.title = element_text(size = 10))
  plots_s1[[i]] <- p
}
grid.arrange(grobs = plots_s1, ncol = 4)
dev.off()

# --- 3b. Spatial trends - Season 2 ---
png("plots/07_spatial_trends_season2.png", width = 1600, height = 1400, res = 150)
plots_s2 <- list()
for (i in seq_along(traits)) {
  t <- traits[i]
  p <- ggplot(season2, aes_string(x = "Row", y = t)) +
    geom_point(alpha = 0.4, color = "#E94F37", size = 1.5) +
    geom_smooth(method = "loess", se = TRUE, color = "#2E86AB", fill = "#2E86AB", alpha = 0.2, linewidth = 1) +
    labs(title = trait_labels[t], x = "Row Number", y = NULL) +
    theme_publication(base_size = 10) +
    theme(plot.title = element_text(size = 10))
  plots_s2[[i]] <- p
}
grid.arrange(grobs = plots_s2, ncol = 4)
dev.off()

# --- 3c. Combined spatial plot ---
png("plots/08_spatial_trends_combined.png", width = 1600, height = 1400, res = 150)
plots_comb <- list()
for (i in seq_along(traits)) {
  t <- traits[i]
  p <- ggplot(whole, aes_string(x = "Row", y = t, color = "Year")) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
    scale_color_manual(values = season_colors) +
    labs(title = trait_labels[t], x = "Row Number", y = NULL) +
    theme_publication(base_size = 10) +
    theme(plot.title = element_text(size = 10),
          legend.position = "none")
  plots_comb[[i]] <- p
}

# Add legend to one plot
plots_comb[[1]] <- plots_comb[[1]] + theme(legend.position = c(0.8, 0.2))

grid.arrange(grobs = plots_comb, ncol = 4)
dev.off()

# =============================================================================
#                    4. HERITABILITY ANALYSIS (using lme4)
# =============================================================================

cat("\n========== HERITABILITY ANALYSIS ==========\n")

# Heritability formula: H = σ²_G / (σ²_G + (1/n) * σ²_R)
# n = 2 replications

calc_heritability <- function(data, trait, n_reps = 2) {
  formula <- as.formula(paste(trait, "~ (1|InbCode)"))

  tryCatch({
    model <- lmer(formula, data = data, REML = TRUE)
    var_comp <- as.data.frame(VarCorr(model))
    var_geno <- var_comp$vcov[var_comp$grp == "InbCode"]
    var_resid <- var_comp$vcov[var_comp$grp == "Residual"]

    heritability <- var_geno / (var_geno + (1/n_reps) * var_resid)

    return(list(trait = trait, var_geno = var_geno, var_resid = var_resid, heritability = heritability))
  }, error = function(e) {
    return(list(trait = trait, var_geno = NA, var_resid = NA, heritability = NA))
  })
}

# Season 1 heritability
herit_s1 <- lapply(traits, function(t) calc_heritability(season1, t, n_reps = 2))
herit_s1_df <- do.call(rbind, lapply(herit_s1, function(x) {
  data.frame(Trait = x$trait, Var_Geno = x$var_geno, Var_Resid = x$var_resid, Heritability = x$heritability)
}))

# Season 2 heritability
herit_s2 <- lapply(traits, function(t) calc_heritability(season2, t, n_reps = 2))
herit_s2_df <- do.call(rbind, lapply(herit_s2, function(x) {
  data.frame(Trait = x$trait, Var_Geno = x$var_geno, Var_Resid = x$var_resid, Heritability = x$heritability)
}))

# Across seasons heritability
calc_heritability_across <- function(data, trait, n_reps = 2, n_years = 2) {
  formula <- as.formula(paste(trait, "~ Year + (1|InbCode)"))

  tryCatch({
    model <- lmer(formula, data = data, REML = TRUE)
    var_comp <- as.data.frame(VarCorr(model))
    var_geno <- var_comp$vcov[var_comp$grp == "InbCode"]
    var_resid <- var_comp$vcov[var_comp$grp == "Residual"]

    heritability <- var_geno / (var_geno + var_resid / (n_reps * n_years))

    return(list(trait = trait, var_geno = var_geno, var_resid = var_resid, heritability = heritability))
  }, error = function(e) {
    return(list(trait = trait, var_geno = NA, var_resid = NA, heritability = NA))
  })
}

herit_across <- lapply(traits, function(t) calc_heritability_across(whole, t))
herit_across_df <- do.call(rbind, lapply(herit_across, function(x) {
  data.frame(Trait = x$trait, Var_Geno = x$var_geno, Var_Resid = x$var_resid, Heritability = x$heritability)
}))

# Combined summary
herit_summary <- data.frame(
  Trait = herit_s1_df$Trait,
  TraitLabel = trait_labels[herit_s1_df$Trait],
  H2_Season1 = round(herit_s1_df$Heritability, 3),
  H2_Season2 = round(herit_s2_df$Heritability, 3),
  H2_Across = round(herit_across_df$Heritability, 3)
)

cat("\n=== HERITABILITY SUMMARY ===\n")
print(herit_summary)
write.csv(herit_summary, "results/heritability_summary.csv", row.names = FALSE)

# Plot heritability - VERTICAL BARS
herit_long <- herit_summary %>%
  pivot_longer(cols = starts_with("H2"), names_to = "Season", values_to = "Heritability") %>%
  mutate(Season = recode(Season,
                         "H2_Season1" = "Season 1",
                         "H2_Season2" = "Season 2",
                         "H2_Across" = "Across Seasons")) %>%
  mutate(Season = factor(Season, levels = c("Season 1", "Season 2", "Across Seasons")))

# Order traits by mean heritability (low to high)
trait_order <- herit_long %>%
  group_by(TraitLabel) %>%
  summarise(mean_h2 = mean(Heritability, na.rm = TRUE)) %>%
  arrange(mean_h2) %>%
  pull(TraitLabel)

herit_long$TraitLabel <- factor(herit_long$TraitLabel, levels = trait_order)

png("plots/09_heritability_comparison.png", width = 1400, height = 800, res = 150)
ggplot(herit_long, aes(x = TraitLabel, y = Heritability, fill = Season)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c("Season 1" = "#1E4D2B", "Season 2" = "#D4A84B", "Across Seasons" = "#5B7C99")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = NULL, y = expression(Heritability~(H^2))) +
  theme_publication() +
  theme(legend.position = "top",
        legend.justification = "center",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank())
dev.off()

# =============================================================================
#                    5. VARIANCE PARTITIONING
# =============================================================================

cat("\n========== VARIANCE PARTITIONING ==========\n")

variance_partition <- function(data, trait) {
  data$InbCode <- as.factor(data$InbCode)
  data$Year <- as.factor(data$Year)
  data$Row <- as.factor(data$Row)

  formula <- as.formula(paste(trait, "~ (1|InbCode) + (1|Year) + (1|Row) + (1|InbCode:Year)"))

  tryCatch({
    model <- lmer(formula, data = data, REML = TRUE,
                  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

    var_comp <- as.data.frame(VarCorr(model))

    var_geno <- var_comp$vcov[var_comp$grp == "InbCode"]
    var_year <- var_comp$vcov[var_comp$grp == "Year"]
    var_row <- var_comp$vcov[var_comp$grp == "Row"]
    var_gxe <- var_comp$vcov[var_comp$grp == "InbCode:Year"]
    var_resid <- var_comp$vcov[var_comp$grp == "Residual"]

    # Handle NAs
    var_geno <- ifelse(is.null(var_geno) || length(var_geno) == 0, 0, var_geno)
    var_year <- ifelse(is.null(var_year) || length(var_year) == 0, 0, var_year)
    var_row <- ifelse(is.null(var_row) || length(var_row) == 0, 0, var_row)
    var_gxe <- ifelse(is.null(var_gxe) || length(var_gxe) == 0, 0, var_gxe)

    total_var <- var_geno + var_year + var_row + var_gxe + var_resid

    return(data.frame(
      Trait = trait,
      Pct_Genotype = round(100 * var_geno / total_var, 2),
      Pct_Year = round(100 * var_year / total_var, 2),
      Pct_Row = round(100 * var_row / total_var, 2),
      Pct_GxE = round(100 * var_gxe / total_var, 2),
      Pct_Residual = round(100 * var_resid / total_var, 2)
    ))
  }, error = function(e) {
    cat("Error for", trait, ":", e$message, "\n")
    return(data.frame(Trait = trait, Pct_Genotype = NA, Pct_Year = NA,
                      Pct_Row = NA, Pct_GxE = NA, Pct_Residual = NA))
  })
}

var_part_results <- do.call(rbind, lapply(traits, function(t) variance_partition(whole, t)))
var_part_results$TraitLabel <- trait_labels[var_part_results$Trait]

cat("\n=== VARIANCE PARTITIONING (%) ===\n")
print(var_part_results)
write.csv(var_part_results, "results/variance_partitioning.csv", row.names = FALSE)

# Plot variance partitioning - VERTICAL BARS with green color scheme
var_part_long <- var_part_results %>%
  select(Trait, TraitLabel, starts_with("Pct")) %>%
  pivot_longer(cols = starts_with("Pct"), names_to = "Component", values_to = "Percentage") %>%
  mutate(Component = recode(Component,
                            "Pct_Genotype" = "Genotype",
                            "Pct_Year" = "Year",
                            "Pct_Row" = "Row",
                            "Pct_GxE" = "G x E",
                            "Pct_Residual" = "Residual")) %>%
  mutate(Component = factor(Component, levels = c("Residual", "G x E", "Row", "Year", "Genotype")))

# Premium color scheme: Distinct, muted, sophisticated colors
var_colors <- c(
  "Genotype" = "#1E4D2B",    # Deep forest green
  "Year" = "#D4A84B",        # Muted gold
  "Row" = "#5B7C99",         # Slate blue
  "G x E" = "#C67D5E",       # Terracotta
  "Residual" = "#8C8C8C"     # Sophisticated grey
)

png("plots/10_variance_partitioning.png", width = 1400, height = 800, res = 150)
ggplot(var_part_long, aes(x = TraitLabel, y = Percentage, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = var_colors) +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = NULL, y = "Variance Explained (%)") +
  theme_publication() +
  theme(legend.position = "top",
        legend.justification = "center",
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major.x = element_blank()) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE))
dev.off()

# =============================================================================
#                    6. SUMMARY OUTPUT
# =============================================================================

cat("\n\n")
cat("==============================================================\n")
cat("                    ANALYSIS COMPLETE                          \n")
cat("==============================================================\n")

cat("\nPLOTS SAVED IN: plots/\n")
cat("  01_correlation_season1.png\n")
cat("  02_correlation_season2.png\n")
cat("  03_correlation_combined.png\n")
cat("  04_between_season_correlations.png\n")
cat("  05_between_season_scatter.png\n")
cat("  06_spatial_trends_season1.png\n")
cat("  07_spatial_trends_season2.png\n")
cat("  08_spatial_trends_combined.png\n")
cat("  09_heritability_comparison.png\n")
cat("  10_variance_partitioning.png\n")

cat("\nRESULTS SAVED IN: results/\n")
cat("  between_season_correlations.csv\n")
cat("  heritability_summary.csv\n")
cat("  variance_partitioning.csv\n")

cat("\n==============================================================\n")
