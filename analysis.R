# Resource limitation shapes the evolution of sexual dimorphism and condition dependence
# in Drosophila melanogaster
#
# Reproducible analysis for Kashyap et al.
# R version used for the revised manuscript: 4.4.1
#
# Run this script from the repository root.
# It reproduces the principal analyses, model comparisons, effect sizes,
# diagnostics, and main/supplementary figures described in the manuscript.

# Packages
required_packages <- c(
  "car", "DHARMa", "dplyr", "emmeans", "ggplot2", "lme4",
  "patchwork", "reshape2", "tibble", "tidyr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running the analysis: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(car)
library(DHARMa)
library(dplyr)
library(emmeans)
library(ggplot2)
library(lme4)
library(patchwork)
library(reshape2)
library(tibble)
library(tidyr)

options(contrasts = c("contr.sum", "contr.poly"))

dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# Read and clean data

dat <- read.csv("data/cdsd_dat.csv", stringsAsFactors = FALSE) |>
  rename(Treatment = Treatement) |>
  mutate(
    Block = factor(Block),
    Selection = factor(Selection, levels = c("CU", "MB")),
    Treatment = factor(Treatment, levels = c("HD", "ID", "LD")),
    Sex = factor(Sex, levels = c("Female", "Male"))
  ) |>
  filter(ID != "CImF10")


pca_traits <- c("Femur", "Tibia", "Thorax", "LW", "RW")

#Composite body size (PCA) 

body_pca <- prcomp(dat[, pca_traits], center = TRUE, scale. = TRUE)
analysis_dat <- dat |>
  mutate(PC1 = body_pca$x[, 1])

pca_loadings <- as.data.frame(body_pca$rotation) |>
  rownames_to_column("Trait")

pca_variance <- data.frame(
  PC = paste0("PC", seq_along(body_pca$sdev)),
  Variance_percent = 100 * body_pca$sdev^2 / sum(body_pca$sdev^2)
)

write.csv(pca_loadings, "outputs/tables/pca_loadings.csv", row.names = FALSE)
write.csv(pca_variance, "outputs/tables/pca_variance_explained.csv", row.names = FALSE)

# Figure S1: correlation matrix including PC1.
cor_dat <- analysis_dat |>
  select(all_of(pca_traits), PC1)

cor_long <- reshape2::melt(cor(cor_dat, use = "complete.obs"))
names(cor_long) <- c("Trait_1", "Trait_2", "Correlation")

fig_s1 <- ggplot(cor_long, aes(Trait_1, Trait_2, fill = Correlation)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0.5, limits = c(0, 1), name = "Correlation"
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("outputs/figures/Figure_S1_correlation_heatmap.pdf", fig_s1,
       width = 7, height = 6)

# Composite body size mixed model 

# AIC comparisons must be made with maximum likelihood (REML = FALSE).
model_full_ml <- lmer(
  PC1 ~ Treatment * Sex * Selection + (1 | Block),
  data = analysis_dat, REML = FALSE
)
model_treatment_selection_ml <- lmer(
  PC1 ~ Treatment * Selection + Sex + (1 | Block),
  data = analysis_dat, REML = FALSE
)
model_treatment_sex_ml <- lmer(
  PC1 ~ Treatment * Sex + Selection + (1 | Block),
  data = analysis_dat, REML = FALSE
)
model_additive_ml <- lmer(
  PC1 ~ Treatment + Sex + Selection + (1 | Block),
  data = analysis_dat, REML = FALSE
)

model_comparison <- data.frame(
  Model = c(
    "Full factorial model",
    "Treatment * Selection + Sex",
    "Sex * Treatment + Selection",
    "Sex + Treatment + Selection"
  ),
  AIC = c(
    AIC(model_full_ml),
    AIC(model_treatment_selection_ml),
    AIC(model_treatment_sex_ml),
    AIC(model_additive_ml)
  ),
  Log_likelihood = c(
    as.numeric(logLik(model_full_ml)),
    as.numeric(logLik(model_treatment_selection_ml)),
    as.numeric(logLik(model_treatment_sex_ml)),
    as.numeric(logLik(model_additive_ml))
  )
)
write.csv(model_comparison, "outputs/tables/Table_S2_model_comparison.csv",
          row.names = FALSE)

# The inferential model reported in the manuscript is fitted by REML.
model_pc1 <- lmer(
  PC1 ~ Treatment * Sex * Selection + (1 | Block),
  data = analysis_dat, REML = TRUE
)

pc1_anova <- car::Anova(model_pc1, type = 3) |>
  as.data.frame() |>
  rownames_to_column("Effect")
write.csv(pc1_anova, "outputs/tables/Table_1_PC1_typeIII_ANOVA.csv",
          row.names = FALSE)

# Model diagnostics and outlier sensitivity 

png("outputs/figures/PC1_model_diagnostics.png", width = 1800, height = 600,
    res = 180)
par(mfrow = c(1, 3))
hist(resid(model_pc1), main = "Residuals", xlab = "Residual")
qqnorm(resid(model_pc1)); qqline(resid(model_pc1))
plot(fitted(model_pc1), resid(model_pc1),
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2)
dev.off()

set.seed(1)
dharma_pc1 <- DHARMa::simulateResiduals(model_pc1, plot = FALSE)
png("outputs/figures/PC1_DHARMa_diagnostics.png", width = 1200, height = 900,
    res = 150)
plot(dharma_pc1)
dev.off()

standardized_residual <- as.numeric(scale(resid(model_pc1)))
outlier_rows <- which(abs(standardized_residual) > 3)
outliers <- analysis_dat[outlier_rows, c("ID", "Block", "Selection",
                                         "Treatment", "Sex", "PC1")]
outliers$standardized_residual <- standardized_residual[outlier_rows]
write.csv(outliers, "outputs/tables/PC1_residual_outliers.csv", row.names = FALSE)

# Sensitivity check described in the manuscript: refit after removing residual
# tail observations (> 3 SD). Primary results above always use the full dataset.
if (length(outlier_rows) > 0) {
  model_pc1_no_outliers <- lmer(
    PC1 ~ Treatment * Sex * Selection + (1 | Block),
    data = analysis_dat[-outlier_rows, ], REML = TRUE
  )
  no_outlier_anova <- car::Anova(model_pc1_no_outliers, type = 3) |>
    as.data.frame() |>
    rownames_to_column("Effect")
  write.csv(
    no_outlier_anova,
    "outputs/tables/PC1_typeIII_ANOVA_outlier_sensitivity.csv",
    row.names = FALSE
  )
}

# Pairwise comparisons and sexual dimorphism effect sizes 

emm_sex_selection <- emmeans(
  model_pc1, ~ Sex * Selection | Treatment
)
table_s3 <- as.data.frame(pairs(emm_sex_selection, adjust = "tukey"))
write.csv(table_s3, "outputs/tables/Table_S3_pairwise_contrasts_PC1.csv",
          row.names = FALSE)

emm_sex <- emmeans(model_pc1, ~ Sex | Treatment * Selection)
ssd_pc1 <- as.data.frame(
  eff_size(emm_sex, sigma = sigma(model_pc1), edf = df.residual(model_pc1))
)
write.csv(ssd_pc1, "outputs/tables/Table_S4_sexual_dimorphism_effect_sizes.csv",
          row.names = FALSE)

fig_s2 <- ggplot(
  ssd_pc1,
  aes(Selection, effect.size, shape = Selection)
) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.25, linewidth = 0.4
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~ Treatment) +
  scale_x_discrete(labels = c(CU = "Crowding evolved", MB = "Control")) +
  labs(
    x = "Experimental populations",
    y = "Sexual dimorphism in composite body size (Cohen's d for Female - Male contrast)"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none")

ggsave("outputs/figures/Figure_S2_sexual_dimorphism_effect_sizes.pdf",
       fig_s2, width = 8, height = 4.5)

#  Density sensitivity in females and males 

emm_density <- emmeans(model_pc1, ~ Treatment | Sex * Selection)
density_effects <- as.data.frame(
  eff_size(emm_density, sigma = sigma(model_pc1), edf = df.residual(model_pc1))
)
density_effects$contrast <- factor(
  density_effects$contrast,
  levels = c("HD - ID", "HD - LD", "ID - LD")
)
write.csv(density_effects, "outputs/tables/PC1_density_effect_sizes.csv",
          row.names = FALSE)

# Compare the Cohen's d estimates between selection regimes for each sex and
# density contrast. This reproduces the z-tests reported in Table S5.
density_wide <- density_effects |>
  select(contrast, Sex, Selection, effect.size, SE) |>
  pivot_wider(
    names_from = Selection,
    values_from = c(effect.size, SE)
  ) |>
  mutate(
    difference = effect.size_CU - effect.size_MB,
    SE_pooled = sqrt(SE_CU^2 + SE_MB^2),
    z = difference / SE_pooled,
    p_value = 2 * pnorm(abs(z), lower.tail = FALSE)
  )
write.csv(density_wide, "outputs/tables/Table_S5_density_effect_comparisons.csv",
          row.names = FALSE)

# Figure 1: model-estimated composite body size across larval density.
pc1_means <- as.data.frame(
  emmeans(model_pc1, ~ Selection | Sex * Treatment)
)
names(pc1_means)[names(pc1_means) == "lower.CL"] <- "lcl"
names(pc1_means)[names(pc1_means) == "upper.CL"] <- "ucl"

selection_labels <- c(MB = "Control", CU = "Crowding evolved")
selection_palette <- c(CU = "#0072B2", MB = "#D55E00")

fig1 <- ggplot(
  pc1_means,
  aes(Treatment, emmean, colour = Selection, shape = Selection,
      group = Selection)
) +
  geom_line(linewidth = 0.6, position = position_dodge(width = 0.2)) +
  geom_errorbar(
    aes(ymin = lcl, ymax = ucl),
    width = 0.12, position = position_dodge(width = 0.2)
  ) +
  geom_point(size = 2.5, position = position_dodge(width = 0.2)) +
  facet_wrap(~ Sex) +
  scale_colour_manual(values = selection_palette, labels = selection_labels) +
  scale_shape_manual(values = c(CU = 16, MB = 17), labels = selection_labels) +
  labs(
    x = "Larval density",
    y = "Composite body size (PC1)",
    colour = NULL, shape = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

ggsave("outputs/figures/Figure_1_composite_body_size.pdf", fig1,
       width = 7, height = 4.5)

# Figure 2: standardized within-sex density contrasts.
fig2 <- ggplot(
  density_effects,
  aes(contrast, effect.size, colour = Selection, shape = Selection)
) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    position = position_dodge(0.4), width = 0.5, linewidth = 0.4
  ) +
  geom_point(position = position_dodge(0.4), size = 2.3) +
  facet_wrap(~ Sex) +
  scale_colour_manual(values = selection_palette, labels = selection_labels) +
  scale_shape_manual(values = c(CU = 16, MB = 17), labels = selection_labels) +
  labs(x = "Density contrast", y = "Cohen's d",
       colour = "Selection", shape = "Selection") +
  theme_bw(base_size = 12)

ggsave("outputs/figures/Figure_2_density_effect_sizes.pdf", fig2,
       width = 7, height = 4.5)

#  Trait-wise body-size correction 

# The manuscript corrects each focal trait against a PCA body-size axis that
# excludes that focal trait. The six plotted variables include left and right
# wing measurements separately and their average.
trait_vars <- c("Femur", "LW", "RW", "Thorax", "Tibia", "Wings")

body_size_correct <- function(data, traits) {
  corrected <- data

  for (trait in traits) {
    other_traits <- setdiff(traits, trait)
    pca <- prcomp(data[, other_traits], center = TRUE, scale. = TRUE)
    size_axis <- pca$x[, 1]
    corrected[[paste0(trait, "_resid")]] <-
      resid(lm(data[[trait]] ~ size_axis))
  }

  corrected
}

trait_dat <- body_size_correct(analysis_dat, trait_vars)

fit_trait_models <- function(data, traits) {
  setNames(
    lapply(traits, function(trait) {
      response <- paste0(trait, "_resid")
      lmer(
        as.formula(
          paste(response, "~ Treatment * Sex * Selection + (1 | Block)")
        ),
        data = data, REML = TRUE
      )
    }),
    traits
  )
}

trait_models <- fit_trait_models(trait_dat, trait_vars)

table_s6 <- bind_rows(lapply(names(trait_models), function(trait) {
  car::Anova(trait_models[[trait]], type = 3) |>
    as.data.frame() |>
    rownames_to_column("Effect") |>
    mutate(Trait = trait, .before = 1)
}))
write.csv(table_s6, "outputs/tables/Table_S6_traitwise_typeIII_ANOVA.csv",
          row.names = FALSE)

# Figure S3: estimated means for body-size-corrected traits.
trait_means <- bind_rows(lapply(names(trait_models), function(trait) {
  as.data.frame(
    emmeans(trait_models[[trait]], ~ Sex * Selection | Treatment)
  ) |>
    mutate(Trait = trait)
}))

fig_s3 <- ggplot(
  trait_means,
  aes(Selection, emmean, colour = Sex, shape = Sex, group = Sex)
) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    position = position_dodge(0.45), width = 0.25, linewidth = 0.35
  ) +
  geom_point(position = position_dodge(0.45), size = 1.8) +
  facet_grid(Trait ~ Treatment, scales = "free_y") +
  scale_x_discrete(labels = c(CU = "Crowding\nevolved", MB = "Control")) +
  labs(x = NULL, y = "Size-corrected residual", colour = NULL, shape = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

ggsave("outputs/figures/Figure_S3_traitwise_size_corrected.pdf",
       fig_s3, width = 9, height = 10)

# Wing-length based correction sensitivity analysis

wing_correct_traits <- c("Femur", "Thorax", "Tibia")
wing_dat <- analysis_dat

for (trait in wing_correct_traits) {
  wing_dat[[paste0(trait, "_resid")]] <-
    resid(lm(wing_dat[[trait]] ~ wing_dat$LW))
}

wing_models <- fit_trait_models(wing_dat, wing_correct_traits)

table_s7 <- bind_rows(lapply(names(wing_models), function(trait) {
  car::Anova(wing_models[[trait]], type = 3) |>
    as.data.frame() |>
    rownames_to_column("Effect") |>
    mutate(Trait = trait, .before = 1)
}))
write.csv(table_s7, "outputs/tables/Table_S7_wing_corrected_typeIII_ANOVA.csv",
          row.names = FALSE)

wing_means <- bind_rows(lapply(names(wing_models), function(trait) {
  as.data.frame(emmeans(wing_models[[trait]], ~ Sex * Selection | Treatment)) |>
    mutate(Trait = trait)
}))

fig_s4 <- ggplot(
  wing_means,
  aes(Selection, emmean, colour = Sex, shape = Sex, group = Sex)
) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    position = position_dodge(0.45), width = 0.25, linewidth = 0.35
  ) +
  geom_point(position = position_dodge(0.45), size = 1.8) +
  facet_grid(Trait ~ Treatment, scales = "free_y") +
  scale_x_discrete(labels = c(CU = "Crowding\nevolved", MB = "Control")) +
  labs(x = NULL, y = "Wing-size-corrected residual",
       colour = NULL, shape = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom")

ggsave("outputs/figures/Figure_S4_wing_corrected_traits.pdf",
       fig_s4, width = 9, height = 6)

# Trait-wise effect sizes and condition dependence 

trait_ssd <- bind_rows(lapply(names(trait_models), function(trait) {
  emm <- emmeans(trait_models[[trait]], ~ Sex | Treatment * Selection)
  as.data.frame(
    eff_size(
      emm,
      sigma = sigma(trait_models[[trait]]),
      edf = df.residual(trait_models[[trait]])
    )
  ) |>
    mutate(Trait = trait)
}))

trait_density <- bind_rows(lapply(names(trait_models), function(trait) {
  emm <- emmeans(trait_models[[trait]], ~ Treatment | Sex * Selection)
  as.data.frame(
    eff_size(
      emm,
      sigma = sigma(trait_models[[trait]]),
      edf = df.residual(trait_models[[trait]])
    )
  ) |>
    mutate(Trait = trait)
}))

write.csv(trait_ssd, "outputs/tables/traitwise_sexual_dimorphism_effect_sizes.csv",
          row.names = FALSE)
write.csv(trait_density, "outputs/tables/traitwise_density_effect_sizes.csv",
          row.names = FALSE)

# Figure S5.
fig_s5a <- ggplot(
  trait_ssd,
  aes(Treatment, effect.size, colour = Selection, shape = Selection)
) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    position = position_dodge(0.5), width = 0.4, linewidth = 0.35
  ) +
  geom_point(position = position_dodge(0.5), size = 1.7) +
  facet_wrap(~ Trait, nrow = 1) +
  scale_colour_manual(values = selection_palette, labels = selection_labels) +
  scale_shape_manual(values = c(CU = 16, MB = 17), labels = selection_labels) +
  labs(x = "Larval density treatment", y = "Cohen's d (Female - Male)",
       tag = "A", colour = NULL, shape = NULL) +
  theme_bw(base_size = 9)

fig_s5b <- ggplot(
  trait_density,
  aes(contrast, effect.size, colour = Selection, shape = Selection)
) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    position = position_dodge(0.5), width = 0.4, linewidth = 0.35
  ) +
  geom_point(position = position_dodge(0.5), size = 1.7) +
  facet_grid(Sex ~ Trait) +
  scale_colour_manual(values = selection_palette, labels = selection_labels) +
  scale_shape_manual(values = c(CU = 16, MB = 17), labels = selection_labels) +
  labs(x = "Density contrast", y = "Cohen's d", tag = "B",
       colour = NULL, shape = NULL) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig_s5 <- fig_s5a / fig_s5b +
  plot_layout(heights = c(1, 2), guides = "collect") &
  theme(legend.position = "bottom")

ggsave("outputs/figures/Figure_S5_trait_effect_sizes.pdf",
       fig_s5, width = 12, height = 8)

# Association between sexual dimorphism and density sensitivity.
ssd_summary <- trait_ssd |>
  group_by(Trait) |>
  summarise(ssd = mean(abs(effect.size), na.rm = TRUE), .groups = "drop")

density_summary <- trait_density |>
  group_by(Trait) |>
  summarise(density_sensitivity = mean(abs(effect.size), na.rm = TRUE),
            .groups = "drop")

condition_dependence <- left_join(ssd_summary, density_summary, by = "Trait")
spearman_test <- cor.test(
  condition_dependence$ssd,
  condition_dependence$density_sensitivity,
  method = "spearman",
  exact = FALSE
)

write.csv(condition_dependence,
          "outputs/tables/trait_dimorphism_density_sensitivity.csv",
          row.names = FALSE)
capture.output(
  spearman_test,
  file = "outputs/tables/trait_dimorphism_density_spearman.txt"
)

# Session information 
capture.output(sessionInfo(), file = "outputs/sessionInfo.txt")
