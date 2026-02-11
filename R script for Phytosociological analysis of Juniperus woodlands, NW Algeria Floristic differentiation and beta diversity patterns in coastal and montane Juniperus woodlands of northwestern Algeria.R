# Phytosociological analysis of Juniperus woodlands, NW Algeria
# Floristic differentiation and beta diversity patterns in coastal
# and montane Juniperus woodlands of northwestern Algeria

rm(list = ls())
RNGkind(kind = "Mersenne-Twister", sample.kind = "Rejection")
set.seed(2024)

data_path <- "C:/Users/berber.medelmehdi/Downloads/phytosociologie/"
if (!endsWith(data_path, "/") && !endsWith(data_path, "\\")) {
  data_path <- paste0(data_path, "/")
}

packages <- c("readxl", "dplyr", "tidyr", "vegan", "cluster",
              "ggplot2", "ggrepel", "indicspecies", "betapart",
              "writexl", "permute")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# 1. Data preparation
# ============================================================================

ad  <- read_excel(paste0(data_path, "ad_final.xlsx"))
pa  <- read_excel(paste0(data_path, "pa_final.xlsx"))
inv <- read_excel(paste0(data_path, "inv_final.xlsx"))
env_raw <- read_excel(paste0(data_path, "env_final.xlsx"))

releve_cols <- paste0("R", 1:23)

# Station-level aggregation by maximum abundance-dominance
mat_station <- ad %>%
  rowwise() %>%
  mutate(abondance_max = max(c_across(all_of(releve_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(code_station, code_taxon, abondance_max) %>%
  group_by(code_station, code_taxon) %>%
  summarise(abondance = max(abondance_max, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = code_taxon, values_from = abondance, values_fill = 0)

station_names <- mat_station$code_station
mat <- mat_station %>% select(-code_station) %>% as.data.frame()
rownames(mat) <- station_names
mat <- mat[, colSums(mat) > 0]

n_stations <- nrow(mat)
n_species  <- ncol(mat)

# Environmental data
if ("code_station" %in% names(env_raw)) {
  env <- env_raw %>%
    filter(code_station %in% rownames(mat)) %>%
    distinct(code_station, .keep_all = TRUE)
} else {
  env <- ad %>%
    distinct(code_station, .keep_all = TRUE) %>%
    filter(code_station %in% rownames(mat))
}
env <- env[match(rownames(mat), env$code_station), ]

# A priori bioclimatic zones
stations_montaneane <- c("ZRF", "MFR", "TRN")
env$zone <- ifelse(env$code_station %in% stations_montaneane, "montaneane", "coastal")
env$zone <- factor(env$zone, levels = c("coastal", "montaneane"))

# Transformations
mat_hellinger <- decostand(mat, method = "hellinger")
mat_pa <- (mat > 0) * 1

dict_taxons <- ad %>% distinct(code_taxon, taxons) %>% arrange(code_taxon)
richesse <- rowSums(mat_pa)

# ============================================================================
# 1bis. Sensitivity analysis — aggregation method
# ============================================================================

aggregate_station_matrix <- function(ad_data, releve_cols, method = "max") {
  ad_agg <- ad_data %>%
    rowwise() %>%
    mutate(agg_value = switch(method,
      "max"    = max(c_across(all_of(releve_cols)), na.rm = TRUE),
      "mean"   = mean(c_across(all_of(releve_cols)), na.rm = TRUE),
      "median" = median(c_across(all_of(releve_cols)), na.rm = TRUE)
    )) %>%
    ungroup() %>%
    select(code_station, code_taxon, agg_value) %>%
    group_by(code_station, code_taxon) %>%
    summarise(abondance = switch(method,
      "max"    = max(agg_value, na.rm = TRUE),
      "mean"   = mean(agg_value, na.rm = TRUE),
      "median" = median(agg_value, na.rm = TRUE)
    ), .groups = "drop") %>%
    pivot_wider(names_from = code_taxon, values_from = abondance, values_fill = 0)

  sn <- ad_agg$code_station
  m  <- ad_agg %>% select(-code_station) %>% as.data.frame()
  rownames(m) <- sn
  m <- m[, colSums(m) > 0]
  return(m)
}

compute_ari <- function(c1, c2) {
  tab <- table(c1, c2)
  n <- sum(tab)
  sum_ni2  <- sum(choose(rowSums(tab), 2))
  sum_nj2  <- sum(choose(colSums(tab), 2))
  sum_nij2 <- sum(choose(tab, 2))
  expected  <- sum_ni2 * sum_nj2 / choose(n, 2)
  max_index <- 0.5 * (sum_ni2 + sum_nj2)
  if (max_index == expected) return(1)
  (sum_nij2 - expected) / (max_index - expected)
}

grp_ref <- cutree(
  hclust(dist(decostand(mat, method = "hellinger"), method = "euclidean"),
         method = "ward.D2"), k = 2)

sensitivity_results <- list()
for (m in c("max", "mean", "median")) {
  mat_m      <- aggregate_station_matrix(ad, releve_cols, m)
  mat_m_hell <- decostand(mat_m, method = "hellinger")
  dist_m     <- dist(mat_m_hell, method = "euclidean")
  hc_m       <- hclust(dist_m, method = "ward.D2")
  grp_m      <- cutree(hc_m, k = 2)
  sil_m      <- mean(silhouette(grp_m, dist_m)[, 3])
  ari_m      <- compute_ari(grp_ref, grp_m)

  zone_m <- factor(ifelse(rownames(mat_m) %in% stations_montaneane, "montaneane", "coastal"),
                   levels = c("coastal", "montaneane"))
  env_m  <- data.frame(zone = zone_m, row.names = rownames(mat_m))

  set.seed(2024)
  perm_m <- adonis2(mat_m ~ zone, data = env_m, method = "bray", permutations = 9999)

  sensitivity_results[[m]] <- list(
    n_species = ncol(mat_m), silhouette = sil_m, ari = ari_m,
    R2 = perm_m$R2[1], F_stat = perm_m$F[1], p_value = perm_m$`Pr(>F)`[1])
}

sensitivity_table <- data.frame(
  Method       = c("Maximum", "Mean", "Median"),
  N_species    = sapply(sensitivity_results, function(x) x$n_species),
  Silhouette   = sapply(sensitivity_results, function(x) round(x$silhouette, 3)),
  ARI_vs_max   = sapply(sensitivity_results, function(x) round(x$ari, 3)),
  PERMANOVA_R2 = sapply(sensitivity_results, function(x) round(x$R2, 3)),
  PERMANOVA_F  = sapply(sensitivity_results, function(x) round(x$F_stat, 2)),
  PERMANOVA_p  = sapply(sensitivity_results, function(x) round(x$p_value, 4)),
  row.names = NULL)
print(sensitivity_table)

rm(mat_m, mat_m_hell, dist_m, hc_m, grp_m, sil_m, ari_m, zone_m, env_m, perm_m, m)

# ============================================================================
# 2. Hierarchical classification
# ============================================================================

dist_eucl <- dist(mat_hellinger, method = "euclidean")
hc_ward   <- hclust(dist_eucl, method = "ward.D2")

sil_scores <- numeric(5)
for (k in 2:5) {
  grp <- cutree(hc_ward, k = k)
  sil_scores[k] <- mean(silhouette(grp, dist_eucl)[, 3])
}
k_optimal   <- which.max(sil_scores)
sil_optimal <- sil_scores[k_optimal]

groupes   <- cutree(hc_ward, k = k_optimal)
tab_corr  <- table(Cluster = groupes, Zone = env$zone)
pct_corr  <- sum(diag(tab_corr)) / sum(tab_corr) * 100
sil_final <- silhouette(groupes, dist_eucl)

# ============================================================================
# 3. Ordination (DCA, NMDS, CCA)
# ============================================================================

# 3.1 DCA
dca_result      <- decorana(mat)
dca_site_scores <- scores(dca_result, display = "sites", choices = 1:4)
gradient_length <- diff(range(dca_site_scores[, 1]))
if (!is.null(attr(dca_result, "axlen"))) {
  gradient_length <- attr(dca_result, "axlen")[1]
}
dca_eig <- dca_result$evals
dca_var <- dca_eig / sum(dca_eig) * 100

# 3.2 NMDS
set.seed(2024)
nmds_result <- metaMDS(mat, distance = "bray", k = 2, trymax = 100,
                       autotransform = FALSE, trace = 0)

# 3.3 CCA
if ("altitude_m" %in% names(env)) {
  altitude <- env$altitude_m
} else if ("altitude" %in% names(env)) {
  altitude <- env$altitude
} else {
  altitude <- c(SS = 60, BSP1 = 45, BSP2 = 65, MD = 65, SFX = 40,
                RCHP = 2, RCH22 = 35, SSF = 128, MSC = 35, MCTR = 60, CR = 55,
                ZRF = 1105, MFR = 1150, TRN = 1335)
  altitude <- altitude[rownames(mat)]
}
if ("pente_mean" %in% names(env)) {
  pente <- env$pente_mean
} else if ("pente" %in% names(env)) {
  pente <- env$pente
} else {
  pente <- c(SS = 25, BSP1 = 22.5, BSP2 = 27.5, MD = 10, SFX = 15,
             RCHP = 2.5, RCH22 = 32.5, SSF = 13, MSC = 25, MCTR = 32.5, CR = 35,
             ZRF = 21, MFR = 10, TRN = 40)
  pente <- pente[rownames(mat)]
}

env_num <- data.frame(altitude_m = as.numeric(altitude),
                      pente_mean = as.numeric(pente),
                      row.names = rownames(mat))

set.seed(2024)
cca_result <- cca(mat ~ altitude_m + pente_mean, data = env_num)

cca_total         <- cca_result$tot.chi
cca_constrained   <- sum(cca_result$CCA$eig)
cca_var_explained <- cca_constrained / cca_total * 100
cca_eig           <- cca_result$CCA$eig
cca_var_axes      <- cca_eig / sum(cca_eig) * 100

set.seed(2024)
cca_anova_global <- anova(cca_result, permutations = 9999)
cca_p_global     <- cca_anova_global$`Pr(>F)`[1]

set.seed(2024)
cca_anova_margin <- anova(cca_result, by = "margin", permutations = 9999)

vif_values   <- vif.cca(cca_result)
cor_alt_zone <- cor(env_num$altitude_m, as.numeric(env$zone == "montaneane"))

# ============================================================================
# 4. Multivariate tests (a priori zones)
# ============================================================================

dist_bray <- vegdist(mat, method = "bray")

# ANOSIM
set.seed(2024)
anosim_result <- anosim(mat, env$zone, distance = "bray", permutations = 9999)

# PERMANOVA
set.seed(2024)
permanova_result <- adonis2(mat ~ zone, data = env, method = "bray", permutations = 9999)

# BETADISPER
betadisper_result <- betadisper(dist_bray, env$zone)
set.seed(2024)
betadisper_test <- permutest(betadisper_result, permutations = 9999)
betadisper_p    <- betadisper_test$tab$`Pr(>F)`[1]
betadisper_F    <- betadisper_test$tab$F[1]

# ============================================================================
# 5. Indicator species (IndVal on a priori zones)
# ============================================================================

set.seed(2024)
indval_result <- multipatt(mat, env$zone, func = "IndVal.g",
                           control = how(nperm = 9999))

indval_summary <- indval_result$sign
indval_summary$code_taxon <- rownames(indval_summary)
indval_fdr <- indval_summary %>% filter(p.value < 0.05) %>% arrange(p.value)

n_indic_total     <- nrow(indval_fdr)
n_indic_coastal   <- sum(indval_fdr$index == 1)
n_indic_montaneane <- sum(indval_fdr$index == 2)

get_taxon_name <- function(code) {
  nom <- dict_taxons$taxons[match(code, dict_taxons$code_taxon)]
  if (is.na(nom) || length(nom) == 0) return(code)
  nom[1]
}

# FDR correction (Benjamini-Hochberg)
indval_all <- indval_result$sign
indval_fdr <- indval_all[!is.na(indval_all$p.value) & indval_all$p.value < 0.05, ]
indval_fdr$p.adj      <- p.adjust(indval_fdr$p.value, method = "BH")
indval_fdr             <- indval_fdr[order(indval_fdr$p.value), ]
indval_fdr$code_taxon <- rownames(indval_fdr)
indval_fdr             <- merge(indval_fdr, dict_taxons, by = "code_taxon", all.x = TRUE)
indval_fdr$zone <- ifelse(indval_fdr$s.coastal == 1 & indval_fdr$s.montaneane == 0, "coastal",
                   ifelse(indval_fdr$s.montaneane == 1 & indval_fdr$s.coastal == 0, "montaneane", "both"))

write_xlsx(indval_fdr, paste0(data_path, "indval_results_with_padj.xlsx"))

# ============================================================================
# 6. Diversity indices
# ============================================================================

# 6.1 Alpha diversity
div_alpha <- data.frame(
  station = rownames(mat), zone = env$zone,
  S = rowSums(mat_pa),
  H = diversity(mat, index = "shannon"),
  D = diversity(mat, index = "simpson"))

div_summary <- div_alpha %>%
  group_by(zone) %>%
  summarise(n = n(),
    S_mean = round(mean(S), 1), S_sd = round(sd(S), 1),
    H_mean = round(mean(H), 2), H_sd = round(sd(H), 2),
    D_mean = round(mean(D), 3), D_sd = round(sd(D), 3),
    .groups = "drop")

# 6.2 Beta diversity — multi-site partitioning (Baselga)
gamma_div     <- ncol(mat_pa)
alpha_mean    <- mean(rowSums(mat_pa))
beta_whittaker <- gamma_div / alpha_mean - 1

beta_pair <- beta.pair(mat_pa, index.family = "sorensen")
beta_sor  <- mean(as.dist(beta_pair$beta.sor))
beta_sim  <- mean(as.dist(beta_pair$beta.sim))
beta_sne  <- mean(as.dist(beta_pair$beta.sne))

turnover_pct   <- beta_sim / beta_sor * 100
nestedness_pct <- beta_sne / beta_sor * 100

# 6.3 LCBD (Jaccard — Legendre & De Cáceres, 2013)
dist_jaccard <- vegdist(mat_pa, method = "jaccard", binary = TRUE)
D_mat    <- as.matrix(dist_jaccard)
D2       <- D_mat^2
n_lcbd   <- nrow(D2)
I_mat    <- diag(n_lcbd)
ones_mat <- matrix(1, n_lcbd, n_lcbd)
A_mat    <- -0.5 * D2
G_mat    <- (I_mat - ones_mat/n_lcbd) %*% A_mat %*% (I_mat - ones_mat/n_lcbd)

lcbd_values <- diag(G_mat) / sum(diag(G_mat))
names(lcbd_values) <- rownames(mat_pa)

# Permutation test
set.seed(2024)
nperm_lcbd    <- 9999
lcbd_perm_mat <- matrix(NA, nrow = nperm_lcbd, ncol = n_lcbd)
for (k in 1:nperm_lcbd) {
  mat_perm <- mat_pa[sample(n_lcbd), ]
  rownames(mat_perm) <- rownames(mat_pa)
  dist_perm <- vegdist(mat_perm, method = "jaccard", binary = TRUE)
  D2_perm   <- as.matrix(dist_perm)^2
  A_perm    <- -0.5 * D2_perm
  G_perm    <- (I_mat - ones_mat/n_lcbd) %*% A_perm %*% (I_mat - ones_mat/n_lcbd)
  lcbd_perm_mat[k, ] <- diag(G_perm) / sum(diag(G_perm))
}

lcbd_pvalues <- sapply(1:n_lcbd, function(i) {
  (sum(lcbd_perm_mat[, i] >= lcbd_values[i]) + 1) / (nperm_lcbd + 1)
})
names(lcbd_pvalues) <- rownames(mat_pa)

lcbd_df <- data.frame(
  station = names(lcbd_values), zone = env$zone, S = rowSums(mat_pa),
  LCBD = lcbd_values, p_value = round(lcbd_pvalues, 4),
  signif = ifelse(lcbd_pvalues < 0.05, "*", "")) %>%
  arrange(desc(LCBD))

lcbd_montaneane_total <- sum(lcbd_df$LCBD[lcbd_df$zone == "montaneane"])

cor_lcbd_rich <- cor.test(lcbd_values, richesse, method = "spearman")

rm(D_mat, D2, n_lcbd, I_mat, ones_mat, A_mat, G_mat, lcbd_perm_mat, dist_jaccard, nperm_lcbd)

# ============================================================================
# 7. Figures
# ============================================================================

couleurs <- c("coastal" = "#E41A1C", "montaneane" = "#377EB8")

# Fig 1: Dendrogram
png(paste0(data_path, "Fig1_dendrogram.png"), width = 1200, height = 800, res = 150)
par(mar = c(5, 4, 4, 2))
plot(hc_ward,
     main = paste0("Hierarchical clustering (k = ", k_optimal,
                   ", silhouette = ", round(sil_optimal, 3), ")"),
     sub = "Hellinger transformation + Euclidean distance + Ward.D2",
     xlab = "Stations", ylab = "Distance", hang = -1)
rect.hclust(hc_ward, k = k_optimal, border = c("#E41A1C", "#377EB8"))
dev.off()

# Fig 2: NMDS
nmds_scores <- data.frame(scores(nmds_result, display = "sites"),
                          station = rownames(mat), zone = env$zone)

fig_nmds <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = zone, shape = zone)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text_repel(aes(label = station), size = 4, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = couleurs, name = "Bioclimatic zone") +
  scale_shape_manual(values = c(coastal = 16, montaneane = 17), name = "Bioclimatic zone") +
  labs(title = "NMDS ordination of stations",
       subtitle = paste("Stress =", round(nmds_result$stress, 3), "(Bray-Curtis)"),
       x = "NMDS1", y = "NMDS2") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") + coord_fixed()
ggsave(paste0(data_path, "Fig2_NMDS.png"), fig_nmds, width = 10, height = 8, dpi = 300)

# Fig 3: CCA biplot
cca_sites <- data.frame(scores(cca_result, display = "sites", choices = 1:2),
                        station = rownames(mat), zone = env$zone)
cca_bp <- data.frame(scores(cca_result, display = "bp", choices = 1:2),
                     variable = rownames(scores(cca_result, display = "bp", choices = 1:2)))

fig_cca <- ggplot() +
  geom_point(data = cca_sites, aes(x = CCA1, y = CCA2, color = zone, shape = zone),
             size = 5, alpha = 0.8) +
  geom_text_repel(data = cca_sites, aes(x = CCA1, y = CCA2, label = station, color = zone),
                  size = 4, max.overlaps = 20, show.legend = FALSE) +
  geom_segment(data = cca_bp,
               aes(x = 0, y = 0, xend = CCA1 * 1.5, yend = CCA2 * 1.5),
               arrow = arrow(length = unit(0.3, "cm")), color = "darkred", linewidth = 1) +
  geom_text(data = cca_bp,
            aes(x = CCA1 * 1.7, y = CCA2 * 1.7, label = variable),
            color = "darkred", fontface = "bold", size = 4) +
  scale_color_manual(values = couleurs, name = "Bioclimatic zone") +
  scale_shape_manual(values = c(coastal = 16, montaneane = 17), name = "Bioclimatic zone") +
  labs(title = "CCA — Vegetation-environment relationships",
       subtitle = paste0("Constrained variance = ", round(cca_var_explained, 1),
                         "% (p = ", cca_p_global, ")"),
       x = paste0("CCA1 (", round(cca_var_axes[1], 1), "%)"),
       y = paste0("CCA2 (", round(cca_var_axes[2], 1), "%)")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) + coord_fixed()
ggsave(paste0(data_path, "Fig3_CCA.png"), fig_cca, width = 10, height = 8, dpi = 300)

# Fig 4: Species richness by zone
fig_richesse <- ggplot(div_alpha, aes(x = zone, y = S, fill = zone)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
  scale_fill_manual(values = couleurs) +
  labs(title = "Species richness by bioclimatic zone",
       x = "Bioclimatic zone", y = "Number of species (S)") +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(paste0(data_path, "Fig4_richesse.png"), fig_richesse, width = 8, height = 6, dpi = 300)

# Fig 5: Beta diversity partitioning
beta_data <- data.frame(
  Component  = c("Turnover (βSIM)", "Nestedness (βSNE)"),
  Value      = c(beta_sim, beta_sne),
  Percentage = c(turnover_pct, nestedness_pct))

fig_beta <- ggplot(beta_data, aes(x = "", y = Value, fill = Component)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Turnover (βSIM)" = "#E41A1C", "Nestedness (βSNE)" = "#377EB8")) +
  labs(title = "Beta diversity partitioning",
       subtitle = paste0("Total Sørensen = ", round(beta_sor, 3))) +
  theme_void() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
ggsave(paste0(data_path, "Fig5_beta_partition.png"), fig_beta, width = 8, height = 8, dpi = 300)

# Fig 6: LCBD
lcbd_df$station <- factor(lcbd_df$station, levels = lcbd_df$station)

fig_lcbd <- ggplot(lcbd_df, aes(x = station, y = LCBD, fill = zone)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = couleurs, name = "Zone") +
  labs(title = "Local Contribution to Beta Diversity (LCBD)",
       subtitle = paste0("Montaneane stations: ", round(lcbd_montaneane_total * 100, 1),
                         "% of total beta diversity"),
       x = "Station", y = "LCBD") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 1/n_stations, linetype = "dashed", color = "red")
ggsave(paste0(data_path, "Fig6_LCBD.png"), fig_lcbd, width = 12, height = 6, dpi = 300)

# ============================================================================
# 8. Export
# ============================================================================

save(
  mat, mat_hellinger, mat_pa, env, env_num, dict_taxons,
  hc_ward, dist_eucl, dist_bray, groupes, k_optimal, sil_optimal, sil_final,
  dca_result, nmds_result, cca_result, gradient_length,
  cca_var_explained, cca_p_global, cca_anova_margin, vif_values, cor_alt_zone,
  anosim_result, permanova_result, betadisper_result, betadisper_p,
  indval_result, indval_fdr, n_indic_total, n_indic_coastal, n_indic_montaneane,
  div_alpha, div_summary, gamma_div, alpha_mean, beta_whittaker,
  beta_sor, beta_sim, beta_sne, turnover_pct, nestedness_pct,
  lcbd_df, lcbd_montaneane_total,
  file = paste0(data_path, "resultats_article.RData"))

write_xlsx(list(
  Stations = data.frame(
    station = rownames(mat), zone = as.character(env$zone),
    altitude_m = env_num$altitude_m, pente_mean = env_num$pente_mean,
    S = rowSums(mat_pa),
    H = round(diversity(mat, "shannon"), 2),
    D = round(diversity(mat, "simpson"), 3)),
  Tests = data.frame(
    Test = c("ANOSIM R", "ANOSIM p", "PERMANOVA R²", "PERMANOVA p",
             "BETADISPER p", "DCA gradient (SD)",
             "CCA constrained variance (%)", "CCA p-value",
             "Correlation altitude-zone"),
    Value = c(round(anosim_result$statistic, 3), anosim_result$signif,
              round(permanova_result$R2[1], 3), permanova_result$`Pr(>F)`[1],
              round(betadisper_p, 3), round(gradient_length, 2),
              round(cca_var_explained, 1), cca_p_global,
              round(cor_alt_zone, 3))),
  Beta_Diversity = data.frame(
    Index = c("Gamma", "Alpha mean", "Beta Whittaker",
              "Beta Sørensen", "Turnover (βSIM)", "Nestedness (βSNE)",
              "Turnover (%)", "Nestedness (%)"),
    Value = c(gamma_div, round(alpha_mean, 1), round(beta_whittaker, 2),
              round(beta_sor, 3), round(beta_sim, 3), round(beta_sne, 3),
              round(turnover_pct, 1), round(nestedness_pct, 1))),
  LCBD = lcbd_df,
  Indicators_Coastal = indval_fdr %>%
    filter(index == 1) %>%
    mutate(taxon = sapply(code_taxon, get_taxon_name)) %>%
    select(code_taxon, taxon, stat, p.value) %>% arrange(desc(stat)),
  Indicators_Montaneane = indval_fdr %>%
    filter(index == 2) %>%
    mutate(taxon = sapply(code_taxon, get_taxon_name)) %>%
    select(code_taxon, taxon, stat, p.value) %>% arrange(desc(stat)),
  Silhouette = data.frame(k = 2:5, silhouette = round(sil_scores[2:5], 3))
), path = paste0(data_path, "resultats_article.xlsx"))


# ============================================================================
# Supplementary analyses
# ============================================================================

# ============================================================================
# A. Raup-Crick (βRC) — determinism vs. stochasticity
# ============================================================================

set.seed(2024)
rc_matrix <- raupcrick(mat_pa, nsimul = 9999, chase = TRUE)
rc_mat    <- as.matrix(rc_matrix)

rc_values        <- as.dist(rc_matrix)
rc_mean_global   <- mean(rc_values)
rc_sd_global     <- sd(rc_values)
rc_median_global <- median(rc_values)

idx_coastal    <- which(env$zone == "coastal")
idx_montaneane <- which(env$zone == "montaneane")

get_pairwise_values <- function(mat, idx1, idx2 = NULL) {
  if (is.null(idx2)) idx2 <- idx1
  values <- c()
  for (i in idx1) for (j in idx2) if (i < j) values <- c(values, mat[i, j])
  values
}

rc_intra_coastal    <- get_pairwise_values(rc_mat, idx_coastal)
rc_intra_montaneane <- get_pairwise_values(rc_mat, idx_montaneane)
rc_inter <- c()
for (i in idx_coastal) for (j in idx_montaneane) rc_inter <- c(rc_inter, rc_mat[i, j])

rc_intra_all <- c(rc_intra_coastal, rc_intra_montaneane)
if (length(rc_intra_all) > 2 && length(rc_inter) > 2) {
  wilcox_rc <- wilcox.test(rc_inter, rc_intra_all, alternative = "greater")
}

rc_summary <- data.frame(
  Comparison = c("Intra-coastal", "Intra-montaneane", "Inter-zone", "Global"),
  n_pairs = c(length(rc_intra_coastal), length(rc_intra_montaneane),
              length(rc_inter), length(rc_values)),
  Mean = c(round(mean(rc_intra_coastal), 3),
           ifelse(length(rc_intra_montaneane) > 0, round(mean(rc_intra_montaneane), 3), NA),
           round(mean(rc_inter), 3), round(rc_mean_global, 3)),
  Median = c(round(median(rc_intra_coastal), 3),
             ifelse(length(rc_intra_montaneane) > 0, round(median(rc_intra_montaneane), 3), NA),
             round(median(rc_inter), 3), round(rc_median_global, 3)),
  SD = c(round(sd(rc_intra_coastal), 3),
         ifelse(length(rc_intra_montaneane) > 1, round(sd(rc_intra_montaneane), 3), NA),
         round(sd(rc_inter), 3), round(rc_sd_global, 3)))

# ============================================================================
# B. Biological spectrum (Raunkiaer life forms)
# ============================================================================

# Requires a "type_biologique" column in dict_taxons.
# If unavailable, life forms must be coded from Quézel & Santa.
# Template code provided below for use once data are available.

# [Template omitted — see supplementary materials if needed]

# ============================================================================
# C. Evenness and dominance
# ============================================================================

H_shannon    <- diversity(mat, index = "shannon")
S_richness   <- specnumber(mat)
J_pielou     <- H_shannon / log(S_richness)
berger_parker <- apply(mat, 1, function(x) max(x) / sum(x))
D_simpson    <- diversity(mat, index = "simpson")
lambda_simpson <- 1 - D_simpson
N1_hill      <- exp(H_shannon)
N2_hill      <- 1 / (1 - D_simpson)
hill_ratio   <- N2_hill / N1_hill

diversity_indices <- data.frame(
  station = rownames(mat), zone = env$zone,
  S = S_richness, H = round(H_shannon, 3), J = round(J_pielou, 3),
  D = round(D_simpson, 4), BP = round(berger_parker, 4),
  N1 = round(N1_hill, 1), N2 = round(N2_hill, 1), N2_N1 = round(hill_ratio, 3))

diversity_summary <- diversity_indices %>%
  group_by(zone) %>%
  summarise(n = n(),
    S_mean = round(mean(S), 1), S_sd = round(sd(S), 1),
    J_mean = round(mean(J), 3), J_sd = round(sd(J), 3),
    BP_mean = round(mean(BP), 4), BP_sd = round(sd(BP), 4),
    N1_mean = round(mean(N1), 1), N2_mean = round(mean(N2), 1),
    .groups = "drop")

# Wilcoxon tests
indices_to_test <- c("S", "H", "D", "J", "BP", "N1", "N2", "N2_N1")
test_results <- data.frame(Index = character(), Coastal_mean = numeric(),
  Montaneane_mean = numeric(), W = numeric(), p_value = numeric(),
  stringsAsFactors = FALSE)

for (idx in indices_to_test) {
  cv <- diversity_indices[[idx]][diversity_indices$zone == "coastal"]
  mv <- diversity_indices[[idx]][diversity_indices$zone == "montaneane"]
  tt <- wilcox.test(cv, mv, exact = FALSE, conf.int = FALSE)
  test_results <- rbind(test_results, data.frame(
    Index = idx, Coastal_mean = round(mean(cv), 3),
    Montaneane_mean = round(mean(mv), 3),
    W = round(tt$statistic, 1), p_value = round(tt$p.value, 3)))
}

# ============================================================================
# D. Supplementary figures
# ============================================================================

# Raup-Crick heatmap
station_order <- c(rownames(mat)[env$zone == "coastal"],
                   rownames(mat)[env$zone == "montaneane"])
rc_df <- expand.grid(Station1 = rownames(mat), Station2 = rownames(mat))
rc_df$betaRC  <- as.vector(rc_mat)
rc_df$Zone1   <- env$zone[match(rc_df$Station1, rownames(mat))]
rc_df$Zone2   <- env$zone[match(rc_df$Station2, rownames(mat))]
rc_df$Station1 <- factor(rc_df$Station1, levels = station_order)
rc_df$Station2 <- factor(rc_df$Station2, levels = station_order)

fig_rc_heatmap <- ggplot(rc_df, aes(x = Station1, y = Station2, fill = betaRC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#377EB8", mid = "white", high = "#E41A1C",
                       midpoint = 0.5, limits = c(0, 1), name = "βRC") +
  labs(title = "Raup-Crick dissimilarity matrix",
       subtitle = paste0("Mean βRC = ", round(rc_mean_global, 3),
                         "; Inter-zone mean = ", round(mean(rc_inter), 3)),
       x = "", y = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  coord_fixed() +
  geom_vline(xintercept = sum(env$zone == "coastal") + 0.5, linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = sum(env$zone == "coastal") + 0.5, linetype = "dashed", linewidth = 1)
ggsave(paste0(data_path, "Fig_RaupCrick_heatmap.png"), fig_rc_heatmap, width = 10, height = 9, dpi = 300)

# Raup-Crick boxplot by comparison type
rc_comparison <- data.frame(
  Type = factor(c(rep("Intra-coastal", length(rc_intra_coastal)),
                  rep("Intra-montaneane", length(rc_intra_montaneane)),
                  rep("Inter-zones", length(rc_inter))),
                levels = c("Intra-coastal", "Intra-montaneane", "Inter-zones")),
  betaRC = c(rc_intra_coastal, rc_intra_montaneane, rc_inter))

fig_rc_boxplot <- ggplot(rc_comparison, aes(x = Type, y = betaRC, fill = Type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("Intra-coastal" = "#E41A1C",
                               "Intra-montaneane" = "#377EB8",
                               "Inter-zones" = "#984EA3")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  labs(title = "Raup-Crick index by comparison type",
       subtitle = "Values > 0.5 indicate deterministic assembly",
       x = "Comparison type", y = "βRC") +
  theme_bw(base_size = 12) + theme(legend.position = "none") + ylim(0, 1)
ggsave(paste0(data_path, "Fig_RaupCrick_boxplot.png"), fig_rc_boxplot, width = 8, height = 6, dpi = 300)

# Evenness by zone
fig_equitability <- ggplot(diversity_indices, aes(x = zone, y = J, fill = zone)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("coastal" = "#E41A1C", "montaneane" = "#377EB8")) +
  labs(title = "Pielou's evenness (J') by bioclimatic zone",
       x = "Bioclimatic zone", y = "Pielou's J'") +
  theme_bw(base_size = 12) + theme(legend.position = "none") + ylim(0.85, 1)
ggsave(paste0(data_path, "Fig_equitability.png"), fig_equitability, width = 8, height = 6, dpi = 300)

# Richness vs Evenness
fig_S_J <- ggplot(diversity_indices, aes(x = S, y = J, color = zone, shape = zone)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", aes(group = 1), color = "gray50") +
  geom_text_repel(aes(label = station), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("coastal" = "#E41A1C", "montaneane" = "#377EB8"), name = "Zone") +
  scale_shape_manual(values = c(coastal = 16, montaneane = 17), name = "Zone") +
  labs(title = "Relationship between species richness and evenness",
       x = "Species richness (S)", y = "Pielou's evenness (J')") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")
ggsave(paste0(data_path, "Fig_richness_evenness.png"), fig_S_J, width = 10, height = 8, dpi = 300)

# ============================================================================
# E. Export supplementary results
# ============================================================================

write_xlsx(list(
  RaupCrick_Summary = rc_summary,
  RaupCrick_Wilcoxon = data.frame(
      Test = "Inter-zone vs Intra-zone (one-sided)",
      W = wilcox_rc$statistic,
      p_value = wilcox_rc$p.value),
  RaupCrick_Matrix  = data.frame(Station = rownames(rc_mat), as.data.frame(round(rc_mat, 3))),
  Diversity_Indices = diversity_indices,
  Diversity_Summary = diversity_summary,
  Diversity_Tests   = test_results
), path = paste0(data_path, "resultats_complementaires.xlsx"))

save(rc_matrix, rc_mat, rc_mean_global, rc_summary,
     rc_intra_coastal, rc_intra_montaneane, rc_inter, wilcox_rc,
     J_pielou, berger_parker, N1_hill, N2_hill, hill_ratio,
     diversity_indices, diversity_summary, test_results,
     file = paste0(data_path, "resultats_complementaires.RData"))


# ============================================================================
# TRY database and CSI analysis
# ============================================================================

packages_csi <- c("readxl", "writexl", "dplyr", "tidyr", "rtry",
                  "ade4", "vegan", "ggplot2")
for (pkg in packages_csi) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

set.seed(2024)

# ============================================================================
# F. TRY data import and processing
# ============================================================================

data_try <- rtry_import(paste0(data_path, "46739.txt"))

data_traits <- data_try[!is.na(data_try$TraitID) &
                          (is.na(data_try$ErrorRisk) | data_try$ErrorRisk <= 4), ]

traits_summary <- data_traits %>%
  filter(!is.na(StdValue)) %>%
  group_by(AccSpeciesName, TraitID, TraitName) %>%
  summarise(mean_value = mean(StdValue, na.rm = TRUE), n_obs = n(), .groups = "drop")

traits_wide <- traits_summary %>%
  select(AccSpeciesName, TraitName, mean_value) %>%
  pivot_wider(names_from = TraitName, values_from = mean_value)

traits_clean <- traits_wide %>%
  rename(
    Species = AccSpeciesName,
    Seed_mass_mg = `Seed dry mass`,
    Height_veg_m = `Plant height vegetative`,
    Height_gen_m = `Plant height generative`,
    Leaf_N_mg_g = `Leaf nitrogen (N) content per leaf dry mass`,
    LDMC_g_g = `Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)`,
    SLA_petiole_excl = `Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded`,
    SLA_petiole_incl = `Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included`,
    SLA_undefined = `Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded`,
    Wood_density = `Stem specific density (SSD, stem dry mass per stem fresh volume) or wood density`
  ) %>%
  mutate(SLA_mm2_mg = coalesce(SLA_petiole_excl, SLA_petiole_incl, SLA_undefined)) %>%
  select(Species, Height_veg_m, Height_gen_m, SLA_mm2_mg, Leaf_N_mg_g,
         LDMC_g_g, Seed_mass_mg, Wood_density)

write_xlsx(traits_clean, paste0(data_path, "traits_TRY_clean.xlsx"))
write.csv(traits_clean, paste0(data_path, "traits_TRY_clean.csv"), row.names = FALSE)

# ============================================================================
# G. Community data preparation
# ============================================================================

ad  <- read_excel(paste0(data_path, "ad_final.xlsx"))
env <- read_excel(paste0(data_path, "env_final.xlsx"))

ad$taxons <- trimws(gsub("\\s+", " ", ad$taxons))

map_to_station <- function(station_name) {
  sl <- tolower(trimws(station_name))
  if (grepl("^beni saf plage [0-9]+$", sl) && !grepl("^beni saf plage 1$", sl))
    return("beni saf plage 1")
  sl
}

ad$station_mapped  <- sapply(ad$station, map_to_station)
env$station_mapped <- tolower(trimws(env$station))

releve_cols <- grep("^R[0-9]+$", names(ad), value = TRUE)

ad_station <- ad %>%
  rowwise() %>%
  mutate(abondance_max = max(c_across(all_of(releve_cols)), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(station_mapped, taxons) %>%
  summarise(abondance = max(abondance_max, na.rm = TRUE), .groups = "drop")

mat_station <- ad_station %>%
  pivot_wider(names_from = taxons, values_from = abondance, values_fill = 0)

stations <- mat_station$station_mapped
mat_comm <- as.matrix(mat_station[, -1])
rownames(mat_comm) <- stations

env_aligned <- env[match(stations, env$station_mapped), ]

# ============================================================================
# H. OMI analysis (Outlying Mean Index — Dolédec et al., 2000)
# ============================================================================

env_vars <- as.data.frame(env_aligned[, c("altitude_m", "pente_mean")])
rownames(env_vars) <- stations
env_scaled <- scale(env_vars)
rownames(env_scaled) <- stations

mat_species <- t(mat_comm)
freq_species <- rowSums(mat_species > 0)
mat_species_filtered <- mat_species[freq_species >= 2, ]

pca_env    <- dudi.pca(as.data.frame(env_scaled), scannf = FALSE, nf = 2)
omi_result <- niche(pca_env, as.data.frame(t(mat_species_filtered)), scannf = FALSE)
niche_params <- niche.param(omi_result)

# ============================================================================
# I. SSI (Species Specialization Index)
# ============================================================================

tolerance <- niche_params[, "Tol"]
tolerance[tolerance == 0] <- min(tolerance[tolerance > 0]) / 10

ssi_raw        <- 1 / tolerance
ssi_normalized <- (ssi_raw - min(ssi_raw)) / (max(ssi_raw) - min(ssi_raw))

ssi_df <- data.frame(
  Species   = rownames(niche_params),
  OMI       = niche_params[, "OMI"],
  Tolerance = niche_params[, "Tol"],
  SSI_raw   = ssi_raw,
  SSI       = ssi_normalized)
ssi_df <- ssi_df[order(-ssi_df$SSI), ]

# ============================================================================
# J. CSI (Community Specialization Index)
# ============================================================================

ssi_vector <- setNames(ssi_df$SSI, ssi_df$Species)

csi_results <- data.frame(station = rownames(mat_comm),
                          CSI = NA, richness = NA, n_species_with_ssi = NA)

for (i in 1:nrow(mat_comm)) {
  station_abd     <- mat_comm[i, ]
  sp_present      <- names(station_abd)[station_abd > 0]
  sp_with_ssi     <- sp_present[sp_present %in% names(ssi_vector)]
  csi_results$richness[i]          <- length(sp_present)
  csi_results$n_species_with_ssi[i] <- length(sp_with_ssi)
  if (length(sp_with_ssi) > 0) {
    abd <- station_abd[sp_with_ssi]
    ssi <- ssi_vector[sp_with_ssi]
    csi_results$CSI[i] <- sum(abd * ssi) / sum(abd)
  }
}

csi_results$region     <- env_aligned$region
csi_results$altitude_m <- env_aligned$altitude_m

# ============================================================================
# K. Statistical tests on CSI
# ============================================================================

csi_by_region <- csi_results %>%
  group_by(region) %>%
  summarise(n = n(),
    CSI_mean = mean(CSI, na.rm = TRUE), CSI_sd = sd(CSI, na.rm = TRUE),
    CSI_min = min(CSI, na.rm = TRUE), CSI_max = max(CSI, na.rm = TRUE),
    richness_mean = mean(richness, na.rm = TRUE), .groups = "drop")

test_mw  <- wilcox.test(CSI ~ region, data = csi_results)
cor_alt  <- cor.test(csi_results$CSI, csi_results$altitude_m, method = "spearman")
cor_rich <- cor.test(csi_results$CSI, csi_results$richness, method = "spearman")

# Deconfounding CSI–richness
model_csi_s            <- lm(CSI ~ richness, data = csi_results)
csi_results$CSI_residuals <- residuals(model_csi_s)
test_resid             <- wilcox.test(CSI_residuals ~ region, data = csi_results)

resid_summary <- csi_results %>%
  group_by(region) %>%
  summarise(resid_mean = round(mean(CSI_residuals), 5),
            resid_sd = round(sd(CSI_residuals), 5), .groups = "drop")

# ============================================================================
# L. CSI figures
# ============================================================================

output_dir <- paste0(data_path, "resultats_CSI/")
dir.create(output_dir, showWarnings = FALSE)

# SSI distribution
p1 <- ggplot(ssi_df, aes(x = SSI)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_vline(xintercept = mean(ssi_df$SSI), linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Distribution of the Species Specialization Index (SSI)",
       subtitle = paste("n =", nrow(ssi_df), "species | Mean =", round(mean(ssi_df$SSI), 3)),
       x = "SSI (0 = generalist, 1 = specialist)", y = "Number of species") +
  theme_minimal(base_size = 12)
ggsave(paste0(output_dir, "Fig1_SSI_distribution.png"), p1, width = 8, height = 6, dpi = 300)

# CSI boxplot by region
p2 <- ggplot(csi_results, aes(x = region, y = CSI, fill = region)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  labs(title = "Community Specialization Index by bioclimatic zone",
       subtitle = paste("Mann-Whitney: W =", test_mw$statistic,
                        ", p =", round(test_mw$p.value, 4)),
       x = "Bioclimatic zone", y = "CSI") +
  scale_fill_manual(values = c("coastal" = "#2166AC", "montane" = "#B2182B")) +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(paste0(output_dir, "Fig2_CSI_boxplot_region.png"), p2, width = 7, height = 6, dpi = 300)

# CSI vs Altitude
p3 <- ggplot(csi_results, aes(x = altitude_m, y = CSI, color = region)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", alpha = 0.3) +
  labs(title = "CSI–Altitude relationship",
       subtitle = paste("Spearman rho =", round(cor_alt$estimate, 3),
                        ", p =", round(cor_alt$p.value, 4)),
       x = "Altitude (m)", y = "CSI", color = "Zone") +
  scale_color_manual(values = c("coastal" = "#2166AC", "montane" = "#B2182B")) +
  theme_minimal(base_size = 12)
ggsave(paste0(output_dir, "Fig3_CSI_altitude.png"), p3, width = 8, height = 6, dpi = 300)

# CSI vs Richness
p4 <- ggplot(csi_results, aes(x = richness, y = CSI, color = region)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", alpha = 0.3) +
  labs(title = "CSI–Species richness relationship",
       subtitle = paste("Spearman rho =", round(cor_rich$estimate, 3),
                        ", p =", round(cor_rich$p.value, 4)),
       x = "Species richness", y = "CSI", color = "Zone") +
  scale_color_manual(values = c("coastal" = "#2166AC", "montane" = "#B2182B")) +
  theme_minimal(base_size = 12)
ggsave(paste0(output_dir, "Fig4_CSI_richesse.png"), p4, width = 8, height = 6, dpi = 300)

# CSI barplot by station
p5 <- ggplot(csi_results, aes(x = reorder(station, CSI), y = CSI, fill = region)) +
  geom_bar(stat = "identity", alpha = 0.8) + coord_flip() +
  labs(title = "Community Specialization Index by station",
       x = "Station", y = "CSI", fill = "Zone") +
  scale_fill_manual(values = c("coastal" = "#2166AC", "montane" = "#B2182B")) +
  theme_minimal(base_size = 12)
ggsave(paste0(output_dir, "Fig5_CSI_stations.png"), p5, width = 10, height = 7, dpi = 300)

# ============================================================================
# M. Export CSI results
# ============================================================================

write_xlsx(ssi_df, paste0(output_dir, "SSI_species.xlsx"))
write_xlsx(csi_results, paste0(output_dir, "CSI_stations.xlsx"))
write_xlsx(as.data.frame(csi_by_region), paste0(output_dir, "CSI_summary_regions.xlsx"))

save(ssi_df, csi_results, omi_result, niche_params, mat_comm, env_aligned,
     file = paste0(output_dir, "CSI_analysis.RData"))

# ============================================================================
# N. Print results summary
# ============================================================================

cat("\n====== DATA ======\n")
n_montaneane <- sum(rownames(mat) %in% stations_montaneane)
cat("Stations:", n_stations, "(", n_stations - n_montaneane, "coastal,",
    n_montaneane, "montaneane )\n")
cat("Species:", n_species, "\n")
cat("Mean richness:", round(mean(richesse), 1), "\n")

cat("\n====== SENSITIVITY ANALYSIS ======\n")
print(sensitivity_table)

cat("\n====== CLASSIFICATION ======\n")
cat("Optimal k:", k_optimal, "| Silhouette:", round(sil_optimal, 3), "\n")
cat("Cluster-zone correspondence:", round(pct_corr, 1), "%\n")
print(tab_corr)

cat("\n====== SILHOUETTE (k = 2 to 5) ======\n")
print(data.frame(k = 2:5, silhouette = round(sil_scores[2:5], 3)))

cat("\n====== DCA ======\n")
cat("Gradient length DCA1:", round(gradient_length, 2), "SD\n")
print(dca_result)

cat("\n====== NMDS ======\n")
cat("Stress:", round(nmds_result$stress, 3), "\n")

cat("\n====== CCA ======\n")
cat("Constrained variance:", round(cca_var_explained, 1), "% | p =", cca_p_global, "\n")
cat("Eigenvalues:", round(cca_eig, 4), "\n")
cat("VIF:", round(vif_values, 2), "\n")
cat("Correlation altitude-zone:", round(cor_alt_zone, 3), "\n")
cat("\nMarginal tests:\n")
print(cca_anova_margin)

cat("\n====== ANOSIM ======\n")
cat("R =", round(anosim_result$statistic, 3), "| p =", round(anosim_result$signif, 4), "\n")

cat("\n====== PERMANOVA ======\n")
print(permanova_result)
cat("PERMANOVA R² =", round(permanova_result$R2[1], 4),
    "| F =", round(permanova_result$F[1], 4),
    "| p =", round(permanova_result$`Pr(>F)`[1], 4), "\n")

cat("\n====== BETADISPER ======\n")
cat("F =", round(betadisper_F, 2), "| p =", round(betadisper_p, 3), "\n")
cat("Distance to centroid:\n")
print(round(betadisper_result$group.distances, 4))
n_coastal  <- n_stations - length(stations_montaneane)
n_montane  <- length(stations_montaneane)
n_configs  <- choose(n_stations, n_montane)
cat("\n--- Permutation space diagnostic ---\n")
cat("Design:", n_coastal, "coastal,", n_montane, "montane\n")
cat("Unique group-label configurations: C(", n_stations, ",", n_montane, ") =",
    n_configs, "\n")
cat("Theoretical minimum p (exhaustive): 1/", n_configs, " = ",
    round(1/n_configs, 4), "\n")
cat("Observed permutation p (9999 perms): PERMANOVA =",
    round(permanova_result$`Pr(>F)`[1], 4),
    "| ANOSIM =", round(anosim_result$signif, 4), "\n")
cat("Note: observed p < theoretical min because random permutations\n",
    "     do not exhaustively enumerate all configurations.\n")

cat("\n====== INDICATOR SPECIES ======\n")
cat("Total (p < 0.05):", n_indic_total, "| Coastal:", n_indic_coastal,
    "| Montaneane:", n_indic_montaneane, "\n")

cat("\n====== ALPHA DIVERSITY ======\n")
print(div_alpha)
cat("\nSummary by zone:\n")
print(div_summary)

cat("\n====== BETA DIVERSITY ======\n")
cat("Gamma:", gamma_div, "| Alpha mean:", round(alpha_mean, 1),
    "| Whittaker:", round(beta_whittaker, 2), "\n")
cat("Sorensen:", round(beta_sor, 3), "| Turnover:", round(turnover_pct, 1),
    "% | Nestedness:", round(nestedness_pct, 1), "%\n")

cat("\n====== LCBD ======\n")
print(lcbd_df)
cat("Montaneane contribution:", round(lcbd_montaneane_total * 100, 1), "%\n")
cat("LCBD ~ Richness (Spearman): rho =", round(cor_lcbd_rich$estimate, 3),
    "| p =", round(cor_lcbd_rich$p.value, 3), "\n")

cat("\n====== RAUP-CRICK ======\n")
print(rc_summary)
if (exists("wilcox_rc")) {
  cat("\nWilcoxon inter-zone vs intra-zone: W =", wilcox_rc$statistic,
      "| p =", format.pval(wilcox_rc$p.value, digits = 4), "\n")
}

cat("\n====== EVENNESS & DOMINANCE ======\n")
print(diversity_indices)
cat("\nSummary by zone:\n")
print(diversity_summary)
cat("\nWilcoxon tests:\n")
print(test_results)

cat("\n====== CSI ======\n")
cat("CSI by region:\n")
print(as.data.frame(csi_by_region))
cat("\nMann-Whitney: W =", test_mw$statistic, "| p =", round(test_mw$p.value, 4), "\n")
cat("CSI ~ Altitude: rho =", round(cor_alt$estimate, 3), "| p =", round(cor_alt$p.value, 4), "\n")
cat("CSI ~ Richness: rho =", round(cor_rich$estimate, 3), "| p =", round(cor_rich$p.value, 4), "\n")
cat("Deconfounding (residuals ~ region): W =", test_resid$statistic,
    "| p =", round(test_resid$p.value, 4), "\n")
print(as.data.frame(resid_summary))

cat("\n====== DONE ======\n")
