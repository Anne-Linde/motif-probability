### This script was used to perform the statistical analyses corresponding to the following article:
## Unveiling hidden temporal physical behavior patterns: 
# A forward algorithm of the hidden semi-Markov model to capture motif probabilty beyond total volume and complexity
## Annelinde Lettink, Mai JM Chin A Paw, Eva Corpeleijn, Teatske M Altenburg, & Wessel N van Wieringen

### User input:
dataDir <- "/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/analyses"
load(paste0(dataDir, "/merged_data.RData")) # Load the data file (obtained by running script: 2_motif_probability_create_dataset.R)

# load packages
library(factoextra, cluster)
library(rstatix)
library(rags2ridges)
library(vioplot)
library(robust)
library(MASS)
library(robustbase)
library(effsize)
library(xtable)
library(devtools)

################################################################################
# Structure data
#################################################################################
# Probabilities
pBGt <- merged_data[,2:11]
pBG <- numeric()
for (u in 1:ncol(pBGt)){
  pBG <- rbind(pBG, as.numeric(pBGt[,u]))
}
pBG <- t(pBG)
pBG <- -log(-log(pBG + sort(unique(as.numeric(pBG)))[2]))
# Re-order motif probabilities: PA 1 min - PA 30 min -- SB 5 min - SB 60 min
pBG <- pBG[, c(6:1, 7:10)] 
motifs <- c("motif.PA_1min", "motif.PA_5min", "motif.PA_7min", "motif.PA_10min", "motif.PA_15min",  "motif.PA_30min", "motif.SB_5min", "motif.SB_10min", "motif.SB_30min", "motif.SB_60min")

Xmot <- pBG

# Volume-based estimates
Xvol <- merged_data[, c("avg_acc_mg", "avg_acc_VMcts", "M30", "M60", "SB", "LPA", "MVPA")]
# Complexity metrics
Xcom <- merged_data[, c("sample_entropy", "lempel_ziv")]

## Standardize covariates
# Centering around mean
Xmot <- sweep(Xmot, 2, FUN="-", apply(Xmot, 2, mean))
Xvol <- sweep(Xvol, 2, FUN = "-", apply(Xvol, 2, mean, na.rm = TRUE))
Xcom <- sweep(Xcom, 2, FUN="-", apply(Xcom, 2, mean))
# Standardized covariate information, per category: motif prob, volume, complexity
Xmot <- sweep(Xmot, 2, FUN="/", apply(Xmot,  2, mad))
Xvol <- sweep(Xvol, 2, FUN="/", apply(Xvol, 2, mad, na.rm = TRUE))
Xcom <- sweep(Xcom, 2, FUN="/", apply(Xcom, 2, mad))

# Response variables: sex (boys2girls) & BMI z-score 
Ybg   <- 1*(merged_data$sex == "1") # 1 = boy, 0 = girl
Yzbmi <- as.numeric(merged_data$zbmi)

## Combine data
data <- cbind(Xmot, Xvol, Xcom, Yzbmi, Ybg)
colnames(data) <- c(motifs, names(data)[11:21])

X <- data[,-c(20:21)]

################################################################################
# Consensus clustering of similar physical behavior patterns (individuals)
################################################################################
devtools::source_url("https://raw.githubusercontent.com/Anne-Linde/motif-probability/refs/heads/test_functions/R/ConsensusClusterPlus.R")
if(!dir.exists(paste0(dataDir, "/consensus_clustering_individuals"))){
  dir.create(paste0(dataDir, "/consensus_clustering_individuals"))
}
dist_matrix <- pairwise_dist(t(data))
if(!file.exists(paste0(dataDir, "/consensus_clustering_individuals/cluster_results.RData"))){ # only run if not performed before as it is computationally expensive

    cluster_kmeans <- ConsensusClusterPlus(
    d=dist_matrix, maxK = 10, reps=1000, pItem=0.8, pFeature = 1, clusterAlg="km",title="consensus_clustering_individuals",
    distance="euclidean", seed = 12345, plot = "png", verbose = TRUE)
  save(cluster_kmeans, file = paste0(dataDir, "/consensus_clustering_individuals/cluster_results.RData"))
} else{
  load(paste0(dataDir, "/consensus_clustering_individuals/cluster_results.RData"))
}

# Cluster analyses statistics
icl = calcICL(cluster_kmeans,title="consensus_clustering_individuals",plot="png",writeTable=TRUE)
inter_cluster_consensus <- aggregate(clusterConsensus ~ k, data = icl$clusterConsensus, mean)
intra_cluster_consensus <- aggregate(itemConsensus ~ k, data = icl$itemConsensus, mean)
pac <- aggregate(itemConsensus ~ k, data = icl$itemConsensus, function(x) mean(x > 0.05 & x < 0.95))

### Check optimal number of clusters
# Perform k-means clustering
png(paste0(dataDir, "/consensus_clustering_individuals/kmeans_Elbow.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_nbclust(dist_matrix, kmeans, method = "wss") # Elbow method
dev.off()
png(paste0(dataDir, "/consensus_clustering_individuals/kmeans_Silhouette.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_nbclust(dist_matrix, kmeans, method = "silhouette") # Silhouette method
dev.off()
png(paste0(dataDir, "/consensus_clustering_individuals/kmeans_Gap.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_nbclust(dist_matrix, kmeans, method = "gap_stat") # Gap statistic
dev.off()
# Perform hierarchical clustering
hc <- hclust(as.dist(dist_matrix), method = "ward.D2")
png(paste0(dataDir, "/consensus_clustering_individuals/hierarchical_dendrogram.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_dend(hc, k = 5, rect = TRUE, show_labels = FALSE) + guides(color = "none", fill = "none")
dev.off()

# 5 clusters is optimal according most methods - Assign cluster labels for 5 clusters
cluster_labels <- cluster_kmeans[[5]]$consensusClass
data_clusters <- cbind(data, cluster_labels)

### Principal Components Analysis (PCA)
################################################################################
if(!dir.exists(paste0(dataDir, "/consensus_clustering_individuals/PCA"))){
  dir.create(paste0(dataDir, "/consensus_clustering_individuals/PCA"))
}
data_clean <- na.omit(data_clusters)

# Perform PCA
pca_result <- prcomp(data_clean[,1:21], scale. = TRUE)
summary(pca_result)

# Plot eigenvalues of the components -- 4 components
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/Eigenvalues_components.png"), width = 170, height = 170, units = "mm", res = 300)
screeplot(pca_result, col = "blue", type = "lines", main = "") 
dev.off()

# Plot Cluster components
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_components12.png"), width = 170, height = 170, units = "mm", res = 300)
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/Figure4.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(1, 2), 
             xlab = paste("Component 1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()
dev.off()

png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_components13.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(1, 3), 
             xlab = paste("Component 1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 3 (", round(summary(pca_result)$importance[2, 3] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()
dev.off()

png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_components14.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(1, 4), 
             xlab = paste("Component 1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 4 (", round(summary(pca_result)$importance[2, 4] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()
dev.off()

png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_components23.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(2, 3), 
             xlab = paste("Component 2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 3 (", round(summary(pca_result)$importance[2, 3] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()
dev.off()

png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_components24.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(2, 4), 
             xlab = paste("Component 2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 4 (", round(summary(pca_result)$importance[2, 4] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()
dev.off()

png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_components34.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(3, 4), 
             xlab = paste("Component 3 (", round(summary(pca_result)$importance[2, 3] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 4 (", round(summary(pca_result)$importance[2, 4] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()
dev.off()

# Plot component loadings
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_component_loadings_12.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(1,2))
dev.off()
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_component_loadings_13.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(1,3))
dev.off()
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_component_loadings_14.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(1,4))
dev.off()
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_component_loadings_23.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(2,3))
dev.off()
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_component_loadings_24.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(2,4))
dev.off()
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/clusters_component_loadings_34.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(3,4))
dev.off()

# Visualize the contribution of the variables contributing to these dimensions (Figure 5)
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/component1_contributions.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_contrib(pca_result, choice = "var", axes = 1) + ggtitle("")  + theme_minimal() +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Physical behavior estimate") # Dimension 1 (x-axis)
dev.off()
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/component2_contributions.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_contrib(pca_result, choice = "var", axes = 2)+ ggtitle("")  + theme_minimal()+   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Physical behavior estimate")# Dimension 2 (y-axis)
dev.off()
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/component3_contributions.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_contrib(pca_result, choice = "var", axes = 3)+ ggtitle("")  + theme_minimal()+   theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("Physical behavior estimate")# Dimension 3
dev.off()
png(paste0(dataDir, "/consensus_clustering_individuals/PCA/component4_contributions.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_contrib(pca_result, choice = "var", axes = 4)+ ggtitle("")  + theme_minimal()+   theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("Physical behavior estimate")# Dimension 4
dev.off()

### Cluster differences
################################################################################
# Calculate the mean of each variable for each cluster
cluster_summary <- aggregate(data_clusters, by = list(cluster_labels), function(x) mean(x, na.rm = TRUE))
colnames(cluster_summary) <- c("cluster", names(data))

table(data$Ybg, cluster_labels) # 1 = boys, 0 = girls

# Heatmap (Figure 6)
data_matrix <- as.matrix(cluster_summary[,-1])
rownames(data_matrix) <- cluster_summary$cluster
x <- rep(0, 5)
y <- rep(0, 5)
z <- rep(0, 5)
xy <- rep(0, 5)
matrix <- cbind(as.matrix(cluster_summary[,2:7]), x, as.matrix(cluster_summary[,8:11]), y, 
                as.matrix(cluster_summary[,12:18]), z, as.matrix(cluster_summary[,19:20]), xy, as.matrix(cluster_summary[,21:22]))
colnames(matrix)[c(7,12,20, 23)] <- c("", " ", "  ", "   ") # Remove labels for x, y, z, xy

png(paste0(dataDir, "/consensus_clustering_individuals/Figure6.png"), width = 170, height = 170, units = "mm", res = 300)
edgeHeat(t(matrix)) + theme_minimal()
dev.off()

# Test for differences in variables between clusters using ANOVA 
for (var in colnames(data)) {
  print(var)
  aov_result <- aov(data[[var]] ~ as.factor(cluster_labels))
  print(summary(aov_result))
  print(eta_squared(aov_result))
  if(summary(aov_result)[[1]]$`Pr(>F)`[1] < .05){
    tukey_result <- TukeyHSD(aov_result, data = data)
    print(tukey_result)
    # Calculate Cohen's d for each pair of clusters
    for (i in 1:(length(unique(cluster_labels)) - 1)) {
      for (j in (i + 1):length(unique(cluster_labels))) {
        group1 <- data[[var]][cluster_labels == unique(cluster_labels)[i]]
        group2 <- data[[var]][cluster_labels == unique(cluster_labels)[j]]
        cohens_d <- cohen.d(group1, group2, na.rm=TRUE)
        print(paste("Cohen's d between cluster", unique(cluster_labels)[i], "and cluster", unique(cluster_labels)[j], ":", cohens_d$estimate))
      }
    }
  }
}

summary(glm(as.factor(Ybg) ~ as.factor(cluster_labels), data = data_clusters, family = "binomial")) #logistic regression differences boys girls
anova(glm(as.factor(Ybg) ~ as.factor(cluster_labels), data = data_clusters, family = "binomial"), test="Chisq")

################################################################################
# Correlation between physical behavior estimates
################################################################################
if(!dir.exists(paste0(dataDir, "/correlation_analyses"))){
  dir.create(paste0(dataDir, "/correlation_analyses"))
}
# Compute correlation matrix
Xall <- cbind(Xmot[,1:6], rnorm(nrow(Xmot), sd=0.000001), Xmot[,7:10], rnorm(nrow(Xmot), sd=0.000001), 
              Xvol[,c(1,2,3,4,6,7,5)], rnorm(nrow(Xmot), sd=0.000001), Xcom)
colnames(Xall) <- c(motifs[1:6], " ", motifs[7:10], "  ", names(Xvol[,c(1,2,3,4,6,7,5)]), "   ", names(Xcom))

corMatAll <- cor(Xall, m="k", use = "pairwise.complete.obs")
corMatAll[c(7,12,20),] <- 0
corMatAll[,c(7,12,20)] <- 0

# Figure7
png(paste0(dataDir, "/correlation_analyses/Figure7.png"), width = 170, height = 170, units = "mm", res = 300)
edgeHeat(corMatAll)
dev.off()

# Mask the upper triangle of the matrix
corMatAll[upper.tri(corMatAll)] <- NA
corMatAll[c(7,12,20),] <- NA
corMatAll[,c(7,12,20)] <- NA

# Convert to data frame for nicer LaTeX output
cor_matrix_df <- as.data.frame(corMatAll)
rownames(cor_matrix_df) <- colnames(corMatAll)

xtable_obj <- xtable(cor_matrix_df, caption = "Correlation Matrix", label = "tab:lower_triangle") # Use xtable to generate LaTeX table
print(xtable_obj, type = "latex", booktabs = TRUE, na.string = "") # Print the LaTeX code


################################################################################
# Consensus clustering of physical behavior estimates
################################################################################
if(!dir.exists(paste0(dataDir, "/consensus_clustering_estimates"))){
  dir.create(paste0(dataDir, "/consensus_clustering_estimates"))
}
dist_matrix_vars <- pairwise_dist(t(data[,1:19]))
if(!file.exists(paste0(dataDir, "/consensus_clustering_estimates/cluster_results.RData"))){ # only run if not performed before as it is computationally expensive
  cluster_vars <- ConsensusClusterPlus(
    d=dist_matrix_vars, maxK = 10, reps=1000, pItem=0.8, pFeature = 1, clusterAlg="km",title="consensus_clustering_vars",
    distance="euclidean", seed = 12345, plot = "png", verbose = TRUE)
  save(cluster_vars, file = paste0(dataDir, "/consensus_clustering_estimates/cluster_estimates.RData"))
} else{
  load(paste0(dataDir, "/consensus_clustering_estimates/cluster_estimates.RData"))
}

# Cluster analyses statistics
icl_var = calcICL(cluster_vars,title="consensus_clustering_estimates",plot="png",writeTable=TRUE)
inter_cluster_consensus_var <- aggregate(clusterConsensus ~ k, data = icl_var$clusterConsensus, mean)
intra_cluster_consensus_var <- aggregate(itemConsensus ~ k, data = icl_var$itemConsensus, mean)
pac_var <- aggregate(itemConsensus ~ k, data = icl_var$itemConsensus, function(x) mean(x > 0.05 & x < 0.95))

# Perform k-means clustering
png(paste0(dataDir, "/consensus_clustering_estimates/kmeans_Elbow.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_nbclust(dist_matrix_vars, kmeans, method = "wss") # Elbow method
dev.off()
png(paste0(dataDir, "/consensus_clustering_estimates/kmeans_Silhouette.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_nbclust(dist_matrix_vars, kmeans, method = "silhouette") # Silhouette method
dev.off()
png(paste0(dataDir, "/consensus_clustering_estimates/kmeans_Gap.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_nbclust(dist_matrix_vars, kmeans, method = "gap_stat") # Gap statistic
dev.off()

# Perform hierarchical clustering
hc <- hclust(as.dist(dist_matrix_vars), method = "ward.D2")
png(paste0(dataDir, "/consensus_clustering_estimates/hierarchical_dendrogram.png"), width = 170, height = 170, units = "mm", res = 300)
fviz_dend(hc, k = 6, rect = TRUE, show_labels = FALSE)
dev.off()

# 6 clusters is optimal according to both methods -- Assign cluster labels
cluster_labels_vars <- cluster_vars[[6]]$consensusClass

################################################################################
# regression analysis of sex (boys vs girls)
################################################################################
if(!dir.exists(paste0(dataDir, "/regression_analyses"))){
  dir.create(paste0(dataDir, "/regression_analyses"))
}
pseudo_r_squared <- function(model) {
  1 - sum(resid(model)^2, na.rm = TRUE) / sum((model$model[,1] - mean(model$model[,1], na.rm = TRUE))^2, na.rm = TRUE)
}

glmRes_sex <- numeric()
for (u in 1:ncol(X)){
  model <- glmRob(Ybg ~ X[,u], family="binomial")
  summary_robust <- summary(model)
  estimate <- summary_robust$coefficients[2, "Estimate"]
  p_value <- summary(model)$coefficients[2, "Pr(>|z|)"]
  r_squared <- pseudo_r_squared(model)
  logLik_value <- sum(model$y * log(model$fitted.values) + (1 - model$y) * log(1 - model$fitted.values))
  aic_value <- 2 * length(coef(model)) - 2 * logLik_value
  bic_value <- log(length(model$y)) * length(coef(model)) - 2 * logLik_value
  mse <- mean(resid(model)^2)
  glmRes_sex <- rbind(glmRes_sex, c(cor(Yzbmi, X[,u], m="k", use = "pairwise.complete.obs"), estimate, p_value, r_squared, aic_value, bic_value, mse))
}

colnames(glmRes_sex) <- c("cor(B/G,Covariate)", "Estimate", "p-value", "r_squared", "AIC", "BIC", "MSE")
rownames(glmRes_sex) <- names(X)

# Convert the data frame to an xtable object
table_sex <- xtable(glmRes_sex, caption = "Simple linear regression results with sex as dependent variable", label = "tab:regression_sex")
print(table_sex, type = "latex") # Print the xtable object as LaTeX code

# diagnostic of logistic regression (single probs in model)
nBins <- 9
ps <- matrix(NA, nBins+1, ncol(X))
for (u in 1:ncol(X)){
  qs <- c(min(X[,u], na.rm = TRUE), quantile(X[,u], prob=c(1:nBins)/(nBins+1), na.rm = TRUE), max(X[,u], na.rm = TRUE))
  for (w in 1:(nBins+1)){
    idQs <- intersect(which(X[,u] >= qs[w]), which(X[,u] <= qs[w+1]))	
    ps[w,u] <- sum(Ybg[idQs], na.rm = TRUE) / length(idQs)
  }
}
colnames(ps) <- colnames(X)
ps    <- ps - mean(ps)
psMod <- cbind(ps[,1:6], rep(0, nrow(ps)), ps[,7:10], rep(0, nrow(ps)), ps[,11:17], rep(0, nrow(ps)), ps[,18:19])
psMod <- cbind(ps[,1:6], rnorm(nrow(ps), sd=1/1000), ps[,7:10], rnorm(nrow(ps), sd=1/1000), ps[,11:17], rnorm(nrow(ps), sd=1/1000), ps[,18:19])
colnames(psMod)[c(7,12,20)] <- c("", " ", "  ")
rownames(psMod) <- c("0-10th percentile", "10-20th percentile", "20-30th percentile", 
                     "30-40th percentile", "40-50th percentile", "50-60th percentile", 
                     "60-70th percentile", "70-80th percentile", "80-90th percentile", 
                     "90-100th percentile")
png(paste0(dataDir, "/regression_analyses/Figure8.png"), width = 170, height = 170, units = "mm", res = 300)
edgeHeat(as.matrix(psMod)*-1)
dev.off()

# # Fit the robust linear models
df_Xvolcomp <- cbind(Xvol, Xcom, Ybg)
df_Xvolcompmotifs <- cbind(Xvol, Xcom, Xmot, Ybg)
colnames(df_Xvolcompmotifs) <- c(names(Xvol), names(Xcom), motifs, "Ybg")

baseline_model <- glmRob(Ybg ~ ., data =  df_Xvolcomp, family="binomial")
full_model <-glmRob(Ybg ~ ., data =  df_Xvolcompmotifs, family="binomial")

logLik_value_baseline <- sum(baseline_model$y * log(baseline_model$fitted.values) + (1 - baseline_model$y) * log(1 - baseline_model$fitted.values))
logLik_value_full <- sum(full_model$y * log(full_model$fitted.values) + (1 - full_model$y) * log(1 - full_model$fitted.values))

# Perform the Likelihood Ratio Test
lrt_stat <- -2 * (logLik_value_baseline - logLik_value_full)
df <- length(coef(full_model)) - length(coef(baseline_model))
p_value <- pchisq(lrt_stat, df, lower.tail = FALSE)

aic_value_baseline <- 2 * length(coef(baseline_model)) - 2 * logLik_value_baseline
aic_value_full <- 2 * length(coef(full_model)) - 2 * logLik_value_full

bic_value_baseline <- log(length(baseline_model$y)) * length(coef(baseline_model)) - 2 * logLik_value_baseline
bic_value_full <- log(length(full_model$y)) * length(coef(full_model)) - 2 * logLik_value_full

r_squared_baseline <- pseudo_r_squared(baseline_model)
r_squared_full <- pseudo_r_squared(full_model)

# Deviance comparison
deviance_diff <-model_VolComp$deviance - model_VolCompMotifs$deviance
df_diff <- length(coef(model_VolCompMotifs)) - length(coef(model_VolComp))
p_value <- pchisq(deviance_diff, df = df_diff, lower.tail = FALSE)

################################################################################
# regression analysis of BMI z-score
################################################################################

glmRes_zbmi <- numeric()
for (u in 1:ncol(X)){
  model <- rlm(Yzbmi ~ X[,u])
  summary_robust <- summary(model)
  estimate <- summary_robust$coefficients[2, "Value"]
  p_value <- 2 * pt(-abs(summary_robust$coefficients[2, "t value"]), summary_robust$df[2])
  r_squared <- pseudo_r_squared(model)
  aic_value <- AIC(model)
  bic_value <- BIC(model)
  mse <- mean(resid(model)^2)
  
  glmRes_zbmi <- rbind(glmRes_zbmi, c(cor(Yzbmi, X[,u], m="k", use = "pairwise.complete.obs"), estimate, p_value, r_squared, aic_value, bic_value, mse))
}

colnames(glmRes_zbmi) <- c("cor(zbmi,Covariate)", "Estimate", "p-value", "r_squared", "AIC", "BIC", "MSE")
rownames(glmRes_zbmi) <- names(X)

# Convert the data frame to an xtable object
table_bmi <- xtable(glmRes_zbmi, caption = "Simple linear regression results with BMI z-score as dependent variable", label = "tab:regression_bmi")
print(table_bmi, type = "latex") # Print the xtable object as LaTeX code

