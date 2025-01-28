# Tabula rasa
rm(list=ls())
gc()

### User input required:
dataDir <- "/Users/annelindelettink/GECKO/preprocessing/manuscript/subset/analyses"
#load(paste0(dataDir, "/merged_data.RData"))

load(paste0(dataDir, "/merged_data_rerun.RData"))

#setwd("/home/wessel/Research/AnnelindeMotif/FW_ Motif probability vorderingen")
#load("merged_data.RData")

# load packages
library(factoextra, cluster)
library(rstatix)
library(rags2ridges)
library(vioplot)
library(robust)
library(MASS)
library(robustbase)
library(effsize)

rDir <- "/Users/annelindelettink/Documents/Work MacBook Pro Annelinde/My Little Moves (MLM)/Sequence mapping/Physical behavior patterns/motif-probability/R/"
setwd(rDir)
source("ConsensusClusterPlus.R")
setwd(dataDir)

################################################################################
# structure data
#################################################################################
# Probabilities
pBGt <- merged_data[,2:11]
pBG <- numeric()
for (u in 1:ncol(pBGt)){
  pBG <- rbind(pBG, as.numeric(pBGt[,u]))
}
pBG <- t(pBG)
pBG <- -log(-log(pBG + sort(unique(as.numeric(pBG)))[2]))
# Restructure variables PA 1 min - PA 30 min -- SB 5 min - SB 60 min
pBG <- pBG[, c(6:1, 7:10)] 
motifs <- c("motif.PA_1min", "motif.PA_5min", "motif.PA_7min", "motif.PA_10min", "motif.PA_15min",  "motif.PA_30min", "motif.SB_5min", "motif.SB_10min", "motif.SB_30min", "motif.SB_60min")

Xmot <- pBG # Motif probabilities
Xvol <- merged_data[, c("avg_acc_mg", "avg_acc_VMcts", "M30", "M60", "SB", "LPA", "MVPA")] # Volume metrics

#Xvol <- merged_data[, c("avg_acc_mg", "avg_acc_VMcts", "L5", "M30", "M60", "SB", "LPA", "MVPA")] # Volume metrics
#Xvol <- merged_data[, c("mean_HFEN+", "L5_HFEN+", "M5_HFEN+", "mean_NeishabouriCount_VM", "MVPA_T890_NeishabouriCount_y")] # Volume metrics
Xcom <- merged_data[, c("sample_entropy", "lempel_ziv")] # Complexity scores

# Centering around mean
Xmot <- sweep(Xmot, 2, FUN="-", apply(Xmot, 2, mean))
Xvol <- sweep(Xvol, 2, FUN = "-", apply(Xvol, 2, mean, na.rm = TRUE))
#Xvol <- sweep(Xvol, 2, FUN="-", apply(Xvol, 2, mean))
Xcom <- sweep(Xcom, 2, FUN="-", apply(Xcom, 2, mean))

# Standardized covariate information, per category: motif prob, volume, complexity
Xmot <- sweep(Xmot, 2, FUN="/", apply(Xmot,  2, mad))
#Xvol <- sweep(Xvol, 2, FUN="/", apply(Xvol, 2, mad))
Xvol <- sweep(Xvol, 2, FUN="/", apply(Xvol, 2, mad, na.rm = TRUE))
Xcom <- sweep(Xcom, 2, FUN="/", apply(Xcom, 2, mad))

# Response type: boys2girls & zbmi, added 
Ybg   <- 1*(merged_data$sex == "1") # 1 = boy, 0 = girl
Yzbmi <- as.numeric(merged_data$zbmi)

data <- cbind(Xmot, Xvol, Xcom, Yzbmi, Ybg)
#colnames(data) <- c(motifs, names(data)[11:19])
colnames(data) <- c(motifs, names(data)[11:21])

#X <- data[,-c(18:19)]
X <- data[,-c(20:21)]


################################################################################
# Consensus clustering of individuals
################################################################################
# Cluster Consensus: average consensus index between all pairs of items belonging to the same cluster, for each K.
# Item consensus: average consensus index between item i and all the (other) items in cluster cl, for all i and cl, for each K.
# Intra- and Inter- Cluster consensus: intra consensus statistic (the mean of all cluster consensus for each K) and inter consensus statistic (mean of all item consensus between an item and all clusters to which the item does not belong, for each K).
# Cluster Consensus Plot: This plot highlights the mean pairwise consensus values between a cluster's members for each k. The color scheme follows all previous graphs and sample are stacked bars grouped by K value on the horizontal x-axis. High values show that the clusters hold high stability and likewise low values highlights a clusters instability. 
# Item Consensus Plot: Each stacked bar is a sample. Item-consensus values are indicated by the heights of the colored portion of the bars (using the tracking color scheme). This plot provides a view of item-consensus across all other clusters at a given k. As Wilkerman (2010) explains, with this plot it is possible to see if a sample is a very "pure" member of a cluster or if it shares high consensus to multiple clusters (large rectangles in a column of multiple colors), suggesting that it is an unstable member.
# The Proportion of Ambiguous Clustering (PAC) is the fraction of sample pairs that hold consensus index values within a given sub-interval (x1, x2) in [0,1] (usually, x1 = 0.1 and x2 = 0.9). The CDF values correspond to the fraction of sample pairs with a consensus index values less or equal to the value 'c'. The PAC is then calculated by CDF(x2) - CDF(x1), optimal K should present a low PAC score. 

# cluster_kmeans <- ConsensusClusterPlus(
#   d=as.matrix(t(data)), maxK = 10, reps=1000, pItem=0.8, pFeature = 1, clusterAlg="km",title="kmeans_consensus_cluster",
#   distance="euclidean", seed = 12345, plot = "png", verbose = TRUE)

# Calculate pairwise complete distance matrix
dist_matrix <- pairwise_dist(data)

cluster_kmeans <- ConsensusClusterPlus(
  d=pairwise_dist(data), maxK = 10, reps=1000, pItem=0.8, pFeature = 1, clusterAlg="km",title="kmeans_consensus_cluster_opnieuw",
  distance="euclidean", seed = 12345, plot = "png", verbose = TRUE)

# save(cluster_kmeans, file = paste0(dataDir, "/kmeans_consensus_cluster/cluster_results.RData"))
save(cluster_kmeans, file = paste0(dataDir, "/rerun/kmeans_consensus_cluster/cluster_results.RData"))

#load(paste0(dataDir, "/kmeans_consensus_cluster/cluster_results.RData"))
setwd(paste0(dataDir, "/rerun"))
icl = calcICL(cluster_kmeans,title="kmeans_consensus_cluster",plot="png",writeTable=TRUE)

# Cluster analyses statistics
inter_cluster_consensus <- aggregate(clusterConsensus ~ k, data = icl$clusterConsensus, mean)
intra_cluster_consensus <- aggregate(itemConsensus ~ k, data = icl$itemConsensus, mean)
pac <- aggregate(itemConsensus ~ k, data = icl$itemConsensus, function(x) mean(x > 0.05 & x < 0.95))

#Check number of clusters
fviz_nbclust(dist_matrix, kmeans, method = "wss") # Elbow method
fviz_nbclust(dist_matrix, kmeans, method = "silhouette") # Silhouette method
fviz_nbclust(dist_matrix, kmeans, method = "gap_stat") # Gap statistic

# Perform hierarchical clustering
hc <- hclust(as.dist(dist_matrix), method = "ward.D2")
fviz_dend(hc, k = 5, rect = TRUE, show_labels = FALSE)
# 5 clusters is optimal according most methods

# Analyze clusters (descriptive statistics and differences)
################################################################################
# Calculate the mean of each variable for each cluster
cluster_labels <- cluster_kmeans[[5]]$consensusClass
data_clusters <- cbind(data, cluster_labels)
cluster_summary <- aggregate(data_clusters, by = list(cluster_labels), function(x) mean(x, na.rm = TRUE))
colnames(cluster_summary) <- c("cluster", names(data))

table(data$Ybg, cluster_labels) # 1 = boys, 0 = girls

# # Heatmap of cluster summary (eerder gebruikte figuur) vervang door code die volgt
# data_matrix <- as.matrix(cluster_summary[,-1])
# rownames(data_matrix) <- cluster_summary$cluster
# pheatmap(data_matrix, 
#          cluster_rows = TRUE,
#          cluster_cols = TRUE, 
#          display_numbers = TRUE,
#          color = colorRampPalette(c("blue", "white", "red"))(50),
#          main = "Heatmap of motif probabilities, volume and complexity metrics per cluster")
# Heatmap of cluster summary
data_matrix <- as.matrix(cluster_summary[,-1])
rownames(data_matrix) <- cluster_summary$cluster

x <- rep(0, 5)
y <- rep(0, 5)
z <- rep(0, 5)
xy <- rep(0, 5)
matrix <- cbind(as.matrix(cluster_summary[,2:7]), x, as.matrix(cluster_summary[,8:11]), y, 
                as.matrix(cluster_summary[,12:18]), z, as.matrix(cluster_summary[,19:20]), xy, as.matrix(cluster_summary[,21:22]))
colnames(matrix)[c(7,12,20, 23)] <- c("", " ", "  ", "   ") # Remove labels for x, y, z, xy
edgeHeat(t(matrix)) + theme_minimal()

## PCA clusters
data_clean <- na.omit(data_clusters)

pca_result <- prcomp(data_clean[,1:21], scale. = TRUE)
summary(pca_result)
screeplot(pca_result, col = "blue", type = "lines") 
# 4 main components

fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(1, 2), 
             xlab = paste("Component 1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()

fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(1, 3), 
             xlab = paste("Component 1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 3 (", round(summary(pca_result)$importance[2, 3] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()

fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(1, 4), 
             xlab = paste("Component 1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 4 (", round(summary(pca_result)$importance[2, 4] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()

fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(2, 3), 
             xlab = paste("Component 2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 3 (", round(summary(pca_result)$importance[2, 3] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()

fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(2, 4), 
             xlab = paste("Component 2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 4 (", round(summary(pca_result)$importance[2, 4] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()

fviz_cluster(list(data = data_clean[,1:21], cluster = data_clean[,22]), 
             geom = "point", 
             ellipse.type = "convex", 
             show.clust.cent = TRUE, 
             axes = c(3, 4), 
             xlab = paste("Component 3 (", round(summary(pca_result)$importance[2, 3] * 100, 1), "%)", sep = ""), 
             ylab = paste("Component 4 (", round(summary(pca_result)$importance[2, 4] * 100, 1), "%)", sep = "")) + 
  ggtitle("") + 
  theme_minimal()


# Inspect component loadings
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(1,2))
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(1,3))
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(1,4))
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(2,3))
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(2,4))
fviz_pca_var(pca_result, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), axes = c(3,4))

# Visualize the contribution of the variables contributing to these dimensions
fviz_contrib(pca_result, choice = "var", axes = 1) + ggtitle("Contributions to Component 1")  + theme_minimal() +   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Physical behavior estimate") # Dimension 1 (x-axis)
fviz_contrib(pca_result, choice = "var", axes = 2)+ ggtitle("Contributions to Component 2")  + theme_minimal()+   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Physical behavior estimate")# Dimension 2 (y-axis)
fviz_contrib(pca_result, choice = "var", axes = 3)+ ggtitle("Contributions to Component 3")  + theme_minimal()+   theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("Physical behavior estimate")# Dimension 3
fviz_contrib(pca_result, choice = "var", axes = 4)+ ggtitle("Contributions to Component 4")  + theme_minimal()+   theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab("Physical behavior estimate")# Dimension 4

# Calculate the importance of each variable in distinguishing clusters
importance <- cluster_summary %>% 
  summarise(across(-cluster, ~var(.))) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Variance") %>%
  arrange(desc(Variance))


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
# Coherence between variables
################################################################################
# Compute correlation matrix
Xall <- cbind(Xmot[,1:6], rnorm(nrow(Xmot), sd=0.000001), Xmot[,7:10], rnorm(nrow(Xmot), sd=0.000001), 
              Xvol[,c(1,2,3,4,6,7,5)], rnorm(nrow(Xmot), sd=0.000001), Xcom)
colnames(Xall) <- c(motifs[1:6], " ", motifs[7:10], "  ", names(Xvol[,c(1,2,3,4,6,7,5)]), "   ", names(Xcom))

corMatAll <- cor(Xall, m="k", use = "pairwise.complete.obs")
corMatAll[c(7,12,20),] <- 0
corMatAll[,c(7,12,20)] <- 0

pdf("heatmap_correlationAmongMetrics.pdf")
png("heatmap_correlationAmongMetrics.png")

edgeHeat(corMatAll)
dev.off()

# Consensus clustering variables

dist_matrix_vars <- pairwise_dist(t(data[,1:19]))

 cluster_vars <- ConsensusClusterPlus(
  d=pairwise_dist(t(data[,1:19])), maxK = 10, reps=1000, pItem=0.8, pFeature = 1, clusterAlg="km",title="consensus_clustering_vars",
   distance="euclidean", seed = 12345, plot = "png", verbose = TRUE)
 save(cluster_vars, file = paste0(dataDir, "rerun/consensus_clustering_vars/cluster_results.RData"))

load(paste0(dataDir, "/consensus_clustering_vars/cluster_results.RData"))
icl_var = calcICL(cluster_vars,title="consensus_clustering_vars",plot="png",writeTable=TRUE)

# Cluster analyses statistics
inter_cluster_consensus_var <- aggregate(clusterConsensus ~ k, data = icl_var$clusterConsensus, mean)
intra_cluster_consensus_var <- aggregate(itemConsensus ~ k, data = icl_var$itemConsensus, mean)
pac_var <- aggregate(itemConsensus ~ k, data = icl_var$itemConsensus, function(x) mean(x > 0.05 & x < 0.95))

#Check number of clusters
fviz_nbclust(dist_matrix_vars, kmeans, method = "wss") # Elbow method
fviz_nbclust(dist_matrix_vars, kmeans, method = "silhouette") # Silhouette method
fviz_nbclust(dist_matrix_vars, kmeans, method = "gap_stat") # Gap statistic

# Perform hierarchical clustering
hc <- hclust(as.dist(dist_matrix_vars), method = "ward.D2")
fviz_dend(hc, k = 6, rect = TRUE, show_labels = FALSE)
# 5 clusters is optimal according to both methods

# Analyze clusters (descriptive statistics and differences)
################################################################################
# Calculate the mean of each variable for each cluster
cluster_labels_vars <- cluster_vars[[6]]$consensusClass

################################################################################
# regression analysis of boys vs girls
################################################################################

pseudo_r_squared <- function(model) {
  1 - sum(resid(model)^2, na.rm = TRUE) / sum((model$model[,1] - mean(model$model[,1], na.rm = TRUE))^2, na.rm = TRUE)
}

# # logistic regression of sex, one covariate at the time
# glmRes <- numeric()
# for (u in 1:ncol(X)){
#   glmRes <- rbind(glmRes, c(cor(Ybg, X[,u], m="k"), coef(summary(glmRob(Ybg ~ X[,u], family="binomial")))[2,c(1,4)]))
# }
# colnames(glmRes)[1] <- c("cor(B/G,Covariate)")
# rownames(glmRes) <- names(X)

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
pdf("heatmap_binBGratio2variable.png")
edgeHeat(as.matrix(psMod))
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

# model_VolCompMotifs <- lmRob(Yzbmi ~ ., data = df_Xvolcompmotifs_clean)
# 
# # Get the summary of the model
# summary_model_vol <- summary(model_Vol)
# summary_model_volcomp <- summary(model_VolComp)
# summary_model_volcompmotifs <- summary(model_VolCompMotifs)
# 
# #Added value
# 
# # Extract estimates and p-values
# estimates <- summary_model_vol$coefficients[, "Value"]
# p_values <- 2 * pt(-abs(summary_model$coefficients[, "t value"]), summary_model$df[2])
# 
# # Print estimates and p-values
# print(estimates)
# print(p_values)

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

# # Add Yzbmi to the DataFrame
# df_Xvol <- cbind(Xvol, Yzbmi)
# df_Xvolcomp <- cbind(Xvol, Xcom, Yzbmi)
# df_Xvolcompmotifs <- cbind(Xvol, Xcom, Xmot, Yzbmi)
# df_Xvol_clean <- na.omit(df_Xvol)
# df_Xvolcomp_clean <- na.omit(df_Xvolcomp)
# df_Xvolcompmotifs_clean <- na.omit(df_Xvolcompmotifs)
# colnames(df_Xvolcompmotifs_clean) <- c(names(Xvol), names(Xcom), motifs, "Yzbmi")
# 
# # Fit the robust linear models
# model_Vol <- lmRob(Yzbmi ~ ., data = Xvol_clean)
# model_VolComp <- lmRob(Yzbmi ~ ., data = df_Xvolcomp_clean)
# model_VolCompMotifs <- lmRob(Yzbmi ~ ., data = df_Xvolcompmotifs_clean)
# 
# # Get the summary of the model
# summary_model_vol <- summary(model_Vol)
# summary_model_volcomp <- summary(model_VolComp)
# summary_model_volcompmotifs <- summary(model_VolCompMotifs)
# 
# #Added value
# 
# # Extract estimates and p-values
# estimates <- summary_model_vol$coefficients[, "Value"]
# p_values <- 2 * pt(-abs(summary_model$coefficients[, "t value"]), summary_model$df[2])
# 
# # Print estimates and p-values
# print(estimates)
# print(p_values)
# 
# # Add Yzbmi to the DataFrame
# dataframe <- cbind(Xvol, Xcom, Yzbmi)
# 
# # Fit the extended linear model using rlm from the MASS package
# dataframe_ext <- cbind(Xvol, Xcom, Xmot, Yzbmi)
# dataframe_ext_clean <- na.omit(dataframe_ext) # Handle NA values by dropping rows with NA values
# names(dataframe_ext_clean) <- c(names(Xvol), names(Xcom), motifs, "Yzbmi")
# 
# model_ext <- rlm(Yzbmi ~ ., data = dataframe_ext_clean)
# 
# # Get the summary of the model
# summary_model <- summary(model_ext)
# 
# # Extract estimates and p-values
# estimates <- summary_model$coefficients[, "Value"]
# p_values <- 2 * pt(-abs(summary_model$coefficients[, "t value"]), summary_model$df[2])


