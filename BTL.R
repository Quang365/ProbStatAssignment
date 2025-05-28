# Cài đặt gói
if (!require(psych)) install.packages("psych")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(factoextra)) install.packages("factoextra")
library(psych)
library(ggplot2)
library(factoextra)

# Data Preprocessing
getwd()
BRCA <- read.csv("BRCA.csv", na.strings = c("", " ", "NA", "missing"), stringsAsFactors = FALSE)
head(BRCA)
BRCA[BRCA == ""] <- NA
BRCA <- BRCA[!apply(is.na(BRCA), 1, all), ]
colSums(is.na(BRCA))  # Kiểm tra NA

# Xử lý NA
BRCA$Age[is.na(BRCA$Age)] <- mean(BRCA$Age, na.rm = TRUE)
BRCA$Protein1[is.na(BRCA$Protein1)] <- mean(BRCA$Protein1, na.rm = TRUE)
BRCA$Protein2[is.na(BRCA$Protein2)] <- mean(BRCA$Protein2, na.rm = TRUE)
BRCA$Protein3[is.na(BRCA$Protein3)] <- mean(BRCA$Protein3, na.rm = TRUE)
BRCA$Protein4[is.na(BRCA$Protein4)] <- mean(BRCA$Protein4, na.rm = TRUE)
BRCA <- BRCA[!is.na(BRCA$Tumour_Stage) & !is.na(BRCA$Patient_Status), ]
colSums(is.na(BRCA))  # Kiểm tra lại NA

# Chuyển đổi biến
BRCA$Tumour_Stage <- as.factor(BRCA$Tumour_Stage)
BRCA$Patient_Status <- as.factor(BRCA$Patient_Status)

# Tạo tập dữ liệu
new_data <- BRCA[, c("Age", "Protein1", "Protein2", "Protein3", "Protein4", 
                     "Tumour_Stage", "Patient_Status")]
head(new_data)
write.csv(new_data, "BRCA_processed.csv", row.names = FALSE)

# Thống kê biến liên tục
desc_stats <- describe(new_data[, c("Age", "Protein1", "Protein2", "Protein3", "Protein4")], fast = TRUE)
head(desc_stats)
write.csv(desc_stats, "descriptive_stats.csv", row.names = TRUE)

# Histogram
vars <- c("Age", "Protein1", "Protein2", "Protein3", "Protein4")
for (v in vars) {
  p <- ggplot(new_data, aes_string(x = v, fill = "Patient_Status")) +
    geom_histogram(binwidth = 1, alpha = 0.6, position = "identity") +
    labs(title = paste("Histogram of", v, "by Patient Status"))
  print(p)
  ggsave(paste0("histogram_", v, ".png"), plot = p, width = 6, height = 4, dpi = 300)
}

# Bar plot
p_bar <- ggplot(new_data, aes(x = Tumour_Stage, fill = Patient_Status)) +
  geom_bar(position = "dodge") +
  labs(title = "Tumour Stage by Patient Status") +
  theme_minimal()
print(p_bar)
ggsave("bar_plot_tumour_stage.png", plot = p_bar, width = 6, height = 4, dpi = 300)

p_bar2 <- ggplot(new_data, aes(x = Patient_Status, fill = Tumour_Stage)) +
  geom_bar(position = "dodge") +
  labs(title = "Patient Status by Tumour Stage") +
  theme_minimal()
print(p_bar2)
ggsave("bar_plot_patient_status.png", plot = p_bar2, width = 6, height = 4, dpi = 300)

# Thống kê biến phân loại
tumour_table <- table(new_data$Tumour_Stage)
status_table <- table(new_data$Patient_Status)
head(tumour_table)
head(status_table)
write.csv(tumour_table, "tumour_stage_freq.csv")
write.csv(status_table, "patient_status_freq.csv")

# Chuẩn bị dữ liệu cho K-means
data_for_clustering <- new_data[, c("Age", "Protein1", "Protein2", "Protein3", "Protein4")]

# Kiểm tra NA và sd
colSums(is.na(data_for_clustering))  # Phải trả về 0
sds <- apply(data_for_clustering, 2, sd)
print("Degree of standard deviation:")
print(sds)
if (any(sds == 0)) {
  cat("Warning: Columns with sd = 0 detected. Removing them.\n")
  data_for_clustering <- data_for_clustering[, sds != 0]
}

# Chuẩn hóa dữ liệu
data_scaled <- scale(data_for_clustering)

# Kiểm tra NA/NaN/Inf trong data_scaled
cat("Check NA in data_scaled:\n")
print(colSums(is.na(data_scaled)))
cat("Check NaN in data_scaled:\n")
print(any(is.nan(data_scaled)))
cat("Check Inf in data_scaled:\n")
print(any(is.infinite(data_scaled)))

# Nếu vẫn có NA/NaN, loại bỏ hàng chứa NA
if (any(is.na(data_scaled)) || any(is.nan(data_scaled))) {
  valid_rows <- complete.cases(data_scaled)
  data_scaled <- data_scaled[valid_rows, ]
  new_data <- new_data[valid_rows, ]
  cat("Removed", sum(!valid_rows), "rows with NA/NaN.\n")
}

# Elbow Method
wss <- function(k) {
  kmeans(data_scaled, centers = k, nstart = 25)$tot.withinss
}
k_values <- 1:10
wss_values <- sapply(k_values, wss)
p_elbow <- fviz_nbclust(data_scaled, kmeans, method = "wss") +
  labs(title = "Elbow Method for K-means Clustering")
print(p_elbow)
ggsave("elbow_plot.png", plot = p_elbow, width = 6, height = 4, dpi = 300)

# K-means với k = 3
set.seed(123)
kmeans_result <- kmeans(data_scaled, centers = 3, nstart = 25)
new_data$Kmeans_Cluster <- as.factor(kmeans_result$cluster)
write.csv(new_data, "BRCA_kmeans_output.csv", row.names = FALSE)  # Lưu dữ liệu với cụm

# Chi-square Test
sink("chisq_results.txt")
print("Chi-square Test: K-means Clusters vs Tumour_Stage")
chisq_tumour <- chisq.test(table(new_data$Kmeans_Cluster, new_data$Tumour_Stage))
print(chisq_tumour)
if (any(chisq_tumour$expected < 5)) {
  print("Expected frequencies < 5, using Fisher’s Exact Test:")
  fisher_tumour <- fisher.test(table(new_data$Kmeans_Cluster, new_data$Tumour_Stage))
  print(fisher_tumour)
}
write.csv(table(new_data$Kmeans_Cluster, new_data$Tumour_Stage), "chisq_tumour_table.csv")

print("Chi-square Test: K-means Clusters vs Patient_Status")
chisq_status <- chisq.test(table(new_data$Kmeans_Cluster, new_data$Patient_Status))
print(chisq_status)
if (any(chisq_status$expected < 5)) {
  print("Expected frequencies < 5, using Fisher’s Exact Test:")
  fisher_status <- fisher.test(table(new_data$Kmeans_Cluster, new_data$Patient_Status))
  print(fisher_status)
}
write.csv(table(new_data$Kmeans_Cluster, new_data$Patient_Status), "chisq_status_table.csv")
sink()

# PCA Plot
pca_result <- prcomp(data_scaled, scale. = FALSE)
summary(pca_result)
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], 
                       Kmeans_Cluster = new_data$Kmeans_Cluster)
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Kmeans_Cluster)) +
  geom_point(size = 3) +
  labs(title = "K-means Clustering Results (PCA)", x = "PC1", y = "PC2") +
  theme_minimal()
print(p_pca)
ggsave("pca_plot.png", plot = p_pca, width = 6, height = 4, dpi = 300)

# Thống kê theo cụm
cluster_means <- aggregate(cbind(Age, Protein1, Protein2, Protein3, Protein4) ~ Kmeans_Cluster, 
                           data = new_data, mean)
write.csv(cluster_means, "cluster_means.csv", row.names = FALSE)

# 4. Extension: Hierarchical Clustering
dist_matrix <- dist(data_scaled, method = "euclidean")
hclust_result <- hclust(dist_matrix, method = "ward.D2")
png("dendrogram.png", width = 800, height = 600)
plot(hclust_result, labels = FALSE, main = "Dendrogram of Hierarchical Clustering")
dev.off()
clusters_hclust <- cutree(hclust_result, k = 3)
new_data$Hclust_Cluster <- as.factor(clusters_hclust)
write.csv(new_data, "BRCA_hclust_output.csv", row.names = FALSE)  # Lưu dữ liệu với cụm

# Chi-square cho Hierarchical Clustering
sink("chisq_hclust_results.txt")
print("Chi-square Test: Hierarchical Clusters vs Tumour_Stage")
chisq_hclust <- chisq.test(table(new_data$Hclust_Cluster, new_data$Tumour_Stage))
print(chisq_hclust)
write.csv(table(new_data$Hclust_Cluster, new_data$Tumour_Stage), "chisq_hclust_table.csv")
sink()

# So sánh K-means và Hierarchical
write.csv(table(new_data$Kmeans_Cluster, new_data$Hclust_Cluster), "kmeans_vs_hclust.csv")
