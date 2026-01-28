library(tidyverse)
library(caret)
library(glmnet)
library(factoextra)
library(Rtsne)
library(umap)
library(corrplot)
library(pheatmap)
library(cluster)
library(randomForest)
library(kohonen)
library(elasticnet)
library(ConsensusClusterPlus)
library(patchwork)
library(ggpubr)
library(CCA)
library(patchwork)

data <- read_csv("data.csv")
labels <- read_csv("labels.csv")

# ==================== 数据预处理函数 ====================
prepare_data_for_analysis <- function(data, min_var = 1e-10, verbose = TRUE) {
  if (verbose) {
    cat("=== 数据预处理 ===\n")
    cat("原始维度:", dim(data), "\n")
  }
  
  # 1. 检查并移除可能的索引列
  if (is.character(data[[1]]) && nrow(data) == length(unique(data[[1]]))) {
    if (verbose) cat("检测到样本ID列，设为行名:", colnames(data)[1], "\n")
    rownames(data) <- data[[1]]
    data <- data[, -1, drop = FALSE]
  }
  
  # 2. 只保留数值列
  numeric_cols <- sapply(data, is.numeric)
  if (!all(numeric_cols)) {
    if (verbose) cat("移除非数值列:", sum(!numeric_cols), "\n")
    data <- data[, numeric_cols, drop = FALSE]
  }
  
  # 3. 移除方差为零或接近零的列（避免PCA错误）
  if (ncol(data) > 0) {
    feature_var <- apply(data, 2, var, na.rm = TRUE)
    zero_var_features <- which(feature_var < min_var)
    
    if (length(zero_var_features) > 0) {
      if (verbose) {
        cat("移除低方差特征 (方差 <", min_var, "):", length(zero_var_features), "\n")
        if (length(zero_var_features) <= 10) {
          cat("移除的特征:", colnames(data)[zero_var_features], "\n")
        } else {
          cat("移除的特征: 前10个 -", paste(colnames(data)[head(zero_var_features, 10)], collapse = ", "), "...\n")
        }
      }
      data <- data[, -zero_var_features, drop = FALSE]
    }
  }
  
  # 4. 检查并处理缺失值
  na_count <- sum(is.na(data))
  if (na_count > 0) {
    if (verbose) cat("发现缺失值:", na_count, "\n")
    # 用列均值填充缺失值
    for (col in colnames(data)) {
      if (any(is.na(data[[col]]))) {
        col_mean <- mean(data[[col]], na.rm = TRUE)
        if (is.na(col_mean)) {  # 如果整个列都是NA
          col_mean <- 0
        }
        data[[col]][is.na(data[[col]])] <- col_mean
      }
    }
  }
  
  # 5. 检查并处理无穷大值
  inf_count <- sum(sapply(data, function(x) any(is.infinite(x))))
  if (inf_count > 0) {
    if (verbose) cat("发现无穷大值:", inf_count, "\n")
    data <- do.call(data.frame, 
                    lapply(data, function(x) replace(x, is.infinite(x), NA)))
    # 再次用均值填充
    for (col in colnames(data)) {
      if (any(is.na(data[[col]]))) {
        col_mean <- mean(data[[col]], na.rm = TRUE)
        if (is.na(col_mean)) {
          col_mean <- 0
        }
        data[[col]][is.na(data[[col]])] <- col_mean
      }
    }
  }
  
  # 6. 验证数据
  if (ncol(data) == 0) {
    stop("错误：处理后没有剩下任何特征！")
  }
  
  if (any(apply(data, 2, var, na.rm = TRUE) < min_var)) {
    warning("警告：仍然存在低方差特征")
  }
  
  if (verbose) {
    cat("处理后维度:", dim(data), "\n")
    cat("特征方差统计:\n")
    cat("  最小值:", min(apply(data, 2, var, na.rm = TRUE)), "\n")
    cat("  中位数:", median(apply(data, 2, var, na.rm = TRUE)), "\n")
    cat("  最大值:", max(apply(data, 2, var, na.rm = TRUE)), "\n")
    cat("=== 数据预处理完成 ===\n\n")
  }
  
  return(data)
}

# 预处理数据
data_clean <- prepare_data_for_analysis(data)

# 方法1：主成分分析（PCA） - 基础降维
cat("=== 开始PCA分析 ===\n")
tryCatch({
  # 检查是否还有低方差特征
  feature_var <- apply(data_clean, 2, var, na.rm = TRUE)
  if (any(feature_var < 1e-10)) {
    cat("警告：仍然存在低方差特征，手动移除\n")
    data_clean <- data_clean[, feature_var >= 1e-10, drop = FALSE]
  }
  
  # 检查数据维度
  cat("PCA输入数据维度:", dim(data_clean), "\n")
  cat("特征数:", ncol(data_clean), "\n")
  cat("样本数:", nrow(data_clean), "\n")
  
  # 运行PCA
  pca_result <- prcomp(data_clean, scale. = TRUE, center = TRUE)
  cat("✓ PCA成功完成！\n")
  cat("  生成的主成分数:", length(pca_result$sdev), "\n")
  cat("  解释的方差比例:\n")
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  cat("    PC1:", round(var_explained[1]*100, 2), "%\n")
  cat("    PC2:", round(var_explained[2]*100, 2), "%\n")
  cat("    前5个PC累计:", round(sum(var_explained[1:5])*100, 2), "%\n")
  cat("    前10个PC累计:", round(sum(var_explained[1:10])*100, 2), "%\n")
  
}, error = function(e) {
  cat("PCA错误:", e$message, "\n")
  
  # 尝试替代方法1：手动标准化
  cat("尝试手动标准化...\n")
  tryCatch({
    data_scaled <- scale(data_clean, center = TRUE, scale = TRUE)
    # 检查是否有NaN
    if (any(is.nan(data_scaled))) {
      data_scaled[is.nan(data_scaled)] <- 0
    }
    pca_result <- prcomp(data_scaled, center = FALSE, scale. = FALSE)
    cat("✓ 使用手动标准化PCA成功！\n")
  }, error = function(e2) {
    cat("手动标准化错误:", e2$message, "\n")
    
    # 尝试替代方法2：不使用scale
    cat("尝试不使用scale的PCA...\n")
    tryCatch({
      pca_result <- prcomp(data_clean, scale. = FALSE, center = TRUE)
      cat("✓ 不使用scale的PCA成功！\n")
    }, error = function(e3) {
      cat("最终PCA错误:", e3$message, "\n")
      stop("无法进行PCA分析，请检查数据质量")
    })
  })
})

# 可视化PCA结果 
if (exists("pca_result")) {
  cat("=== 生成PCA可视化 ===\n")
  tryCatch({
    # 1. 特征值图
    p1 <- fviz_eig(pca_result, 
                   ncp = min(20, length(pca_result$sdev)),
                   addlabels = TRUE,
                   barfill = "#2E9FDF",
                   barcolor = "black",
                   linecolor = "red",
                   ggtheme = theme_minimal()) +
      labs(title = "PCA - Variance Explained by Principal Components",
           x = "Principal Components", y = "Percentage of Explained Variance")
    
    ggsave("pca_eigenvalues.tiff", p1, 
           width = 8, height = 6, dpi = 300, 
           device = "tiff", compression = "lzw")
    cat("✓ PCA特征值图已保存为 pca_eigenvalues.tiff\n")
    
    # 2. 个体因子图
    p2 <- fviz_pca_ind(pca_result, 
                       col.ind = labels$Class,
                       palette = "jco",
                       addEllipses = TRUE,
                       ellipse.type = "confidence",
                       legend.title = "Class",
                       title = "PCA - Individual Factor Map",
                       ggtheme = theme_minimal())
    
    ggsave("pca_individuals.tiff", p2, 
           width = 10, height = 8, dpi = 300,
           device = "tiff", compression = "lzw")
    cat("✓ PCA个体图已保存为 pca_individuals.tiff\n")
    
    # 3. 双标图
    p3 <- fviz_pca_biplot(pca_result, 
                          col.ind = labels$Class,
                          palette = "jco",
                          addEllipses = TRUE,
                          title = "PCA Biplot",
                          ggtheme = theme_minimal())
    
    ggsave("pca_biplot.tiff", p3, 
           width = 12, height = 8, dpi = 300,
           device = "tiff", compression = "lzw")
    cat("✓ PCA双标图已保存为 pca_biplot.tiff\n")
    
  }, error = function(e) {
    cat("PCA可视化错误:", e$message, "\n")
  })
}

# 方法2：t-SNE可视化 
cat("\n=== 开始t-SNE分析 ===\n")
set.seed(123)
tryCatch({
  # 自适应perplexity
  perplexity_val <- min(30, floor((nrow(data_clean)-1)/3))
  if (perplexity_val < 5) perplexity_val <- 5
  
  cat("t-SNE参数:\n")
  cat("  维度: 2\n")
  cat("  Perplexity:", perplexity_val, "\n")
  cat("  最大迭代: 1000\n")
  
  tsne_result <- Rtsne(as.matrix(data_clean), 
                       dims = 2, 
                       perplexity = perplexity_val,
                       pca = TRUE,
                       partial_pca = TRUE,
                       max_iter = 1000,
                       verbose = TRUE)
  
  # 可视化t-SNE
  tsne_df <- data.frame(tsne_x = tsne_result$Y[,1],
                        tsne_y = tsne_result$Y[,2],
                        Class = labels$Class)
  
  p_tsne <- ggplot(tsne_df, aes(x = tsne_x, y = tsne_y, color = Class)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(level = 0.95, alpha = 0.2) +
    theme_minimal(base_size = 12) +
    labs(title = paste("t-SNE Visualization (perplexity =", perplexity_val, ")"),
         x = "t-SNE 1",
         y = "t-SNE 2") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("tsne_plot.tiff", p_tsne, 
         width = 10, height = 8, dpi = 300,
         device = "tiff", compression = "lzw")
  cat("✓ t-SNE可视化已保存为 tsne_plot.tiff\n")
  
}, error = function(e) {
  cat("t-SNE 错误:", e$message, "\n")
})

# 方法3：UMAP降维 
cat("\n=== 开始UMAP分析 ===\n")
tryCatch({
  umap_config <- umap.defaults
  umap_config$n_neighbors <- min(15, nrow(data_clean) - 1)
  if (umap_config$n_neighbors < 5) umap_config$n_neighbors <- 5
  umap_config$min_dist <- 0.1
  umap_config$n_epochs <- 500
  
  cat("UMAP参数:\n")
  cat("  n_neighbors:", umap_config$n_neighbors, "\n")
  cat("  min_dist:", umap_config$min_dist, "\n")
  cat("  n_epochs:", umap_config$n_epochs, "\n")
  
  umap_result <- umap(as.matrix(data_clean), config = umap_config)
  
  # 可视化UMAP
  umap_df <- data.frame(umap_x = umap_result$layout[,1],
                        umap_y = umap_result$layout[,2],
                        Class = labels$Class)
  
  p_umap <- ggplot(umap_df, aes(x = umap_x, y = umap_y, color = Class)) +
    geom_point(size = 2, alpha = 0.7) +
    stat_ellipse(level = 0.95, alpha = 0.2) +
    theme_minimal(base_size = 12) +
    labs(title = "UMAP Visualization",
         x = "UMAP 1",
         y = "UMAP 2") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("umap_plot.tiff", p_umap, 
         width = 10, height = 8, dpi = 300,
         device = "tiff", compression = "lzw")
  cat("✓ UMAP可视化已保存为 umap_plot.tiff\n")
  
}, error = function(e) {
  cat("UMAP错误:", e$message, "\n")
})

# 方法4：LASSO特征选择
cat("\n=== 开始LASSO特征选择 ===\n")
# 准备数据
x <- as.matrix(data_clean)
y <- as.factor(labels$Class)

# 检查类别平衡
cat("类别分布:\n")
print(table(y))

# 如果类别太少，调整参数
if (length(unique(y)) < 2) {
  stop("错误：类别数小于2，无法进行分类分析")
}

# 选择family类型
if (length(unique(y)) == 2) {
  family_type <- "binomial"
  cat("使用二项分布(binomial)\n")
} else {
  family_type <- "multinomial"
  cat("使用多项分布(multinomial)\n")
}

# 交叉验证寻找最优lambda
set.seed(123)
tryCatch({
  cv_lasso <- cv.glmnet(x, y, 
                        family = family_type,
                        alpha = 1,  # LASSO惩罚
                        nfolds = min(10, length(y)),
                        parallel = FALSE,
                        type.measure = "class")
  
  # 可视化交叉验证结果 
  tiff("lasso_cv.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
  plot(cv_lasso, main = "LASSO Cross-Validation")
  abline(v = log(cv_lasso$lambda.min), lty = 2, col = "red")
  mtext(paste("Optimal lambda:", round(cv_lasso$lambda.min, 5)), 
        side = 3, line = 0.5)
  dev.off()
  cat("✓ LASSO交叉验证图已保存为 lasso_cv.tiff\n")
  
  # 提取重要特征
  lasso_model <- glmnet(x, y, 
                        family = family_type,
                        alpha = 1,
                        lambda = cv_lasso$lambda.min)
  
  # 获取非零系数特征
  if (family_type == "multinomial") {
    coef_list <- coef(lasso_model)
    important_features <- unique(unlist(lapply(coef_list, 
                                               function(coef_mat) {
                                                 which(coef_mat[-1, ] != 0)  # 排除截距
                                               })))
  } else {
    coef_mat <- as.matrix(coef(lasso_model))
    important_features <- which(coef_mat[-1, ] != 0)  # 排除截距
  }
  
  cat("LASSO选择的特征数:", length(important_features), "\n")
  
  # 如果选择的特征太少，使用高方差特征
  if (length(important_features) < 10) {
    cat("LASSO选择的特征太少，使用高方差特征\n")
    feature_var <- apply(x, 2, var)
    important_features <- order(feature_var, decreasing = TRUE)[1:min(1000, ncol(x))]
  }
  
  # 保存特征选择结果
  important_data <- data_clean[, important_features, drop = FALSE]
  cat("重要特征数据维度:", dim(important_data), "\n")
  
}, error = function(e) {
  cat("LASSO错误:", e$message, "\n")
  cat("使用所有特征进行分析\n")
  important_features <- 1:ncol(data_clean)
  important_data <- data_clean
})

# 方法5：层次聚类
cat("\n=== 开始层次聚类分析 ===\n")
tryCatch({
  # 随机抽样一部分样本，避免计算量太大
  if (nrow(important_data) > 1000) {
    cat("数据太大，随机抽样1000个样本\n")
    set.seed(123)
    sample_idx <- sample(1:nrow(important_data), 1000)
    cluster_data <- important_data[sample_idx, ]
    cluster_labels <- labels$Class[sample_idx]
  } else {
    cluster_data <- important_data
    cluster_labels <- labels$Class
  }
  
  # 标准化
  cluster_data_scaled <- scale(cluster_data)
  
  # 计算距离矩阵
  cat("计算距离矩阵...\n")
  dist_matrix <- dist(cluster_data_scaled)
  
  # 层次聚类
  cat("进行层次聚类...\n")
  hclust_result <- hclust(dist_matrix, method = "ward.D2")
  
  # 可视化聚类树 
  tiff("hierarchical_dendrogram.tiff", width = 12, height = 8, units = "in", res = 300, compression = "lzw")
  plot(hclust_result, cex = 0.6, 
       main = "Hierarchical Clustering Dendrogram",
       xlab = "Samples", sub = "",
       labels = FALSE)
  dev.off()
  cat("✓ 聚类树已保存为 hierarchical_dendrogram.tiff\n")
  
}, error = function(e) {
  cat("层次聚类错误:", e$message, "\n")
})

# 方法6：随机森林特征重要性
cat("\n=== 开始随机森林特征重要性分析 ===\n")
tryCatch({
  set.seed(123)
  
  # 使用LASSO选择的重要特征
  cat("训练随机森林模型...\n")
  rf_model <- randomForest(x = important_data, 
                           y = y,
                           ntree = 500,
                           importance = TRUE,
                           do.trace = FALSE)
  
  # 可视化特征重要性 
  tiff("rf_importance.tiff", width = 12, height = 8, units = "in", res = 300, compression = "lzw")
  varImpPlot(rf_model, 
             main = "Random Forest - Feature Importance",
             cex = 0.8)
  dev.off()
  cat("✓ 随机森林特征重要性已保存为 rf_importance.tiff\n")
  
  # 提取重要性数据
  imp_df <- as.data.frame(rf_model$importance) %>%
    rownames_to_column("Gene") %>%
    arrange(desc(MeanDecreaseAccuracy))
  
  # 保存重要性排名
  write.csv(imp_df, "feature_importance_ranking.csv", row.names = FALSE)
  cat("✓ 特征重要性排名已保存为 feature_importance_ranking.csv\n")
  
  # 可视化前20个重要特征 
  p_rf_top20 <- ggplot(imp_df[1:min(20, nrow(imp_df)), ], 
                       aes(x = reorder(Gene, MeanDecreaseAccuracy), 
                           y = MeanDecreaseAccuracy)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    labs(title = paste("Top", min(20, nrow(imp_df)), "Important Features (Random Forest)"),
         x = "Gene",
         y = "Mean Decrease Accuracy") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave("rf_top20_features.tiff", p_rf_top20, 
         width = 10, height = 8, dpi = 300,
         device = "tiff", compression = "lzw")
  cat("✓ 前20个重要特征可视化已保存为 rf_top20_features.tiff\n")
  
}, error = function(e) {
  cat("随机森林错误:", e$message, "\n")
})

# 方法7：热图可视化
cat("\n=== 开始热图分析（安全版本） ===\n")
tryCatch({
  if (exists("important_data") && ncol(important_data) > 1) {
    
    # 1. 检查并准备数据
    cat("检查数据质量...\n")
    
    # 确保important_data是数值矩阵
    if (!is.matrix(important_data)) {
      cat("转换数据为矩阵...\n")
      important_data_matrix <- as.matrix(important_data)
    } else {
      important_data_matrix <- important_data
    }
    
    # 检查行列名
    if (is.null(rownames(important_data_matrix))) {
      rownames(important_data_matrix) <- paste0("Sample_", 1:nrow(important_data_matrix))
    }
    if (is.null(colnames(important_data_matrix))) {
      colnames(important_data_matrix) <- paste0("Gene_", 1:ncol(important_data_matrix))
    }
    
    # 2. 选择要展示的基因
    cat("选择基因...\n")
    
    # 如果有重要性数据，使用它
    if (exists("imp_df") && nrow(imp_df) > 0) {
      cat("使用重要性排名选择基因...\n")
      # 获取imp_df中的基因名
      important_genes <- imp_df$Gene
      
      # 找出实际在数据中存在的基因
      available_genes <- intersect(important_genes, colnames(important_data_matrix))
      
      if (length(available_genes) > 5) {
        # 使用前20个可用的基因
        n_genes <- min(20, length(available_genes))
        selected_genes <- available_genes[1:n_genes]
        cat("使用重要性选择", n_genes, "个基因\n")
      } else {
        cat("重要性基因太少，使用高方差基因\n")
        # 计算方差
        gene_variance <- apply(important_data_matrix, 2, var, na.rm = TRUE)
        selected_genes <- names(sort(gene_variance, decreasing = TRUE))[1:min(20, length(gene_variance))]
      }
    } else {
      cat("没有重要性数据，使用高方差基因...\n")
      # 计算方差
      gene_variance <- apply(important_data_matrix, 2, var, na.rm = TRUE)
      selected_genes <- names(sort(gene_variance, decreasing = TRUE))[1:min(20, length(gene_variance))]
    }
    
    # 3. 准备热图数据
    cat("准备热图数据...\n")
    heatmap_data <- important_data_matrix[, selected_genes, drop = FALSE]
    
    # 检查数据
    cat("热图数据维度:", dim(heatmap_data), "\n")
    cat("基因数:", length(selected_genes), "\n")
    
    # 移除任何NA或Inf值
    if (any(is.na(heatmap_data))) {
      cat("发现NA值，用列均值填充...\n")
      for (j in 1:ncol(heatmap_data)) {
        col_mean <- mean(heatmap_data[, j], na.rm = TRUE)
        if (is.na(col_mean)) col_mean <- 0
        heatmap_data[is.na(heatmap_data[, j]), j] <- col_mean
      }
    }
    
    if (any(is.infinite(heatmap_data))) {
      cat("发现无穷大值，用最大值替换...\n")
      for (j in 1:ncol(heatmap_data)) {
        col_max <- max(heatmap_data[!is.infinite(heatmap_data[, j]), j], na.rm = TRUE)
        if (is.infinite(col_max) || is.na(col_max)) col_max <- 1
        heatmap_data[is.infinite(heatmap_data[, j]), j] <- col_max
      }
    }
    
    # 4. 标准化数据（按行）
    cat("标准化数据...\n")
    tryCatch({
      heatmap_data_scaled <- t(scale(t(heatmap_data)))
      # 如果标准化产生NA/Inf，用0替代
      heatmap_data_scaled[is.na(heatmap_data_scaled)] <- 0
      heatmap_data_scaled[is.infinite(heatmap_data_scaled)] <- 0
    }, error = function(e) {
      cat("标准化失败，使用原始数据:", e$message, "\n")
      heatmap_data_scaled <- heatmap_data
    })
    
    # 5. 准备注释
    cat("准备样本注释...\n")
    annotation_df <- data.frame(
      CancerType = as.character(labels$Class),
      row.names = rownames(heatmap_data_scaled)
    )
    
    # 确保行名匹配
    if (!all(rownames(annotation_df) == rownames(heatmap_data_scaled))) {
      cat("警告：行名不匹配，重新对齐...\n")
      # 重新排序annotation_df以匹配heatmap_data_scaled
      annotation_df <- annotation_df[rownames(heatmap_data_scaled), , drop = FALSE]
    }
    
    # 6. 设置颜色
    cancer_colors <- list(
      CancerType = c(
        BRCA = "#E41A1C", 
        COAD = "#377EB8", 
        KIRC = "#4DAF4A", 
        LUAD = "#984EA3", 
        PRAD = "#FF7F00"
      )
    )
    
    # 7. 绘制热图 
    tiff("top_genes_heatmap.tiff", width = 14, height = 10, units = "in", res = 300, compression = "lzw")
    tryCatch({
      pheatmap(heatmap_data_scaled,
               main = paste("Expression Heatmap of", ncol(heatmap_data_scaled), "Important Genes"),
               color = colorRampPalette(c("blue", "white", "red"))(100),
               scale = "row",  # 按行标准化
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = FALSE,  # 不显示样本名
               show_colnames = TRUE,   # 显示基因名
               annotation_row = annotation_df,
               annotation_colors = cancer_colors,
               fontsize_col = 8,
               fontsize_row = 6,
               border_color = NA,
               treeheight_row = 50,
               treeheight_col = 50)
    }, error = function(e) {
      cat("pheatmap错误，尝试简化版本:", e$message, "\n")
      # 简化版本
      pheatmap(heatmap_data_scaled,
               main = paste("Expression Heatmap of", ncol(heatmap_data_scaled), "Genes"),
               color = colorRampPalette(c("blue", "white", "red"))(100),
               scale = "none",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               show_rownames = FALSE,
               show_colnames = TRUE,
               fontsize_col = 8)
    })
    dev.off()
    cat("✓ 热图已保存为 top_genes_heatmap.tiff\n")
    
    # 8. 额外：创建按癌症类型分组的平均表达热图
    cat("\n创建按癌症类型分组的平均表达热图...\n")
    
    # 计算每种癌症类型的平均表达
    cancer_types <- unique(labels$Class)
    mean_expression <- matrix(NA, nrow = length(cancer_types), ncol = ncol(heatmap_data))
    rownames(mean_expression) <- cancer_types
    colnames(mean_expression) <- colnames(heatmap_data)
    
    for (cancer in cancer_types) {
      cancer_samples <- which(labels$Class == cancer)
      if (length(cancer_samples) > 0) {
        mean_expression[cancer, ] <- colMeans(heatmap_data[cancer_samples, , drop = FALSE], na.rm = TRUE)
      }
    }
    
    # 标准化
    mean_expression_scaled <- t(scale(t(mean_expression)))
    mean_expression_scaled[is.na(mean_expression_scaled)] <- 0
    
    # 绘制分组热图 
    tiff("mean_expression_by_cancer.tiff", width = 12, height = 6, units = "in", res = 300, compression = "lzw")
    pheatmap(mean_expression_scaled,
             main = "Mean Gene Expression by Cancer Type",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             scale = "row",
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize_col = 7,
             fontsize_row = 10,
             border_color = "gray")
    dev.off()
    cat("✓ 分组平均表达热图已保存为 mean_expression_by_cancer.tiff\n")
    
  } else {
    cat("重要数据不可用或基因太少\n")
  }
  
}, error = function(e) {
  cat("热图分析主流程错误:", e$message, "\n")
  
  # 最后尝试：最简单的热图
  cat("尝试最基本的heatmap函数...\n")
  tryCatch({
    tiff("basic_heatmap.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")
    
    # 使用基础heatmap函数
    heatmap_data_subset <- as.matrix(important_data[, 1:min(20, ncol(important_data))])
    
    # 移除NA
    heatmap_data_subset[is.na(heatmap_data_subset)] <- 0
    
    heatmap(heatmap_data_subset,
            col = colorRampPalette(c("blue", "white", "red"))(100),
            main = "Basic Gene Expression Heatmap",
            xlab = "Genes",
            ylab = "Samples",
            scale = "row",
            labRow = FALSE,
            cexCol = 0.7)
    
    dev.off()
    cat("✓ 基础热图已保存为 basic_heatmap.tiff\n")
  }, error = function(e2) {
    cat("最终热图失败:", e2$message, "\n")
  })
})

# 方法8：箱线图可视化重要基因表达
cat("\n=== 开始重要基因表达分析 ===\n")
tryCatch({
  if (exists("imp_df") && nrow(imp_df) >= 5) {
    # 使用前5个基因
    n_genes_box <- min(5, nrow(imp_df))
    top_genes <- imp_df$Gene[1:n_genes_box]
    
    cat("分析基因:", paste(top_genes, collapse = ", "), "\n")
    
    # 准备绘图数据
    plot_data <- important_data[, top_genes, drop = FALSE] %>%
      as.data.frame() %>%
      rownames_to_column("Sample") %>%
      pivot_longer(cols = -Sample, 
                   names_to = "Gene", 
                   values_to = "Expression",
                   names_ptypes = list(Gene = factor(levels = top_genes))) %>%
      left_join(data.frame(Sample = rownames(important_data), 
                           Class = labels$Class), 
                by = "Sample")
    
    # 设置因子顺序
    plot_data$Gene <- factor(plot_data$Gene, levels = top_genes)
    
    # 绘制箱线图 
    p_boxplot <- ggplot(plot_data, aes(x = Class, y = Expression, fill = Class)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
      scale_fill_brewer(palette = "Set2") +
      theme_minimal(base_size = 12) +
      labs(title = paste("Expression of Top", n_genes_box, "Important Genes by Cancer Type"),
           x = "Cancer Type",
           y = "Expression Level") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right",
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave("top_genes_boxplot.tiff", p_boxplot, 
           width = 14, height = 10, dpi = 300,
           device = "tiff", compression = "lzw")
    cat("✓ 基因表达箱线图已保存为 top_genes_boxplot.tiff\n")
  } else {
    cat("没有足够的重要性数据进行箱线图分析\n")
  }
}, error = function(e) {
  cat("箱线图绘制错误:", e$message, "\n")
})

# 方法9：多方法降维比较
cat("\n=== 开始多方法降维比较 ===\n")
tryCatch({
  # 准备降维结果
  dim_reduction_results <- list(
    PCA = pca_result$x[, 1:2],
    tSNE = tsne_result$Y,
    UMAP = umap_result$layout
  )
  
  # 分别保存每个降维方法的图
  for (method_name in names(dim_reduction_results)) {
    dim_data <- dim_reduction_results[[method_name]]
    
    plot_df <- data.frame(
      Dim1 = dim_data[, 1],
      Dim2 = dim_data[, 2],
      Class = labels$Class
    )
    
    p <- ggplot(plot_df, aes(x = Dim1, y = Dim2, color = Class)) +
      geom_point(size = 1.5, alpha = 0.6) +
      stat_ellipse(level = 0.95, alpha = 0.2) +
      theme_minimal() +
      labs(title = paste(method_name, "Visualization"),
           x = paste(method_name, "1"),
           y = paste(method_name, "2")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "right")
    
    filename <- paste0(tolower(method_name), "_visualization.tiff")
    ggsave(filename, p, 
           width = 8, height = 6, dpi = 300,
           device = "tiff", compression = "lzw")
    cat("✓ ", method_name, "可视化已保存为 ", filename, "\n", sep = "")
  }
  
  # 创建图例文件
  p_legend <- ggplot(plot_df, aes(x = Dim1, y = Dim2, color = Class)) +
    geom_point(size = 2) +
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
  
  ggsave("dimension_reduction_legend.tiff", p_legend, 
         width = 4, height = 6, dpi = 300,
         device = "tiff", compression = "lzw")
  cat("✓ 图例已保存为 dimension_reduction_legend.tiff\n")
  
}, error = function(e) {
  cat("降维比较错误:", e$message, "\n")
})

# 方法10：简化版机器学习模型评估
cat("\n=== 开始简化版机器学习模型评估 ===\n")
tryCatch({
  set.seed(123)
  
  # 准备数据
  ml_data <- cbind(important_data, Class = y)
  
  # 划分训练测试集
  train_index <- createDataPartition(y, p = 0.7, list = FALSE)
  train_data <- ml_data[train_index, ]
  test_data <- ml_data[-train_index, ]
  
  # 只训练随机森林模型
  cat("训练随机森林模型...\n")
  rf_fit <- train(Class ~ ., 
                  data = train_data,
                  method = "rf",
                  trControl = trainControl(method = "cv", number = 5),
                  tuneLength = 3,
                  ntree = 200)
  
  # 预测
  test_pred_rf <- predict(rf_fit, test_data)
  
  # 混淆矩阵
  cm_rf <- confusionMatrix(test_pred_rf, test_data$Class)
  
  # 保存模型性能
  model_performance <- data.frame(
    Model = "RandomForest",
    Accuracy = cm_rf$overall["Accuracy"],
    Kappa = cm_rf$overall["Kappa"],
    Sensitivity = mean(cm_rf$byClass[,"Sensitivity"], na.rm = TRUE),
    Specificity = mean(cm_rf$byClass[,"Specificity"], na.rm = TRUE)
  )
  
  write.csv(model_performance, "model_performance_simple.csv", row.names = FALSE)
  cat("✓ 模型性能已保存为 model_performance_simple.csv\n")
  
  # 绘制混淆矩阵 
  cm_df <- as.data.frame(cm_rf$table)
  p_cm <- ggplot(data = cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), vjust = 1, size = 4) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme_minimal(base_size = 12) +
    labs(title = "Random Forest Confusion Matrix",
         x = "True Class",
         y = "Predicted Class",
         fill = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("confusion_matrix.tiff", p_cm, 
         width = 10, height = 8, dpi = 300,
         device = "tiff", compression = "lzw")
  cat("✓ 混淆矩阵已保存为 confusion_matrix.tiff\n")
  
  cat("\n模型性能摘要:\n")
  print(model_performance)
  cat("\n混淆矩阵:\n")
  print(cm_rf$table)
  
}, error = function(e) {
  cat("机器学习模型评估错误:", e$message, "\n")
  # 尝试更简单的方法
  tryCatch({
    cat("尝试简单随机森林...\n")
    set.seed(123)
    
    # 直接使用randomForest
    rf_simple <- randomForest(x = important_data[train_index, ],
                              y = y[train_index],
                              xtest = important_data[-train_index, ],
                              ytest = y[-train_index],
                              ntree = 200)
    
    # 计算准确率
    accuracy <- sum(rf_simple$test$predicted == y[-train_index]) / length(y[-train_index])
    cat("随机森林测试准确率:", round(accuracy, 3), "\n")
    
    simple_perf <- data.frame(
      Model = "RandomForest_Simple",
      Accuracy = accuracy
    )
    write.csv(simple_perf, "simple_model_performance.csv", row.names = FALSE)
    
  }, error = function(e2) {
    cat("简单模型也失败:", e2$message, "\n")
  })
})

# 方法11：特征相关性分析
cat("\n=== 开始特征相关性分析 ===\n")
tryCatch({
  if (ncol(important_data) > 1) {
    # 计算特征间的相关性
    cor_matrix <- cor(important_data)
    
    # 可视化相关性矩阵 
    tiff("feature_correlation.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")
    corrplot(cor_matrix, 
             method = "color",
             type = "upper",
             tl.col = "black",
             tl.srt = 45,
             addCoef.col = "black",
             number.cex = 0.5,
             title = "Correlation between Important Features")
    dev.off()
    cat("✓ 特征相关性图已保存为 feature_correlation.tiff\n")
    
    # 提取高度相关的特征对
    cor_threshold <- 0.7
    high_cor_pairs <- which(abs(cor_matrix) > cor_threshold & 
                              upper.tri(cor_matrix), arr.ind = TRUE)
    
    if (nrow(high_cor_pairs) > 0) {
      high_cor_df <- data.frame(
        Gene1 = colnames(cor_matrix)[high_cor_pairs[, 1]],
        Gene2 = colnames(cor_matrix)[high_cor_pairs[, 2]],
        Correlation = cor_matrix[high_cor_pairs]
      ) %>%
        arrange(desc(abs(Correlation)))
      
      write.csv(high_cor_df, "high_correlation_pairs.csv", row.names = FALSE)
      cat("✓ 高相关性基因对已保存为 high_correlation_pairs.csv\n")
      cat("  发现", nrow(high_cor_df), "对高度相关基因 (|r| >", cor_threshold, ")\n")
    }
  }
}, error = function(e) {
  cat("相关性分析错误:", e$message, "\n")
})

# 方法12：PCA loadings分析
cat("\n=== 开始PCA loadings分析 ===\n")
tryCatch({
  # 提取前5个主成分的loadings
  n_pc <- min(5, ncol(pca_result$rotation))
  loadings_df <- as.data.frame(pca_result$rotation[, 1:n_pc])
  colnames(loadings_df) <- paste0("PC", 1:n_pc)
  loadings_df$Gene <- rownames(loadings_df)
  
  # 对每个主成分，提取最重要的基因
  for (i in 1:n_pc) {
    pc_name <- paste0("PC", i)
    pc_loadings <- loadings_df[, c("Gene", pc_name)]
    pc_loadings <- pc_loadings[order(abs(pc_loadings[, pc_name]), decreasing = TRUE), ]
    
    # 保存前20个基因
    top_genes <- head(pc_loadings, 20)
    write.csv(top_genes, 
              paste0("pc", i, "_top_loadings.csv"), 
              row.names = FALSE)
    
    cat("✓ PC", i, "的loadings已保存为 pc", i, "_top_loadings.csv\n", sep = "")
  }
  
  # 可视化PC1 vs PC2 loadings 
  p_loadings <- ggplot(loadings_df, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    theme_minimal() +
    labs(title = "PCA Loadings Plot (PC1 vs PC2)",
         x = "PC1 Loadings",
         y = "PC2 Loadings")
  
  ggsave("pca_loadings_plot.tiff", p_loadings, 
         width = 10, height = 8, dpi = 300,
         device = "tiff", compression = "lzw")
  cat("✓ PCA loadings图已保存为 pca_loadings_plot.tiff\n")
  
  # 添加特征载荷条形图
  p_loadings_bar <- ggplot(head(loadings_df[order(abs(loadings_df$PC1), decreasing = TRUE), ], 20), 
                           aes(x = reorder(Gene, PC1), y = PC1)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    theme_minimal(base_size = 10) +
    labs(title = "Top 20 Gene Loadings on PC1",
         x = "Gene",
         y = "PC1 Loading")
  
  ggsave("pc1_top_loadings_bar.tiff", p_loadings_bar, 
         width = 10, height = 8, dpi = 300,
         device = "tiff", compression = "lzw")
  cat("✓ PC1特征载荷条形图已保存为 pc1_top_loadings_bar.tiff\n")
  
}, error = function(e) {
  cat("PCA loadings分析错误:", e$message, "\n")
})

# 生成最终汇总报告
cat("\n=== 生成最终汇总报告 ===\n")
final_summary <- data.frame(
  Metric = c(
    "Total Samples",
    "Total Genes",
    "Genes after Filtering",
    "Cancer Types",
    "BRCA Samples",
    "COAD Samples", 
    "KIRC Samples",
    "LUAD Samples",
    "PRAD Samples",
    "LASSO Selected Genes",
    "PCA Components (80% Variance)",
    "PCA Components (90% Variance)",
    "Top Gene (Random Forest)",
    "RF MeanDecreaseAccuracy"
  ),
  Value = c(
    nrow(data_clean),
    ncol(data),
    ncol(data_clean),
    length(unique(y)),
    sum(y == "BRCA"),
    sum(y == "COAD"),
    sum(y == "KIRC"),
    sum(y == "LUAD"),
    sum(y == "PRAD"),
    length(important_features),
    which(summary(pca_result)$importance[3,] >= 0.8)[1],
    which(summary(pca_result)$importance[3,] >= 0.9)[1],
    if (exists("imp_df") && nrow(imp_df) > 0) imp_df$Gene[1] else "N/A",
    if (exists("imp_df") && nrow(imp_df) > 0) round(imp_df$MeanDecreaseAccuracy[1], 3) else "N/A"
  )
)

write.csv(final_summary, "final_analysis_summary.csv", row.names = FALSE)
cat("📁 已生成的分析文件:\n")
cat("\n1. 🖼️ 可视化文件 (TIFF格式):\n")
cat("   ├── pca_eigenvalues.tiff               - PCA特征值图\n")
cat("   ├── pca_individuals.tiff               - PCA个体图\n")
cat("   ├── pca_biplot.tiff                    - PCA双标图\n")
cat("   ├── tsne_plot.tiff                     - t-SNE可视化\n")
cat("   ├── umap_plot.tiff                     - UMAP可视化\n")
cat("   ├── pca_visualization.tiff             - PCA可视化\n")
cat("   ├── tsne_visualization.tiff            - t-SNE可视化\n")
cat("   ├── umap_visualization.tiff            - UMAP可视化\n")
cat("   ├── dimension_reduction_legend.tiff    - 降维图例\n")
cat("   ├── lasso_cv.tiff                      - LASSO交叉验证\n")
cat("   ├── hierarchical_dendrogram.tiff       - 层次聚类树\n")
cat("   ├── rf_importance.tiff                 - 随机森林特征重要性\n")
cat("   ├── rf_top20_features.tiff            - 前20个重要特征\n")
cat("   ├── top_genes_heatmap.tiff             - 重要基因热图\n")
cat("   ├── mean_expression_by_cancer.tiff     - 分组平均表达热图\n")
cat("   ├── top_genes_boxplot.tiff             - 基因表达箱线图\n")
cat("   ├── confusion_matrix.tiff              - 混淆矩阵\n")
cat("   ├── feature_correlation.tiff           - 特征相关性图\n")
cat("   ├── pca_loadings_plot.tiff             - PCA loadings图\n")
cat("   └── pc1_top_loadings_bar.tiff          - PC1特征载荷条形图\n")
cat("\n2. 📁 数据文件 (CSV格式):\n")
cat("   ├── feature_importance_ranking.csv     - 特征重要性排名\n")
cat("   ├── model_performance_simple.csv       - 模型性能\n")
cat("   ├── high_correlation_pairs.csv         - 高相关性基因对\n")
cat("   ├── final_analysis_summary.csv         - 最终汇总\n")
cat("   ├── pc1_top_loadings.csv               - PC1 loadings\n")
cat("   ├── pc2_top_loadings.csv               - PC2 loadings\n")
cat("   ├── pc3_top_loadings.csv               - PC3 loadings\n")
cat("   ├── pc4_top_loadings.csv               - PC4 loadings\n")
cat("   └── pc5_top_loadings.csv               - PC5 loadings\n")
cat("\n3. 💾 工作空间文件:\n")
cat("   └── complete_analysis.RData            - 完整R工作空间\n")
cat("\n4. 📈 分析摘要:\n")
cat("   - 样本总数: ", nrow(data_clean), "\n")
cat("   - 基因总数: ", ncol(data_clean), "\n")
cat("   - 癌症类型: ", length(unique(y)), "种\n")
cat("   - LASSO选择基因: ", length(important_features), "个\n")
cat("   - PCA解释80%方差所需主成分: ", which(summary(pca_result)$importance[3,] >= 0.8)[1], "\n")
if (exists("imp_df") && nrow(imp_df) > 0) {
  cat("   - 最重要的基因: ", imp_df$Gene[1], "\n")
  cat("   - 其重要性得分: ", round(imp_df$MeanDecreaseAccuracy[1], 3), "\n")
}