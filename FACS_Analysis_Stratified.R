# ==============================================================================
# FACS GFP分析 - 分层策略（针对批次-组混淆的情况）
# ==============================================================================
#
# 研究目标：比较 Bone_marrow vs Culture vs Sabq 的GFP荧光强度差异
#
# 问题：批次与组完全混淆
#   - Bone_marrow & Culture → BatchA (20250131)
#   - Sabq → BatchC (20250904)
#
# 解决方案：分层分析
#   1. BatchA内部比较：Bone_marrow vs Culture（无批次效应，最可靠）
#   2. 跨批次比较：使用归一化值，但需要谨慎解释
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
})

# 加载数据
cat("加载数据...\n")
if (file.exists("batch_effect_analysis/ALL_CELL_DATA_with_normalization.rds")) {
  data <- readRDS("batch_effect_analysis/ALL_CELL_DATA_with_normalization.rds")
} else if (file.exists("out_exports") && length(list.files("out_exports", pattern = "RData$")) > 0) {
  rdata_file <- list.files("out_exports", pattern = "RData$", full.names = TRUE)[1]
  load(rdata_file)
  data <- ALL_CELL_DATA
} else {
  stop("未找到数据文件。请先运行主分析脚本。")
}

cat(sprintf("数据加载完成：%d个细胞\n\n", nrow(data)))

# 创建输出目录
dir.create("stratified_analysis", showWarnings = FALSE)

# ==============================================================================
# 分析1：BatchA内部比较（Bone_marrow vs Culture）- 最可靠
# ==============================================================================

cat("=== 分析1：BatchA内部比较（无批次效应）===\n\n")

data_batchA <- data %>%
  filter(batch == "BatchA")

cat(sprintf("BatchA细胞数：%d\n", nrow(data_batchA)))
cat("组分布：\n")
print(table(data_batchA$group))

# 统计检验：Bone_marrow vs Culture
cat("\n--- 统计检验：Bone_marrow vs Culture ---\n")

bm_gfp <- data_batchA %>% filter(group == "Bone_marrow") %>% pull(GFP_intensity)
culture_gfp <- data_batchA %>% filter(group == "Culture") %>% pull(GFP_intensity)

# Wilcoxon秩和检验（非参数检验，适合荧光强度数据）
wilcox_result <- wilcox.test(bm_gfp, culture_gfp)

cat(sprintf("Bone_marrow中位数：%.2f\n", median(bm_gfp)))
cat(sprintf("Culture中位数：%.2f\n", median(culture_gfp)))
cat(sprintf("Wilcoxon检验 p值：%.2e\n", wilcox_result$p.value))

if (wilcox_result$p.value < 0.001) {
  cat("*** 差异极显著（p < 0.001）\n")
} else if (wilcox_result$p.value < 0.01) {
  cat("** 差异显著（p < 0.01）\n")
} else if (wilcox_result$p.value < 0.05) {
  cat("* 差异显著（p < 0.05）\n")
} else {
  cat("无显著差异（p >= 0.05）\n")
}

# 可视化
p1 <- ggplot(data_batchA, aes(x = group, y = GFP_intensity, fill = group)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.size = 0.5) +
  scale_y_log10() +
  labs(title = "BatchA内部比较：Bone_marrow vs Culture",
       subtitle = sprintf("Wilcoxon检验 p = %.2e（无批次效应，结果可靠）", wilcox_result$p.value),
       y = "GFP荧光强度 (log10)",
       x = "") +
  theme_bw() +
  theme(legend.position = "none")

print(p1)
ggsave("stratified_analysis/1_BatchA_comparison.png", p1, width = 8, height = 6, dpi = 150)

# ==============================================================================
# 分析2：全局比较（使用归一化值）- 需要谨慎解释
# ==============================================================================

cat("\n\n=== 分析2：全局比较（包括Sabq，使用归一化值）===\n\n")

# 检查是否有归一化值
has_normalized <- "GFP_intensity_normalized" %in% colnames(data) &&
                  sum(!is.na(data$GFP_intensity_normalized)) > 0

if (has_normalized) {
  cat("✓ 使用批次矫正后的归一化值\n")
  cat("⚠️ 注意：由于批次与组混淆，归一化可能影响组间真实差异\n\n")

  # 使用归一化值
  data_plot <- data %>%
    filter(!is.na(GFP_intensity_normalized))

  # 统计汇总
  summary_stats <- data_plot %>%
    group_by(group) %>%
    summarise(
      n = n(),
      median = median(GFP_intensity_normalized),
      q25 = quantile(GFP_intensity_normalized, 0.25),
      q75 = quantile(GFP_intensity_normalized, 0.75),
      .groups = "drop"
    )

  cat("各组归一化GFP强度统计：\n")
  print(summary_stats)

  # 成对比较（使用归一化值）
  cat("\n--- 成对比较（Wilcoxon检验）---\n")

  comparisons <- list(
    c("Bone_marrow", "Culture"),
    c("Bone_marrow", "Sabq"),
    c("Culture", "Sabq")
  )

  for (comp in comparisons) {
    g1 <- comp[1]
    g2 <- comp[2]

    d1 <- data_plot %>% filter(group == g1) %>% pull(GFP_intensity_normalized)
    d2 <- data_plot %>% filter(group == g2) %>% pull(GFP_intensity_normalized)

    test_result <- wilcox.test(d1, d2)

    cat(sprintf("\n%s vs %s:\n", g1, g2))
    cat(sprintf("  中位数：%.2f vs %.2f\n", median(d1), median(d2)))
    cat(sprintf("  p值：%.2e", test_result$p.value))

    if (test_result$p.value < 0.05) {
      cat(" *")
    }

    # 检查是否跨批次
    batch_g1 <- unique(data_plot %>% filter(group == g1) %>% pull(batch))
    batch_g2 <- unique(data_plot %>% filter(group == g2) %>% pull(batch))

    if (length(intersect(batch_g1, batch_g2)) == 0) {
      cat(" ⚠️ 跨批次比较，可能受批次效应影响")
    }
    cat("\n")
  }

  # 可视化（归一化值）
  p2 <- ggplot(data_plot, aes(x = group, y = GFP_intensity_normalized, fill = group)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.size = 0.5) +
    scale_y_log10() +
    labs(title = "全局比较：所有组的GFP荧光强度",
       subtitle = "使用批次矫正后的归一化值（Sabq的比较需谨慎解释）",
       y = "GFP荧光强度（归一化）(log10)",
       x = "") +
    theme_bw() +
    theme(legend.position = "none")

  print(p2)
  ggsave("stratified_analysis/2_Global_comparison_normalized.png", p2,
         width = 8, height = 6, dpi = 150)

} else {
  cat("✗ 未找到归一化值，使用原始值（不推荐跨批次比较）\n\n")

  # 使用原始值（不推荐）
  summary_stats <- data %>%
    group_by(group) %>%
    summarise(
      n = n(),
      median = median(GFP_intensity),
      q25 = quantile(GFP_intensity, 0.25),
      q75 = quantile(GFP_intensity, 0.75),
      .groups = "drop"
    )

  cat("各组原始GFP强度统计：\n")
  print(summary_stats)

  cat("\n⚠️ 警告：使用原始值进行跨批次比较不可靠！\n")
  cat("建议：先成功运行批次矫正，或仅比较BatchA内部的组。\n")
}

# ==============================================================================
# 分析3：批次效应可视化
# ==============================================================================

cat("\n\n=== 分析3：批次效应可视化 ===\n\n")

# 批次间的差异
p3 <- ggplot(data, aes(x = batch, y = GFP_intensity, fill = batch)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.size = 0.5) +
  scale_y_log10() +
  labs(title = "批次效应：不同批次的GFP强度分布",
       subtitle = "BatchA和BatchC的细胞来自不同的实验组，无法区分批次效应和生物学差异",
       y = "GFP荧光强度 (log10)",
       x = "批次") +
  theme_bw() +
  theme(legend.position = "none")

print(p3)
ggsave("stratified_analysis/3_Batch_effect.png", p3, width = 8, height = 6, dpi = 150)

# 组-批次交互图
p4 <- ggplot(data, aes(x = group, y = GFP_intensity, fill = batch)) +
  geom_violin(alpha = 0.7, position = position_dodge(width = 0.9)) +
  scale_y_log10() +
  labs(title = "组与批次的关系",
       subtitle = "每个组只来自一个批次（完全混淆）",
       y = "GFP荧光强度 (log10)",
       x = "",
       fill = "批次") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p4)
ggsave("stratified_analysis/4_Group_batch_interaction.png", p4,
       width = 10, height = 6, dpi = 150)

# ==============================================================================
# 总结报告
# ==============================================================================

cat("\n\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("=== 分析总结与建议 ===\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

cat("数据结构：\n")
cat("  - BatchA (20250131): Bone_marrow + Culture\n")
cat("  - BatchC (20250904): Sabq\n")
cat("  - 批次与组完全混淆\n\n")

cat("分析结果：\n\n")

cat("1. BatchA内部比较（Bone_marrow vs Culture）：\n")
cat(sprintf("   - 无批次效应，结果可靠 ✓\n"))
cat(sprintf("   - Wilcoxon检验 p = %.2e\n", wilcox_result$p.value))
if (wilcox_result$p.value < 0.05) {
  cat("   - 两组GFP强度有显著差异\n\n")
} else {
  cat("   - 两组GFP强度无显著差异\n\n")
}

cat("2. 跨批次比较（Sabq vs 其他组）：\n")
cat("   - ⚠️ 结果需要谨慎解释\n")
cat("   - 观察到的差异可能是：\n")
cat("     a) 真实的生物学差异\n")
cat("     b) 技术性批次效应（不同时间上机）\n")
cat("     c) 两者的混合\n\n")

cat("建议：\n\n")

cat("方案A（推荐）：保守解释\n")
cat("  - 仅报告BatchA内部的比较（Bone_marrow vs Culture）\n")
cat("  - Sabq的比较作为探索性结果，明确说明可能受批次效应影响\n")
cat("  - 建议在结果中注明：'需要独立实验验证'\n\n")

cat("方案B：技术验证\n")
cat("  - 在BatchA中也收集Sabq样本（如果可能）\n")
cat("  - 或在BatchC中也收集Culture/Bone_marrow样本\n")
cat("  - 这样可以区分批次效应和生物学差异\n\n")

cat("方案C：统计模型\n")
cat("  - 使用线性混合模型，将batch作为随机效应\n")
cat("  - 需要足够的生物学重复（每组至少3个独立样本）\n")
cat("  - 仍然存在较大不确定性\n\n")

cat("生成的文件：\n")
cat("  - stratified_analysis/1_BatchA_comparison.png\n")
cat("  - stratified_analysis/2_Global_comparison_normalized.png\n")
cat("  - stratified_analysis/3_Batch_effect.png\n")
cat("  - stratified_analysis/4_Group_batch_interaction.png\n\n")

cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("分析完成！\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

# 保存统计结果
if (has_normalized) {
  write.csv(summary_stats, "stratified_analysis/summary_statistics.csv", row.names = FALSE)
  cat("统计汇总已保存：stratified_analysis/summary_statistics.csv\n")
}
