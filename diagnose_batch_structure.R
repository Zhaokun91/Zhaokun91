# ==============================================================================
# 诊断批次与组的混淆程度
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# 检查是否有保存的数据
if (file.exists("batch_effect_analysis/ALL_CELL_DATA_with_normalization.csv")) {
  cat("加载已保存的数据...\n")
  data <- read.csv("batch_effect_analysis/ALL_CELL_DATA_with_normalization.csv")
} else if (exists("ALL_CELL_DATA")) {
  cat("使用当前环境中的数据...\n")
  data <- ALL_CELL_DATA
} else {
  stop("未找到数据。请先运行主分析脚本。")
}

cat("\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("=== 批次与组结构诊断 ===\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

# 1. 批次与组的交叉表
cat("1. 批次与组的细胞数分布：\n\n")
cross_table <- table(data$batch, data$group)
print(cross_table)

cat("\n比例（%）：\n")
cross_prop <- prop.table(cross_table, margin = 1) * 100
print(round(cross_prop, 2))

# 2. 混淆程度检查
cat("\n\n2. 混淆程度检查：\n\n")
n_groups <- length(unique(data$group))
n_batches <- length(unique(data$batch))
n_combinations <- sum(cross_table > 0)

cat(sprintf("组数：%d\n", n_groups))
cat(sprintf("批次数：%d\n", n_batches))
cat(sprintf("实际的 batch-group 组合数：%d\n", n_combinations))
cat(sprintf("理想的组合数（无混淆）：%d\n", n_groups * n_batches))

if (n_combinations == n_groups && n_groups == n_batches) {
  cat("\n❌ 完全混淆！每个组只来自一个批次，批次效应无法与组间差异区分。\n")
  cat("   建议：不进行批次矫正，或者添加跨批次的对照组。\n")
} else if (n_combinations < n_groups * n_batches * 0.5) {
  cat("\n⚠️ 严重混淆！超过一半的 batch-group 组合缺失。\n")
  cat("   批次矫正效果可能不佳。\n")
} else {
  cat("\n✓ 混淆程度可接受，可以进行批次矫正。\n")
}

# 3. 批次平衡性检查
cat("\n\n3. 批次平衡性检查：\n\n")
batch_counts <- table(data$batch)
print(batch_counts)

cat("\n比例（%）：\n")
batch_prop <- prop.table(batch_counts) * 100
print(round(batch_prop, 2))

min_prop <- min(batch_prop)
max_prop <- max(batch_prop)
balance_ratio <- min_prop / max_prop

cat(sprintf("\n最小批次占比：%.2f%%\n", min_prop))
cat(sprintf("平衡度：%.3f（最小/最大）\n", balance_ratio))

if (balance_ratio < 0.1) {
  cat("\n❌ 批次严重不平衡（< 10%）！小批次的矫正可能不准确。\n")
  cat("   建议：考虑下采样大批次或移除过小的批次。\n")
} else if (balance_ratio < 0.3) {
  cat("\n⚠️ 批次中度不平衡（10-30%）。矫正时需谨慎。\n")
} else {
  cat("\n✓ 批次平衡性良好（> 30%）。\n")
}

# 4. 可视化
cat("\n\n4. 生成诊断图...\n")

# 批次-组分布图
p1 <- ggplot(data, aes(x = group, fill = batch)) +
  geom_bar(position = "dodge") +
  scale_y_log10() +
  labs(title = "各组的批次分布",
       subtitle = "理想情况：每个组在所有批次中都有细胞",
       y = "细胞数（log10）") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 批次比例图
p2 <- ggplot(data, aes(x = batch, fill = group)) +
  geom_bar(position = "fill") +
  labs(title = "各批次的组构成",
       subtitle = "理想情况：各批次有相似的组构成",
       y = "比例") +
  theme_bw()

library(cowplot)
combined <- plot_grid(p1, p2, ncol = 2)
print(combined)

ggsave("batch_group_diagnosis.png", combined, width = 12, height = 5, dpi = 150)
cat("诊断图已保存：batch_group_diagnosis.png\n")

# 5. 建议
cat("\n\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("=== 建议 ===\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

if (n_combinations == n_groups && n_groups == n_batches) {
  cat("由于批次与组完全混淆，建议采取以下措施之一：\n\n")
  cat("方案1（推荐）：添加跨批次的对照数据\n")
  cat("  - 在每个批次中都包含相同的对照组（如野生型细胞）\n")
  cat("  - 这样可以用对照组来估计批次效应\n\n")

  cat("方案2：不进行批次矫正\n")
  cat("  - 在R脚本中设置 DO_BATCH_CORRECTION <- FALSE\n")
  cat("  - 在下游统计分析中将batch作为协变量（如线性模型中）\n\n")

  cat("方案3：仅比较来自同一批次的组\n")
  cat("  - 例如，仅比较BatchA中的Bone_marrow vs Culture\n")
  cat("  - 不与BatchC的Sabq进行比较\n\n")

  cat("方案4：使用covar=NULL的批次矫正（不推荐）\n")
  cat("  - 这会混淆批次效应和真实的生物学差异\n")
  cat("  - 可能导致假阴性（真实差异被矫正掉）\n\n")

} else if (balance_ratio < 0.1) {
  cat("由于批次严重不平衡，建议：\n\n")
  cat("  - 下采样BatchA到与BatchC相近的大小\n")
  cat("  - 或移除BatchC，仅分析BatchA内部的比较\n")
  cat("  - 或增加BatchC的样本量\n\n")
}

cat("\n如需进一步帮助，请提供：\n")
cat("  1. BatchB（如果有）的样本信息\n")
cat("  2. 你的研究问题和比较目标\n")
cat("  3. 实验设计详情\n\n")

cat(paste0(rep("=", 80), collapse = ""), "\n")
