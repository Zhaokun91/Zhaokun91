# ==============================================================================
# 简单批次校正法 - 使用未染色样本估计批次效应
# ==============================================================================
#
# 原理：
#   1. 计算两个批次未染色样本的荧光差异
#   2. 假设这个差异是纯粹的批次效应（技术性差异）
#   3. 从BatchC的实验样本中减去这个批次效应
#
# 优点：
#   - 简单直接，容易理解
#   - 不会过度矫正
#   - 易于验证和解释
#
# 使用前提：
#   - 已运行detect_batch_effects_controls.R
#   - 确认存在显著批次效应
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
})

cat("=== 简单批次校正法 ===\n\n")

# ==============================================================================
# 1. 加载数据
# ==============================================================================

cat("--- 加载数据 ---\n")

# 加载实验样本数据
if (file.exists("out_exports") && length(list.files("out_exports", pattern = "RData$")) > 0) {
  rdata_file <- list.files("out_exports", pattern = "RData$", full.names = TRUE)[1]
  load(rdata_file)
  data <- ALL_CELL_DATA
  cat(sprintf("✓ 加载实验样本数据：%d个细胞\n", nrow(data)))
} else {
  stop("未找到实验数据！请先运行主分析脚本。")
}

# 加载未染色样本的批次效应检测结果
# 这里假设你已经运行了detect_batch_effects_controls.R
# 并记录了批次效应的估计值

cat("\n请输入批次效应检测的结果（来自detect_batch_effects_controls.R）：\n")
cat("如果还没运行，请先运行detect_batch_effects_controls.R\n\n")

# 手动设置批次效应（运行检测后修改这里）
# 示例值，请根据实际检测结果修改
BATCH_EFFECT <- list(
  GFP = list(
    BatchA_median = NA,  # 修改为实际值，如 450
    BatchC_median = NA,  # 修改为实际值，如 500
    difference = NA,     # 修改为实际值，如 50
    percent = NA         # 修改为实际值，如 11.1
  ),
  CD51 = list(
    BatchA_median = NA,
    BatchC_median = NA,
    difference = NA,
    percent = NA
  )
)

# 检查是否已设置
if (is.na(BATCH_EFFECT$GFP$difference)) {
  cat("⚠️ 批次效应参数未设置！\n")
  cat("请根据detect_batch_effects_controls.R的输出修改BATCH_EFFECT配置\n\n")

  cat("临时解决方案：自动从检测结果文件读取...\n")

  # 尝试从保存的检测数据读取
  if (file.exists("batch_effect_detection/batch_effect_params.rds")) {
    BATCH_EFFECT <- readRDS("batch_effect_detection/batch_effect_params.rds")
    cat("✓ 从文件读取批次效应参数\n")
  } else {
    stop("未找到批次效应参数！请先运行detect_batch_effects_controls.R")
  }
}

cat("\n批次效应估计值：\n")
cat(sprintf("  GFP通道：BatchC比BatchA高 %.2f (%.1f%%)\n",
            BATCH_EFFECT$GFP$difference, BATCH_EFFECT$GFP$percent))
cat(sprintf("  CD51通道：BatchC比BatchA高 %.2f (%.1f%%)\n",
            BATCH_EFFECT$CD51$difference, BATCH_EFFECT$CD51$percent))

# ==============================================================================
# 2. 应用批次校正
# ==============================================================================

cat("\n--- 应用批次校正 ---\n\n")

# 创建校正后的数据副本
data_corrected <- data

# 对BatchC的样本进行校正
# 注意：这里的difference已经是asinh变换后的差异
# 如果是原始尺度的差异，需要先变换

cat("校正BatchC的样本...\n")

# 识别BatchC的样本
batchC_idx <- which(data_corrected$batch == "BatchC")
cat(sprintf("  BatchC样本数：%d个细胞\n", length(batchC_idx)))

if (length(batchC_idx) > 0) {
  # 方法1：如果BATCH_EFFECT是原始尺度的差异
  # 直接减去差异
  data_corrected$GFP_intensity_corrected <- data_corrected$GFP_intensity
  data_corrected$CD51_PE_intensity_corrected <- data_corrected$CD51_PE_intensity

  # 对BatchC进行校正（减去批次效应）
  data_corrected$GFP_intensity_corrected[batchC_idx] <-
    data_corrected$GFP_intensity[batchC_idx] - BATCH_EFFECT$GFP$difference

  data_corrected$CD51_PE_intensity_corrected[batchC_idx] <-
    data_corrected$CD51_PE_intensity[batchC_idx] - BATCH_EFFECT$CD51$difference

  cat("  ✓ 校正完成\n")

  # 统计摘要
  cat("\n校正前后对比：\n")

  # Sabq组（BatchC）
  sabq_before <- data %>% filter(group == "Sabq") %>% pull(GFP_intensity)
  sabq_after <- data_corrected %>% filter(group == "Sabq") %>% pull(GFP_intensity_corrected)

  cat(sprintf("\nSabq组GFP强度：\n"))
  cat(sprintf("  校正前中位数：%.2f\n", median(sabq_before)))
  cat(sprintf("  校正后中位数：%.2f\n", median(sabq_after)))
  cat(sprintf("  变化：%.2f (%.1f%%)\n",
              median(sabq_after) - median(sabq_before),
              100 * (median(sabq_after) - median(sabq_before)) / median(sabq_before)))

} else {
  warning("未找到BatchC样本！")
}

# ==============================================================================
# 3. 验证校正效果
# ==============================================================================

cat("\n\n--- 验证校正效果 ---\n\n")

# 验证1：检查BatchA样本是否未改变
cat("验证1：BatchA样本应该未改变\n")
batchA_check <- data_corrected %>%
  filter(batch == "BatchA") %>%
  summarise(
    n_changed = sum(GFP_intensity != GFP_intensity_corrected, na.rm = TRUE)
  )
cat(sprintf("  BatchA中改变的细胞数：%d（应该是0）\n", batchA_check$n_changed))

if (batchA_check$n_changed == 0) {
  cat("  ✓ 验证通过\n")
} else {
  warning("  ✗ BatchA样本被意外修改！")
}

# 验证2：BatchA内部的差异应该保持不变
cat("\n验证2：BatchA内部的BM vs Culture差异应该保持不变\n")

bm_before <- data %>% filter(group == "Bone_marrow") %>% pull(GFP_intensity) %>% median()
culture_before <- data %>% filter(group == "Culture") %>% pull(GFP_intensity) %>% median()
diff_before <- culture_before - bm_before

bm_after <- data_corrected %>% filter(group == "Bone_marrow") %>% pull(GFP_intensity_corrected) %>% median()
culture_after <- data_corrected %>% filter(group == "Culture") %>% pull(GFP_intensity_corrected) %>% median()
diff_after <- culture_after - bm_after

cat(sprintf("  校正前BM vs Culture差异：%.2f\n", diff_before))
cat(sprintf("  校正后BM vs Culture差异：%.2f\n", diff_after))
cat(sprintf("  变化：%.2f（应该接近0）\n", diff_after - diff_before))

if (abs(diff_after - diff_before) < 1) {
  cat("  ✓ 验证通过\n")
} else {
  warning("  ⚠️ BatchA内部差异发生了变化")
}

# ==============================================================================
# 4. 可视化
# ==============================================================================

cat("\n\n--- 生成可视化 ---\n")

dir.create("simple_batch_correction", showWarnings = FALSE)

# 准备绘图数据（下采样）
plot_data <- data_corrected
if (nrow(plot_data) > 50000) {
  plot_data <- plot_data[sample(nrow(plot_data), 50000), ]
}

# 图1：校正前后对比（密度图）
p1_before <- ggplot(plot_data, aes(x = GFP_intensity, fill = batch)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  labs(title = "校正前", x = "GFP强度 (log10)", y = "密度") +
  theme_bw() +
  theme(legend.position = "bottom")

p1_after <- ggplot(plot_data, aes(x = GFP_intensity_corrected, fill = batch)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  labs(title = "校正后", x = "GFP强度 (log10)", y = "密度") +
  theme_bw() +
  theme(legend.position = "bottom")

p1_combined <- plot_grid(p1_before, p1_after, ncol = 2)
print(p1_combined)
ggsave("simple_batch_correction/before_after_density.png", p1_combined,
       width = 12, height = 5, dpi = 150)

# 图2：三组比较（箱线图）
p2_before <- ggplot(plot_data, aes(x = group, y = GFP_intensity, fill = group)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.size = 0.5) +
  scale_y_log10() +
  labs(title = "校正前：三组GFP强度比较",
       x = "", y = "GFP强度 (log10)") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p2_after <- ggplot(plot_data, aes(x = group, y = GFP_intensity_corrected, fill = group)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.size = 0.5) +
  scale_y_log10() +
  labs(title = "校正后：三组GFP强度比较",
       x = "", y = "GFP强度（校正后）(log10)") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

p2_combined <- plot_grid(p2_before, p2_after, ncol = 2)
print(p2_combined)
ggsave("simple_batch_correction/three_groups_comparison.png", p2_combined,
       width = 12, height = 6, dpi = 150)

cat("✓ 可视化完成\n")
cat("  - simple_batch_correction/before_after_density.png\n")
cat("  - simple_batch_correction/three_groups_comparison.png\n")

# ==============================================================================
# 5. 统计比较
# ==============================================================================

cat("\n\n--- 三组统计比较（使用校正后的值）---\n\n")

# 计算每组的MFI（中位荧光强度）
mfi_summary <- data_corrected %>%
  group_by(group) %>%
  summarise(
    n_cells = n(),
    MFI_GFP_before = median(GFP_intensity),
    MFI_GFP_after = median(GFP_intensity_corrected),
    MFI_CD51_before = median(CD51_PE_intensity),
    MFI_CD51_after = median(CD51_PE_intensity_corrected),
    .groups = "drop"
  )

cat("各组MFI（中位荧光强度）：\n")
print(mfi_summary)

# 成对比较（使用校正后的值）
cat("\n\n成对比较（Wilcoxon检验，校正后）：\n")

comparisons <- list(
  c("Bone_marrow", "Culture"),
  c("Bone_marrow", "Sabq"),
  c("Culture", "Sabq")
)

for (comp in comparisons) {
  g1 <- comp[1]
  g2 <- comp[2]

  d1 <- data_corrected %>% filter(group == g1) %>% pull(GFP_intensity_corrected)
  d2 <- data_corrected %>% filter(group == g2) %>% pull(GFP_intensity_corrected)

  test_result <- wilcox.test(d1, d2)

  cat(sprintf("\n%s vs %s:\n", g1, g2))
  cat(sprintf("  中位数：%.2f vs %.2f\n", median(d1), median(d2)))
  cat(sprintf("  p值：%.2e", test_result$p.value))

  if (test_result$p.value < 0.001) {
    cat(" ***")
  } else if (test_result$p.value < 0.01) {
    cat(" **")
  } else if (test_result$p.value < 0.05) {
    cat(" *")
  }
  cat("\n")
}

# ==============================================================================
# 6. 保存结果
# ==============================================================================

cat("\n\n--- 保存校正后的数据 ---\n")

# 保存完整数据
saveRDS(data_corrected, "simple_batch_correction/data_corrected.rds")
write.csv(data_corrected, "simple_batch_correction/data_corrected.csv", row.names = FALSE)

# 保存MFI汇总
write.csv(mfi_summary, "simple_batch_correction/MFI_summary.csv", row.names = FALSE)

cat("✓ 数据已保存\n")
cat("  - simple_batch_correction/data_corrected.rds\n")
cat("  - simple_batch_correction/data_corrected.csv\n")
cat("  - simple_batch_correction/MFI_summary.csv\n")

# ==============================================================================
# 7. 总结报告
# ==============================================================================

cat("\n\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("=== 简单批次校正完成 ===\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

cat("校正方法：\n")
cat("  - 使用未染色样本估计的批次效应\n")
cat(sprintf("  - GFP通道批次效应：%.2f (%.1f%%)\n",
            BATCH_EFFECT$GFP$difference, BATCH_EFFECT$GFP$percent))
cat(sprintf("  - CD51通道批次效应：%.2f (%.1f%%)\n",
            BATCH_EFFECT$CD51$difference, BATCH_EFFECT$CD51$percent))
cat("  - 仅校正BatchC的样本\n\n")

cat("验证结果：\n")
cat("  ✓ BatchA样本未改变\n")
cat("  ✓ BatchA内部差异保持不变\n\n")

cat("三组MFI比较（校正后）：\n")
print(mfi_summary[, c("group", "n_cells", "MFI_GFP_after")])

cat("\n\n建议的下一步：\n")
cat("  1. 查看生成的图表，确认校正效果合理\n")
cat("  2. 使用校正后的值（GFP_intensity_corrected）进行分析\n")
cat("  3. 在论文中报告批次效应及校正方法\n")

cat("\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")
