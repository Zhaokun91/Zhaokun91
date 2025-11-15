# ==============================================================================
# 批次效应检测 - 使用对照样本（避免混淆）
# ==============================================================================
#
# 目的：在进行批次矫正之前，先检测是否真的存在批次效应
# 策略：使用对照样本（未染色、FMO等）进行检测，避免生物学差异的干扰
#
# 参考：https://biosurf.org/cyCombine_detect_batch_effects.html
# ==============================================================================

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(CytoML)
  library(tidyverse)
  library(cyCombine)
  library(cowplot)
})

cat("=== 批次效应检测（使用对照样本）===\n\n")

# ==============================================================================
# 配置：定义对照样本
# ==============================================================================

# 对照样本配置
CONTROL_SAMPLES <- list(
  BatchA_Unstained = list(
    wsp = "20250131 CAR FACS_Miao.wsp",
    samples = c("Specimen_001_NS.fcs"),  # NS = No Stain (未染色)
    batch = "BatchA",
    type = "Unstained",
    gate = "/all cells/Single Cells1/Single Cells2"  # 使用较早的门控，确保有足够细胞
  ),

  BatchC_Unstained = list(
    wsp = "20250904 CAR FACS_Kun.wsp",
    samples = c("Specimen_001_NS.fcs"),  # NS = No Stain (未染色)
    batch = "BatchC",
    type = "Unstained",
    gate = "/all cells/Single Cells1/Single Cells2"
  )

  # 如果有FMO/AMO对照样本，可以添加：
  # BatchC_AMO = list(
  #   wsp = "20250904 CAR FACS_Kun.wsp",
  #   samples = c("Specimen_001_1.fcs"),  # 假设这是除GFP外都染的样本
  #   batch = "BatchC",
  #   type = "AMO_GFP",
  #   gate = "/all cells/Single Cells1/Single Cells2/Live cells subset/CD45-TER119-CD31-"
  # )
)

# 通道定义
MARKERS <- c("GFP_intensity", "CD51_PE_intensity")

# 输出目录
OUTPUT_DIR <- "batch_effect_detection"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 辅助函数：智能选择通道
# ==============================================================================

pick_channel <- function(ff, patterns) {
  cn <- colnames(ff)
  mk <- suppressWarnings(markernames(ff))
  if (is.null(mk)) mk <- rep(NA_character_, length(cn))

  for (pattern in patterns) {
    hit <- which(grepl(pattern, mk, ignore.case = TRUE))
    if (length(hit) > 0) return(cn[hit[1]])
  }

  for (pattern in patterns) {
    hit <- which(grepl(pattern, cn, ignore.case = TRUE))
    if (length(hit) > 0) return(cn[hit[1]])
  }

  return(NA_character_)
}

# ==============================================================================
# 提取对照样本数据
# ==============================================================================

cat("=== 提取对照样本数据 ===\n\n")

all_control_data <- data.frame()
sample_summary <- data.frame()

for (ctrl_name in names(CONTROL_SAMPLES)) {
  info <- CONTROL_SAMPLES[[ctrl_name]]

  cat(sprintf("\n处理：%s\n", ctrl_name))
  cat(sprintf("  批次：%s\n", info$batch))
  cat(sprintf("  类型：%s\n", info$type))

  wsp_path <- file.path(".", info$wsp)
  if (!file.exists(wsp_path)) {
    warning(sprintf("  ✗ WSP文件不存在: %s", wsp_path))
    next
  }

  # 读取WSP并创建GatingSet（简化流程）
  tryCatch({
    ws <- open_flowjo_xml(wsp_path, sample_names_from = "sampleNode")

    # 建立文件链接
    link_dir <- file.path("temp_control_links", ctrl_name)
    dir.create(link_dir, recursive = TRUE, showWarnings = FALSE)

    # 简化的文件匹配
    fcs_files <- list.files("new", pattern = "\\.fcs$", recursive = TRUE,
                           full.names = TRUE, ignore.case = TRUE)
    ws_samples <- CytoML::fj_ws_get_samples(ws)

    for (target in info$samples) {
      target_base <- tools::file_path_sans_ext(basename(target))
      matched <- grep(target_base, fcs_files, value = TRUE, ignore.case = TRUE)

      if (length(matched) > 0) {
        wsp_sample <- ws_samples$name[grep(target_base, ws_samples$name, ignore.case = TRUE)[1]]
        if (!is.na(wsp_sample)) {
          link_path <- file.path(link_dir, wsp_sample)
          if (!file.exists(link_path)) {
            file.copy(matched[1], link_path, overwrite = TRUE)
          }
        }
      }
    }

    # 创建GatingSet
    gs <- flowjo_to_gatingset(
      ws,
      name = "Controls",
      path = normalizePath(link_dir, winslash = "/", mustWork = FALSE),
      leaf.bool = FALSE,
      skip_faulty_gate = TRUE,
      sample_names_from = "sampleNode"
    )

    if (is.null(gs) || length(gs) == 0) {
      warning(sprintf("  ✗ 未创建GatingSet"))
      next
    }

    # 提取数据
    for (i in seq_along(gs)) {
      sample_name <- sampleNames(gs)[i]

      # 提取门控数据
      ff_gate <- gh_pop_get_data(gs[[i]], info$gate)
      ex <- exprs(ff_gate)

      if (nrow(ex) == 0) {
        warning(sprintf("  ✗ 样本%s在门控%s中无细胞", sample_name, info$gate))
        next
      }

      # 下采样（避免内存问题）
      if (nrow(ex) > 10000) {
        ex <- ex[sample(nrow(ex), 10000), ]
      }

      # 识别通道
      gfp_ch <- pick_channel(ff_gate, c("GFP", "FITC", "FL1"))
      cd51_ch <- pick_channel(ff_gate, c("CD51", "PE", "FL2"))

      if (is.na(gfp_ch) || is.na(cd51_ch)) {
        warning(sprintf("  ✗ 样本%s未找到GFP或CD51通道", sample_name))
        next
      }

      # 创建数据框
      sample_df <- data.frame(
        id = seq_len(nrow(ex)),
        sample = sample_name,
        batch = info$batch,
        control_type = info$type,
        GFP_intensity = as.numeric(ex[, gfp_ch]),
        CD51_PE_intensity = as.numeric(ex[, cd51_ch]),
        stringsAsFactors = FALSE
      )

      all_control_data <- rbind(all_control_data, sample_df)

      summary_row <- data.frame(
        control = ctrl_name,
        batch = info$batch,
        sample = sample_name,
        n_cells = nrow(ex),
        stringsAsFactors = FALSE
      )

      sample_summary <- rbind(sample_summary, summary_row)

      cat(sprintf("  ✓ %s: %d个细胞\n", sample_name, nrow(ex)))
    }

  }, error = function(e) {
    warning(sprintf("  ✗ 处理失败: %s", e$message))
  })
}

# 清理临时文件
unlink("temp_control_links", recursive = TRUE)

if (nrow(all_control_data) == 0) {
  stop("未提取到任何对照样本数据！请检查样本名称配置。")
}

cat(sprintf("\n总共提取%d个细胞的对照数据\n", nrow(all_control_data)))
cat("\n样本汇总：\n")
print(sample_summary)

# ==============================================================================
# 数据准备（cyCombine格式）
# ==============================================================================

cat("\n\n=== 准备cyCombine格式数据 ===\n")

# 清理数据
df <- all_control_data %>%
  filter(!is.na(GFP_intensity), !is.na(CD51_PE_intensity)) %>%
  filter(GFP_intensity > 0, CD51_PE_intensity > 0) %>%
  mutate(
    batch = as.character(batch),
    condition = as.character(control_type)  # 对照类型作为condition
  ) %>%
  as_tibble()

cat(sprintf("清理后数据：%d个细胞\n", nrow(df)))

# asinh变换
cofactor <- 150
df <- df %>%
  mutate(
    GFP_intensity = asinh(GFP_intensity / cofactor),
    CD51_PE_intensity = asinh(CD51_PE_intensity / cofactor)
  )

cat("数据范围（asinh变换后）：\n")
cat(sprintf("  GFP: [%.3f, %.3f]\n", min(df$GFP_intensity), max(df$GFP_intensity)))
cat(sprintf("  CD51: [%.3f, %.3f]\n", min(df$CD51_PE_intensity), max(df$CD51_PE_intensity)))

# ==============================================================================
# 批次效应检测
# ==============================================================================

cat("\n\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("=== 批次效应检测 ===\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

# 检查批次分布
cat("批次分布：\n")
print(table(df$batch))
cat("\n")

# 1. 快速检测
cat("--- 方法1：快速批次效应检测 ---\n\n")
tryCatch({
  detect_batch_effect_express(
    df,
    downsample = min(5000, nrow(df)),
    out_dir = file.path(OUTPUT_DIR, "express"),
    markers = c("GFP_intensity", "CD51_PE_intensity"),
    batch_col = "batch"
  )
  cat("✓ 快速检测完成\n")
  cat(sprintf("  结果保存在: %s/express/\n", OUTPUT_DIR))
  cat("  请查看以下图表：\n")
  cat("    - EMD_per_marker.png: 每个marker的批次间距离\n")
  cat("    - Marker_distributions.png: 各批次的marker分布\n")
  cat("    - MDS_plot.png: 基于中位数表达的样本分布\n")
}, error = function(e) {
  cat(sprintf("✗ 快速检测失败: %s\n", e$message))
})

# 2. 手动可视化（作为备选）
cat("\n--- 方法2：手动批次效应可视化 ---\n\n")

tryCatch({
  # 密度图对比
  p1 <- ggplot(df, aes(x = GFP_intensity, fill = batch)) +
    geom_density(alpha = 0.5) +
    labs(title = "GFP强度分布（对照样本）",
         subtitle = "如果两个批次的分布明显不同，说明存在批次效应",
         x = "GFP强度 (asinh变换)") +
    theme_bw()

  p2 <- ggplot(df, aes(x = CD51_PE_intensity, fill = batch)) +
    geom_density(alpha = 0.5) +
    labs(title = "CD51强度分布（对照样本）",
         x = "CD51强度 (asinh变换)") +
    theme_bw()

  # 箱线图
  p3 <- ggplot(df, aes(x = batch, y = GFP_intensity, fill = batch)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.5) +
    labs(title = "GFP强度（按批次）",
         y = "GFP强度") +
    theme_bw() +
    theme(legend.position = "none")

  p4 <- ggplot(df, aes(x = batch, y = CD51_PE_intensity, fill = batch)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.5) +
    labs(title = "CD51强度（按批次）",
         y = "CD51强度") +
    theme_bw() +
    theme(legend.position = "none")

  combined <- plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2)
  print(combined)
  ggsave(file.path(OUTPUT_DIR, "manual_batch_comparison.png"),
         combined, width = 12, height = 10, dpi = 150)

  cat("✓ 手动可视化完成\n")
  cat(sprintf("  图表保存为: %s/manual_batch_comparison.png\n", OUTPUT_DIR))

}, error = function(e) {
  cat(sprintf("✗ 手动可视化失败: %s\n", e$message))
})

# 3. 统计检验
cat("\n--- 方法3：统计检验 ---\n\n")

batches <- unique(df$batch)
if (length(batches) == 2) {
  batch1_gfp <- df %>% filter(batch == batches[1]) %>% pull(GFP_intensity)
  batch2_gfp <- df %>% filter(batch == batches[2]) %>% pull(GFP_intensity)

  batch1_cd51 <- df %>% filter(batch == batches[1]) %>% pull(CD51_PE_intensity)
  batch2_cd51 <- df %>% filter(batch == batches[2]) %>% pull(CD51_PE_intensity)

  # GFP通道
  gfp_test <- wilcox.test(batch1_gfp, batch2_gfp)
  cat("GFP通道批次差异检验：\n")
  cat(sprintf("  %s中位数: %.4f\n", batches[1], median(batch1_gfp)))
  cat(sprintf("  %s中位数: %.4f\n", batches[2], median(batch2_gfp)))
  cat(sprintf("  差异: %.4f (%.1f%%)\n",
              median(batch2_gfp) - median(batch1_gfp),
              100 * (median(batch2_gfp) - median(batch1_gfp)) / median(batch1_gfp)))
  cat(sprintf("  Wilcoxon检验 p值: %.2e\n", gfp_test$p.value))

  if (gfp_test$p.value < 0.01) {
    cat("  *** GFP通道存在显著批次效应 (p < 0.01)\n")
  } else if (gfp_test$p.value < 0.05) {
    cat("  ** GFP通道存在批次效应 (p < 0.05)\n")
  } else {
    cat("  ✓ GFP通道无显著批次效应 (p >= 0.05)\n")
  }

  cat("\n")

  # CD51通道
  cd51_test <- wilcox.test(batch1_cd51, batch2_cd51)
  cat("CD51通道批次差异检验：\n")
  cat(sprintf("  %s中位数: %.4f\n", batches[1], median(batch1_cd51)))
  cat(sprintf("  %s中位数: %.4f\n", batches[2], median(batch2_cd51)))
  cat(sprintf("  差异: %.4f (%.1f%%)\n",
              median(batch2_cd51) - median(batch1_cd51),
              100 * (median(batch2_cd51) - median(batch1_cd51)) / median(batch1_cd51)))
  cat(sprintf("  Wilcoxon检验 p值: %.2e\n", cd51_test$p.value))

  if (cd51_test$p.value < 0.01) {
    cat("  *** CD51通道存在显著批次效应 (p < 0.01)\n")
  } else if (cd51_test$p.value < 0.05) {
    cat("  ** CD51通道存在批次效应 (p < 0.05)\n")
  } else {
    cat("  ✓ CD51通道无显著批次效应 (p >= 0.05)\n")
  }
}

# ==============================================================================
# 总结和建议
# ==============================================================================

cat("\n\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("=== 批次效应检测总结 ===\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

cat("检测方法：\n")
cat("  ✓ 使用对照样本（未染色或FMO）检测批次效应\n")
cat("  ✓ 避免了生物学差异的干扰\n")
cat("  ✓ 估计的是纯粹的技术性批次效应\n\n")

cat("生成的文件：\n")
cat(sprintf("  - %s/express/: cyCombine快速检测结果\n", OUTPUT_DIR))
cat(sprintf("  - %s/manual_batch_comparison.png: 手动可视化\n", OUTPUT_DIR))
cat("\n")

cat("如何判断是否需要批次矫正：\n\n")

cat("1. 查看密度图和箱线图：\n")
cat("   - 如果两个批次的分布明显分离 → 需要批次矫正\n")
cat("   - 如果两个批次高度重叠 → 可能不需要批次矫正\n\n")

cat("2. 查看统计检验：\n")
cat("   - 如果p < 0.01 且差异 > 10% → 强烈建议批次矫正\n")
cat("   - 如果p < 0.05 且差异 > 5% → 建议批次矫正\n")
cat("   - 如果p >= 0.05 或差异 < 5% → 可能不需要批次矫正\n\n")

cat("3. 查看EMD图（如果成功生成）：\n")
cat("   - EMD值越高，批次效应越强\n")
cat("   - EMD > 0.3 通常认为存在显著批次效应\n\n")

cat("下一步建议：\n\n")

if (exists("gfp_test") && gfp_test$p.value < 0.05) {
  cat("⚠️ 检测到显著批次效应！\n\n")
  cat("建议操作：\n")
  cat("  1. 对实验样本进行批次矫正\n")
  cat("  2. 使用对照样本估计的批次效应参数\n")
  cat("  3. 矫正后重新比较三组（Sabq vs Culture vs Bone_marrow）\n")
} else {
  cat("✓ 未检测到显著批次效应\n\n")
  cat("建议操作：\n")
  cat("  1. 可以直接比较三组，无需批次矫正\n")
  cat("  2. 但仍需注意在结果中说明批次情况\n")
  cat("  3. 或进行保守的批次矫正以确保稳健性\n")
}

cat("\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("检测完成！请查看生成的图表。\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
