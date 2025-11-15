# ==============================================================================
# 检查所有门控的细胞数 - 用于批次效应估计
# ==============================================================================
#
# 目的：检查BatchC的Bone_marrow样本在各个门控中的细胞数
#       以确定哪个门控适合用于批次效应估计
# ==============================================================================

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(CytoML)
  library(tidyverse)
})

cat("=== 检查所有门控层级的细胞数 ===\n\n")

# 配置：需要检查的样本
SAMPLES_TO_CHECK <- list(
  BatchA_Bone_marrow = list(
    wsp = "20250131 CAR FACS_Miao.wsp",
    samples = c("Specimen_001_B1.fcs", "Specimen_001_B2.fcs", "Specimen_001_B3.fcs"),
    batch = "BatchA",
    group = "Bone_marrow"
  ),

  BatchC_Bone_marrow = list(
    wsp = "20250904 CAR FACS_Kun.wsp",
    samples = c("Specimen_001_1.fcs"),  # 猜测是这个，需要确认
    batch = "BatchC",
    group = "Bone_marrow"
  ),

  BatchC_Sabq = list(
    wsp = "20250904 CAR FACS_Kun.wsp",
    samples = c("Specimen_001_2.fcs", "Specimen_001_3.fcs", "Specimen_001_7.fcs"),
    batch = "BatchC",
    group = "Sabq"
  )
)

# 存储结果
all_results <- list()

for (group_name in names(SAMPLES_TO_CHECK)) {
  info <- SAMPLES_TO_CHECK[[group_name]]

  cat("\n")
  cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("检查：%s\n", group_name))
  cat(paste0(rep("=", 80), collapse = ""), "\n\n")

  wsp_path <- file.path(".", info$wsp)

  if (!file.exists(wsp_path)) {
    cat(sprintf("WSP文件不存在：%s\n", wsp_path))
    next
  }

  # 读取WSP
  cat("读取WSP文件...\n")
  ws <- open_flowjo_xml(wsp_path, sample_names_from = "sampleNode")

  # 建立文件链接
  link_dir <- file.path("temp_links", group_name)
  dir.create(link_dir, recursive = TRUE, showWarnings = FALSE)

  # 简化版的文件匹配
  fcs_files <- list.files("new", pattern = "\\.fcs$", recursive = TRUE,
                          full.names = TRUE, ignore.case = TRUE)

  ws_samples <- CytoML::fj_ws_get_samples(ws)

  for (target_sample in info$samples) {
    # 尝试按名称匹配
    target_base <- tools::file_path_sans_ext(basename(target_sample))
    matched <- grep(target_base, fcs_files, value = TRUE, ignore.case = TRUE)

    if (length(matched) > 0) {
      # 找到WSP中对应的样本
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
  cat("创建GatingSet...\n")
  tryCatch({
    gs <- flowjo_to_gatingset(
      ws,
      name = "All Samples",
      path = normalizePath(link_dir, winslash = "/", mustWork = FALSE),
      leaf.bool = FALSE,
      skip_faulty_gate = TRUE,
      sample_names_from = "sampleNode"
    )

    if (is.null(gs) || length(gs) == 0) {
      cat("未创建GatingSet\n")
      next
    }

    cat(sprintf("成功创建GatingSet，共%d个样本\n\n", length(gs)))

    # 获取所有门控节点
    nodes <- gs_get_pop_paths(gs)

    cat("检查各个门控的细胞数：\n\n")

    for (node in nodes) {
      # 跳过根节点
      if (node == "root") next

      # 统计所有样本在该门控的细胞数
      total_cells <- 0

      for (i in seq_along(gs)) {
        tryCatch({
          n <- gh_pop_get_count(gs[[i]], node)
          total_cells <- total_cells + n
        }, error = function(e) {
          # 忽略错误
        })
      }

      # 简化节点名显示
      node_short <- basename(node)

      cat(sprintf("  %-40s: %6d 个细胞\n", node_short, total_cells))

      # 保存结果
      result_key <- paste(group_name, node_short, sep = "|")
      all_results[[result_key]] <- list(
        group = group_name,
        node = node_short,
        node_full = node,
        cells = total_cells,
        batch = info$batch,
        bio_group = info$group
      )
    }

  }, error = function(e) {
    cat(sprintf("错误：%s\n", e$message))
  })
}

# ==============================================================================
# 汇总分析
# ==============================================================================

cat("\n\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("=== 批次效应估计可行性分析 ===\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

# 转换为数据框
results_df <- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    group = x$group,
    node = x$node,
    cells = x$cells,
    batch = x$batch,
    bio_group = x$bio_group,
    stringsAsFactors = FALSE
  )
}))

if (nrow(results_df) > 0) {
  # 重点关注Bone_marrow样本
  bm_results <- results_df %>%
    filter(bio_group == "Bone_marrow") %>%
    arrange(node, batch)

  cat("Bone_marrow样本在各门控的细胞数对比：\n\n")

  # 宽格式显示
  bm_wide <- bm_results %>%
    select(node, batch, cells) %>%
    pivot_wider(names_from = batch, values_from = cells, values_fill = 0)

  print(bm_wide, n = 100)

  cat("\n\n建议：\n\n")

  # 找到BatchA和BatchC都有足够细胞的门控
  suitable_gates <- bm_wide %>%
    filter(BatchA > 100 & BatchC > 100) %>%
    pull(node)

  if (length(suitable_gates) > 0) {
    cat("✓ 以下门控可用于批次效应估计（BatchA和BatchC都有>100个细胞）：\n")
    for (gate in suitable_gates) {
      cat(sprintf("  - %s\n", gate))
    }
    cat("\n推荐策略：\n")
    cat("  1. 使用上述门控中的Bone_marrow样本估计批次效应\n")
    cat("  2. 应用批次矫正到CD51+GFP+门控的数据\n")
    cat("  3. 可靠地比较Sabq vs Culture vs Bone_marrow\n")
  } else {
    cat("✗ 没有找到两个批次都有足够细胞(>100)的门控\n\n")
    cat("推荐策略：\n")
    cat("  1. 仅报告BatchA内部的比较（Bone_marrow vs Culture）\n")
    cat("  2. Sabq的比较作为探索性结果，明确说明可能受批次效应影响\n")
    cat("  3. 建议补充实验：在BatchC中也收集GFP+的Bone_marrow样本\n")
  }

  # 保存结果
  write.csv(results_df, "gate_cell_counts.csv", row.names = FALSE)
  cat("\n\n详细结果已保存：gate_cell_counts.csv\n")

} else {
  cat("未获取到有效数据\n")
}

cat("\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("检查完成\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")

# 清理临时文件
unlink("temp_links", recursive = TRUE)
