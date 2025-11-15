# ==============================================================================
# 识别对照样本（未染色、单染管等）
# ==============================================================================
#
# 目的：从WSP文件中识别对照样本，用于批次效应估计
# ==============================================================================

suppressPackageStartupMessages({
  library(CytoML)
  library(flowWorkspace)
  library(flowCore)
})

cat("=== 识别WSP文件中的对照样本 ===\n\n")

# 检查两个WSP文件
wsp_configs <- list(
  BatchA = list(
    file = "20250131 CAR FACS_Miao.wsp",
    batch = "BatchA",
    date = "20250131"
  ),
  BatchC = list(
    file = "20250904 CAR FACS_Kun.wsp",
    batch = "BatchC",
    date = "20250904"
  )
)

all_controls <- list()

for (batch_name in names(wsp_configs)) {
  config <- wsp_configs[[batch_name]]

  cat(paste0(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("%s (%s)\n", batch_name, config$file))
  cat(paste0(rep("=", 80), collapse = ""), "\n\n")

  if (!file.exists(config$file)) {
    cat(sprintf("文件不存在：%s\n\n", config$file))
    next
  }

  # 读取WSP
  ws <- open_flowjo_xml(config$file, sample_names_from = "sampleNode")
  samples <- CytoML::fj_ws_get_samples(ws)

  if (is.null(samples) || nrow(samples) == 0) {
    cat("未找到样本\n\n")
    next
  }

  cat(sprintf("共找到 %d 个样本\n\n", nrow(samples)))

  # 分类样本
  unstained <- c()
  single_stains <- c()
  fmo_amo <- c()
  experimental <- c()

  for (i in seq_len(nrow(samples))) {
    sample_id <- samples$sampleID[i]
    sample_name <- samples$name[i]

    # 获取TUBE NAME
    kw <- tryCatch({
      CytoML::fj_ws_get_keywords(ws, sample_id)
    }, error = function(e) NULL)

    tube_name <- if (!is.null(kw)) {
      as.character(kw[["TUBE NAME"]] %||% kw[["$TUBE"]] %||% kw[["TUBE"]] %||% "")
    } else {
      ""
    }

    # 分类逻辑（基于常见命名模式）
    sample_lower <- tolower(sample_name)
    tube_lower <- tolower(tube_name)

    # 未染色
    if (grepl("unstain|blank|negative|control|ctrl", sample_lower) ||
        grepl("unstain|blank|negative", tube_lower)) {
      unstained <- c(unstained, sample_name)
      type <- "未染色"
    }
    # 单染
    else if (grepl("single|comp|compensation", sample_lower) ||
             grepl("single|comp", tube_lower) ||
             grepl("gfp only|cd51 only|pe only|fitc only", sample_lower)) {
      single_stains <- c(single_stains, sample_name)
      type <- "单染管"
    }
    # FMO/AMO
    else if (grepl("fmo|amo|minus", sample_lower) ||
             grepl("fmo|amo", tube_lower)) {
      fmo_amo <- c(fmo_amo, sample_name)
      type <- "FMO/AMO"
    }
    # 实验样本
    else {
      experimental <- c(experimental, sample_name)
      type <- "实验样本"
    }

    # 打印样本信息
    cat(sprintf("%-3d %-50s [%s]\n", i, sample_name, type))
    if (tube_name != "") {
      cat(sprintf("    TUBE NAME: %s\n", tube_name))
    }
  }

  cat("\n")
  cat(paste0(rep("-", 80), collapse = ""), "\n")
  cat("分类汇总：\n")
  cat(paste0(rep("-", 80), collapse = ""), "\n\n")

  cat(sprintf("未染色样本 (%d):\n", length(unstained)))
  for (s in unstained) cat(sprintf("  - %s\n", s))
  cat("\n")

  cat(sprintf("单染管 (%d):\n", length(single_stains)))
  for (s in single_stains) cat(sprintf("  - %s\n", s))
  cat("\n")

  cat(sprintf("FMO/AMO对照 (%d):\n", length(fmo_amo)))
  for (s in fmo_amo) cat(sprintf("  - %s\n", s))
  cat("\n")

  cat(sprintf("实验样本 (%d):\n", length(experimental)))
  for (s in experimental) cat(sprintf("  - %s\n", s))
  cat("\n\n")

  # 保存结果
  all_controls[[batch_name]] <- list(
    unstained = unstained,
    single_stains = single_stains,
    fmo_amo = fmo_amo,
    experimental = experimental,
    batch = config$batch
  )
}

# ==============================================================================
# 批次矫正建议
# ==============================================================================

cat("\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("=== 批次矫正策略建议 ===\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

batchA_controls <- all_controls[["BatchA"]]
batchC_controls <- all_controls[["BatchC"]]

if (!is.null(batchA_controls) && !is.null(batchC_controls)) {

  # 检查未染色样本
  if (length(batchA_controls$unstained) > 0 && length(batchC_controls$unstained) > 0) {
    cat("✅ 方案1：使用未染色样本估计批次效应（推荐）\n\n")
    cat("   原理：比较两个批次的自发荧光差异\n")
    cat("   BatchA未染色样本：\n")
    for (s in batchA_controls$unstained) cat(sprintf("     - %s\n", s))
    cat("   BatchC未染色样本：\n")
    for (s in batchC_controls$unstained) cat(sprintf("     - %s\n", s))
    cat("\n   步骤：\n")
    cat("   1. 提取未染色样本在各通道的荧光强度\n")
    cat("   2. 计算BatchA和BatchC的中位数差异\n")
    cat("   3. 将差异应用到实验样本\n\n")
  } else {
    cat("⚠️ 未找到两个批次都有的未染色样本\n\n")
  }

  # 检查FMO/AMO
  if (length(batchA_controls$fmo_amo) > 0 || length(batchC_controls$fmo_amo) > 0) {
    cat("✅ 方案2：使用FMO/AMO对照估计批次效应\n\n")
    cat("   原理：在接近目标细胞群的对照中估计批次效应\n")
    if (length(batchA_controls$fmo_amo) > 0) {
      cat("   BatchA FMO/AMO：\n")
      for (s in batchA_controls$fmo_amo) cat(sprintf("     - %s\n", s))
    }
    if (length(batchC_controls$fmo_amo) > 0) {
      cat("   BatchC FMO/AMO：\n")
      for (s in batchC_controls$fmo_amo) cat(sprintf("     - %s\n", s))
    }
    cat("\n   步骤：\n")
    cat("   1. 提取FMO样本在目标门控附近的荧光强度\n")
    cat("   2. 估计批次效应（更精确，因为更接近目标细胞）\n")
    cat("   3. 应用到CD51+GFP+门控\n\n")
  }

  # 检查单染管
  if (length(batchA_controls$single_stains) > 0 && length(batchC_controls$single_stains) > 0) {
    cat("✅ 方案3：使用单染管验证批次效应和补偿一致性\n\n")
    cat("   原理：检查每个通道的批次效应，确保补偿一致\n")
    cat("   BatchA单染管：\n")
    for (s in batchA_controls$single_stains) cat(sprintf("     - %s\n", s))
    cat("   BatchC单染管：\n")
    for (s in batchC_controls$single_stains) cat(sprintf("     - %s\n", s))
    cat("\n   步骤：\n")
    cat("   1. 比较GFP单染在两个批次的强度\n")
    cat("   2. 比较CD51单染在两个批次的强度\n")
    cat("   3. 检查spillover是否一致\n\n")
  }

  cat(paste0(rep("-", 80), collapse = ""), "\n")
  cat("推荐策略：\n")
  cat(paste0(rep("-", 80), collapse = ""), "\n\n")

  if (length(batchA_controls$unstained) > 0 && length(batchC_controls$unstained) > 0) {
    cat("⭐⭐⭐ 优先使用未染色样本\n")
    cat("    - 最稳健的方法\n")
    cat("    - 估计基础的仪器/试剂批次效应\n")
    cat("    - 适用于所有通道\n\n")
  }

  if (length(batchC_controls$fmo_amo) > 0) {
    cat("⭐⭐ 补充使用FMO/AMO对照\n")
    cat("    - 在目标细胞群附近验证批次效应\n")
    cat("    - 更精确地估计CD51+门控的批次效应\n\n")
  }

  cat("建议工作流程：\n")
  cat("  1. 用未染色样本估计整体批次效应\n")
  cat("  2. 用FMO/AMO验证（如果有）\n")
  cat("  3. 用单染管检查补偿一致性（如果有）\n")
  cat("  4. 应用批次矫正到CD51+GFP+门控\n")
  cat("  5. 可靠地比较Sabq vs Culture vs Bone_marrow\n\n")

} else {
  cat("未能获取足够的批次信息\n")
}

cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("识别完成\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

cat("提示：\n")
cat("  如果分类不准确，请手动检查TUBE NAME字段\n")
cat("  或提供样本的准确分类\n")
