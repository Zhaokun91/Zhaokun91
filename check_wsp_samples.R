# ==============================================================================
# 检查WSP文件中的所有样本
# ==============================================================================

suppressPackageStartupMessages({
  library(CytoML)
  library(flowWorkspace)
})

cat("=== 检查WSP文件中的所有样本 ===\n\n")

# 检查所有WSP文件
wsp_files <- c(
  "20250131 CAR FACS_Miao.wsp",
  "20250904 CAR FACS_Kun.wsp"
)

for (wsp_file in wsp_files) {
  if (!file.exists(wsp_file)) {
    cat(sprintf("跳过不存在的文件：%s\n\n", wsp_file))
    next
  }

  cat(paste0(rep("-", 80), collapse = ""), "\n")
  cat(sprintf("WSP文件：%s\n", wsp_file))
  cat(paste0(rep("-", 80), collapse = ""), "\n\n")

  # 读取WSP
  ws <- open_flowjo_xml(wsp_file, sample_names_from = "sampleNode")

  # 获取所有样本
  samples <- CytoML::fj_ws_get_samples(ws)

  if (is.null(samples) || nrow(samples) == 0) {
    cat("未找到样本\n\n")
    next
  }

  cat(sprintf("共有 %d 个样本：\n\n", nrow(samples)))

  # 对每个样本获取详细信息
  for (i in seq_len(nrow(samples))) {
    sample_id <- samples$sampleID[i]
    sample_name <- samples$name[i]

    # 获取关键字
    kw <- tryCatch({
      CytoML::fj_ws_get_keywords(ws, sample_id)
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(kw)) {
      cat(sprintf("%d. %s\n", i, sample_name))
      next
    }

    # 提取有用信息
    tot <- kw[["$TOT"]] %||% kw[["TOT"]] %||% "N/A"
    tube <- kw[["TUBE NAME"]] %||% kw[["$TUBE"]] %||% kw[["TUBE"]] %||% "N/A"

    cat(sprintf("%d. %s\n", i, sample_name))
    cat(sprintf("   TUBE NAME: %s\n", tube))
    cat(sprintf("   总细胞数: %s\n", tot))
    cat("\n")
  }
}

cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("检查完成\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
