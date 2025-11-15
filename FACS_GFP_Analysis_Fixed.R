# ==============================================================================
# FACS GFP+CD51+ 分析脚本 - 修复版本
# ==============================================================================
# 修复内容：
#   1. 批次矫正数据准备阶段的NA值处理
#   2. batch列类型问题
#   3. 归一化值生成和匹配逻辑
#   4. 改进错误处理和调试信息
# ==============================================================================

# --- 0. 控制面板（所有参数都在这里设置） ---
# ============================================

# 路径设置
FCS_ROOT     <- "new"                    # FCS文件根目录
WSP_DIR      <- "."                      # WSP文件目录
OUTPUT_DIR   <- "out_exports"            # 输出目录
METADATA_CSV <- "metadata.csv"           # 元数据文件（可选）

# 并行处理
USE_PARALLEL <- TRUE                     # 是否使用多核并行
MAX_CORES    <- NULL                     # 最大核心数（NULL=自动检测-1）

# 导出控制
EXPORT_DATA  <- TRUE                     # 是否导出结果
EXPORT_EXCEL <- TRUE                     # 导出Excel文件
EXPORT_CSV   <- TRUE                     # 导出CSV文件

# 数据过滤（节省内存）
MAX_CELLS_PER_SAMPLE <- 1000000         # 每个样本最大细胞数（防止内存溢出）

# 批次归一化
DO_BATCH_CORRECTION <- TRUE              # 是否进行批次归一化

# ==============================================================================

# --- 1. 加载必需的R包 ---
# ==========================

cat("=== 加载R包 ===\n")

suppressPackageStartupMessages({
  library(flowCore)       # FCS文件读取
  library(flowWorkspace)  # GatingSet操作
  library(CytoML)         # FlowJo WSP文件解析
  library(tidyverse)      # 数据处理
  library(openxlsx)       # Excel导出
  library(data.table)     # 高效数据处理
  if (DO_BATCH_CORRECTION) {
    library(cyCombine)    # 批次归一化
    library(cowplot)      # 用于组合图表
  }
})

cat("✓ R包加载完成\n\n")

# --- 2. 设置并行处理 ---
# ========================

if (USE_PARALLEL) {
  CORES <- parallel::detectCores() - 1
  if (!is.null(MAX_CORES)) CORES <- min(CORES, MAX_CORES)
  CORES <- max(1, CORES)
  cat("并行处理：使用", CORES, "个CPU核心\n")
} else {
  CORES <- 1
  cat("单线程模式\n")
}

# 创建输出目录
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 3. 定义样本信息 ---
# ========================

cat("\n=== 定义样本信息 ===\n")

# 样本定义：这里定义了所有实验样本的信息
SAMPLE_INFO <- list(
  # 批次A：Bone_marrow
  Bone_marrow = list(
    folder   = "new/20250131",
    samples  = c("Specimen_001_B1.fcs", "Specimen_001_B2.fcs", "Specimen_001_B3.fcs"),
    wsp_file = "20250131 CAR FACS_Miao.wsp",
    batch    = "BatchA",
    target_gate = "/all cells/Single Cells1/Single Cells2/Live cells subset/CD45-TER119-CD31-/CD51+GFP+"
  ),

  # 批次A：Culture
  Culture = list(
    folder   = "new/20250131",
    samples  = c("Specimen_001_C1.fcs", "Specimen_001_C2.fcs",
                 "Specimen_001_C3.fcs", "Specimen_001_C4.fcs"),
    wsp_file = "20250131 CAR FACS_Miao.wsp",
    batch    = "BatchA",
    target_gate = "/all cells/Single Cells1/Single Cells2/Live cells subset/CD45-TER119-CD31-/CD51+GFP+"
  ),

  # 批次C：Sabq
  Sabq_BatchC = list(
    folder   = "new/20250904",
    samples  = c("Specimen_001_2.fcs", "Specimen_001_3.fcs", "Specimen_001_7.fcs"),
    wsp_file = "20250904 CAR FACS_Kun.wsp",
    batch    = "BatchC",
    target_gate = "/all cells/Single Cells1/Single Cells2/Live cells subset/CD45-TER119-CD31-/CD51+GFP+"
  )
)

# 打印样本信息
for (group_name in names(SAMPLE_INFO)) {
  info <- SAMPLE_INFO[[group_name]]
  cat(sprintf("  %s: %d个样本, 批次=%s, WSP=%s\n",
              group_name, length(info$samples), info$batch, info$wsp_file))
}

cat(sprintf("\n总计：%d组，%d个样本\n\n",
            length(SAMPLE_INFO),
            sum(sapply(SAMPLE_INFO, function(x) length(x$samples)))))

# --- 4. 辅助函数定义 ---
# ========================

#' 建立WSP文件名到FCS文件的映射链接
build_wsp_fcs_links <- function(ws, source_root = "new", out_dir = NULL, wsp_name = NULL) {

  if (is.null(out_dir)) {
    if (!is.null(wsp_name)) {
      wsp_basename <- tools::file_path_sans_ext(basename(wsp_name))
    } else {
      wsp_basename <- paste0("wsp_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    }
    out_dir <- file.path("wsp_links", wsp_basename)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # 获取WSP中的样本信息
  ws_samples <- CytoML::fj_ws_get_samples(ws)
  if (is.null(ws_samples) || nrow(ws_samples) == 0) return(out_dir)

  # 索引所有本地FCS文件
  cat("    索引本地FCS文件...\n")
  fcs_files <- list.files(source_root, pattern = "\\.fcs$",
                         recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

  if (length(fcs_files) == 0) {
    warning("未找到FCS文件")
    return(out_dir)
  }

  # 为每个FCS文件提取关键信息
  fcs_index <- list()
  for (fpath in fcs_files) {
    tryCatch({
      kw <- suppressWarnings(keyword(read.FCS(fpath, transformation = FALSE, alter.names = TRUE)))
      fcs_index[[length(fcs_index) + 1]] <- list(
        path = fpath,
        basename = basename(fpath),
        tot  = suppressWarnings(as.numeric(kw[["$TOT"]] %||% kw[["TOT"]] %||% NA)),
        tube = as.character(kw[["TUBE NAME"]] %||% kw[["$TUBE"]] %||% kw[["TUBE"]] %||% NA)
      )
    }, error = function(e) {
      # 静默跳过无法读取的文件
    })
  }

  if (length(fcs_index) == 0) {
    warning("无法读取任何FCS文件")
    return(out_dir)
  }

  cat(sprintf("    找到%d个FCS文件\n", length(fcs_index)))

  # 匹配WSP样本到本地文件
  matches <- 0
  for (i in seq_len(nrow(ws_samples))) {
    wsp_id <- ws_samples$sampleID[i]
    wsp_name <- as.character(ws_samples$name[i])

    # 获取WSP中该样本的关键字
    wsp_kw <- tryCatch({
      CytoML::fj_ws_get_keywords(ws, wsp_id)
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(wsp_kw)) next

    wsp_tot <- suppressWarnings(as.numeric(wsp_kw[["$TOT"]] %||% wsp_kw[["TOT"]] %||% NA))
    wsp_tube <- as.character(wsp_kw[["TUBE NAME"]] %||% wsp_kw[["$TUBE"]] %||% wsp_kw[["TUBE"]] %||% NA)

    # 策略1：通过$TOT和TUBE NAME匹配
    matched_file <- NULL
    if (!is.na(wsp_tot) && !is.na(wsp_tube)) {
      for (idx in seq_along(fcs_index)) {
        fcs_tot <- fcs_index[[idx]]$tot
        fcs_tube <- fcs_index[[idx]]$tube
        if (!is.na(fcs_tot) && !is.na(fcs_tube)) {
          if (abs(fcs_tot - wsp_tot) / max(1, wsp_tot) < 0.02 && fcs_tube == wsp_tube) {
            matched_file <- fcs_index[[idx]]$path
            break
          }
        }
      }
    }

    # 策略2：仅通过$TOT匹配
    if (is.null(matched_file) && !is.na(wsp_tot)) {
      for (idx in seq_along(fcs_index)) {
        fcs_tot <- fcs_index[[idx]]$tot
        if (!is.na(fcs_tot) && abs(fcs_tot - wsp_tot) / max(1, wsp_tot) < 0.02) {
          matched_file <- fcs_index[[idx]]$path
          break
        }
      }
    }

    # 创建符号链接或复制文件
    if (!is.null(matched_file)) {
      tryCatch({
        target_abs <- normalizePath(matched_file, winslash = "/", mustWork = TRUE)
        link_path <- file.path(out_dir, wsp_name)

        if (file.exists(link_path)) unlink(link_path, force = TRUE)

        link_success <- file.symlink(target_abs, link_path)

        if (!link_success || !file.exists(link_path)) {
          file.copy(target_abs, link_path, overwrite = TRUE)
        }

        if (file.exists(link_path)) {
          matches <- matches + 1
        }
      }, error = function(e) {
        # 静默跳过失败的链接
      })
    }
  }

  cat(sprintf("    成功匹配%d/%d个样本\n", matches, nrow(ws_samples)))

  # 清理损坏的链接
  entries <- list.files(out_dir, full.names = TRUE)
  broken <- entries[!file.exists(entries)]
  if (length(broken) > 0) unlink(broken, force = TRUE)

  return(out_dir)
}

#' 从flowFrame中智能选择通道
pick_channel <- function(ff, patterns) {
  cn <- colnames(ff)
  mk <- suppressWarnings(markernames(ff))
  if (is.null(mk)) mk <- rep(NA_character_, length(cn))

  # 首先尝试在marker中匹配
  for (pattern in patterns) {
    hit <- which(grepl(pattern, mk, ignore.case = TRUE))
    if (length(hit) > 0) return(cn[hit[1]])
  }

  # 然后在通道名中匹配
  for (pattern in patterns) {
    hit <- which(grepl(pattern, cn, ignore.case = TRUE))
    if (length(hit) > 0) return(cn[hit[1]])
  }

  return(NA_character_)
}

# --- 5. 读取WSP文件并创建GatingSet ---
# =======================================

cat("\n=== 读取WSP文件并应用门控 ===\n")

# 存储所有GatingSet
ALL_GATING_SETS <- list()

# 按组处理
for (group_name in names(SAMPLE_INFO)) {
  info <- SAMPLE_INFO[[group_name]]

  cat(sprintf("\n处理组：%s\n", group_name))
  cat(sprintf("  批次：%s\n", info$batch))
  cat(sprintf("  WSP文件：%s\n", info$wsp_file))
  cat(sprintf("  期望样本：%s\n", paste(info$samples, collapse = ", ")))

  # 检查WSP文件是否存在
  wsp_path <- file.path(WSP_DIR, info$wsp_file)
  if (!file.exists(wsp_path)) {
    warning(sprintf("  ✗ WSP文件不存在: %s", wsp_path))
    next
  }

  # 读取WSP文件
  cat("  读取WSP文件...\n")
  ws <- open_flowjo_xml(wsp_path, sample_names_from = "sampleNode")

  # 建立WSP到FCS的映射链接
  cat("  建立文件名映射...\n")
  link_dir <- build_wsp_fcs_links(ws, source_root = FCS_ROOT, wsp_name = info$wsp_file)

  # 创建GatingSet
  cat("  创建GatingSet...\n")
  tryCatch({
    gs <- flowjo_to_gatingset(
      ws,
      name = "All Samples",
      path = normalizePath(link_dir, winslash = "/", mustWork = FALSE),
      leaf.bool = FALSE,
      skip_faulty_gate = TRUE,
      sample_names_from = "sampleNode",
      mc.cores = CORES
    )

    if (!is.null(gs) && length(gs) > 0) {

      # 筛选出指定的样本
      cat("  筛选指定样本...\n")

      all_sample_names <- sampleNames(gs)
      cat(sprintf("    GatingSet中共有%d个样本\n", length(all_sample_names)))

      target_patterns <- gsub("\\.fcs$", "", info$samples, ignore.case = TRUE)

      matched_indices <- c()
      matched_names <- c()

      for (pattern in target_patterns) {
        idx <- grep(pattern, all_sample_names, ignore.case = TRUE)
        if (length(idx) > 0) {
          matched_indices <- c(matched_indices, idx)
          matched_names <- c(matched_names, all_sample_names[idx])
          cat(sprintf("    ✓ 匹配到: %s → %s\n", pattern, all_sample_names[idx[1]]))
        } else {
          cat(sprintf("    ✗ 未匹配: %s\n", pattern))
        }
      }

      matched_indices <- unique(matched_indices)

      if (length(matched_indices) == 0) {
        warning(sprintf("  ✗ 未找到任何匹配的样本"))
        next
      }

      cat(sprintf("    最终选中%d个样本\n", length(matched_indices)))

      # 子集化GatingSet
      gs <- gs[matched_indices]

      # 添加组和批次信息到pData
      actual_group <- ifelse(grepl("Sabq", group_name), "Sabq", group_name)
      pData(gs)$group <- actual_group
      pData(gs)$batch <- info$batch
      pData(gs)$display_name <- basename(sampleNames(gs))

      ALL_GATING_SETS[[group_name]] <- gs

      cat(sprintf("  ✓ 成功：%d个样本，%d个门控节点\n",
                  length(gs), length(gs_get_pop_paths(gs))))

    } else {
      warning(sprintf("  ✗ 创建GatingSet失败"))
    }

  }, error = function(e) {
    warning(sprintf("  ✗ 错误：%s", e$message))
  })
}

cat(sprintf("\n总共成功加载%d个组的GatingSet\n", length(ALL_GATING_SETS)))

# --- 6. 提取目标门控数据（每个细胞的GFP强度） ---
# ============================================================

cat("\n=== 提取目标门控数据（单细胞水平）===\n")

# 存储提取的细胞数据
ALL_CELL_DATA <- data.frame()

# 样本汇总统计
SAMPLE_SUMMARY <- data.frame()

for (group_name in names(ALL_GATING_SETS)) {
  gs <- ALL_GATING_SETS[[group_name]]
  info <- SAMPLE_INFO[[group_name]]

  TARGET_GATE <- info$target_gate

  cat(sprintf("\n处理组：%s\n", group_name))
  cat(sprintf("目标门控：%s\n", TARGET_GATE))

  # 检查目标门控是否存在
  nodes <- gs_get_pop_paths(gs)
  if (!(TARGET_GATE %in% nodes)) {
    warning(sprintf("  ✗ 未找到目标门控：%s", TARGET_GATE))
    next
  }

  # 对每个样本提取数据
  for (i in seq_along(gs)) {
    sample_name <- pData(gs)$display_name[i]

    tryCatch({
      # 提取门控后的细胞
      ff_gate <- gh_pop_get_data(gs[[i]], TARGET_GATE)
      ex <- exprs(ff_gate)

      # 限制细胞数
      if (nrow(ex) > MAX_CELLS_PER_SAMPLE) {
        cat(sprintf("    警告：样本%s有%d个细胞，超过限制%d，进行随机采样\n",
                    sample_name, nrow(ex), MAX_CELLS_PER_SAMPLE))
        ex <- ex[sample(nrow(ex), MAX_CELLS_PER_SAMPLE), ]
      }

      # 智能识别GFP和CD51/PE通道
      gfp_ch <- pick_channel(ff_gate, c("GFP", "FITC", "FL1"))
      cd51_ch <- pick_channel(ff_gate, c("CD51", "PE", "FL2"))

      if (is.na(gfp_ch)) {
        warning(sprintf("    ✗ 样本%s：未找到GFP通道", sample_name))
        next
      }

      if (is.na(cd51_ch)) {
        cat(sprintf("    ⚠ 样本%s：未找到CD51/PE通道，仅提取GFP数据\n", sample_name))
      }

      # 提取荧光强度
      gfp_int <- ex[, gfp_ch]
      cd51_int <- if (!is.na(cd51_ch)) ex[, cd51_ch] else rep(NA_real_, length(gfp_int))

      # 创建数据框
      actual_group <- ifelse(grepl("Sabq", group_name), "Sabq", group_name)
      sample_df <- data.frame(
        group = actual_group,
        batch = info$batch,
        sample = sample_name,
        cell_id = seq_len(length(gfp_int)),
        GFP_intensity = as.numeric(gfp_int),
        CD51_PE_intensity = as.numeric(cd51_int),
        stringsAsFactors = FALSE
      )

      ALL_CELL_DATA <- rbind(ALL_CELL_DATA, sample_df)

      # 简单汇总
      actual_group <- ifelse(grepl("Sabq", group_name), "Sabq", group_name)
      summary_row <- data.frame(
        group = actual_group,
        batch = info$batch,
        sample = sample_name,
        n_cells = length(gfp_int),
        stringsAsFactors = FALSE
      )

      SAMPLE_SUMMARY <- rbind(SAMPLE_SUMMARY, summary_row)

      cat(sprintf("  ✓ %s: 提取了%d个细胞的GFP强度\n", sample_name, length(gfp_int)))

    }, error = function(e) {
      warning(sprintf("  ✗ 样本%s处理失败：%s", sample_name, e$message))
    })
  }
}

cat(sprintf("\n总共提取%d个细胞的GFP荧光强度数据\n", nrow(ALL_CELL_DATA)))

# --- 7. cyCombine 批次效应检测和矫正 (修复版本) ---
# ====================================================

if (DO_BATCH_CORRECTION && nrow(ALL_CELL_DATA) > 0) {
  cat("\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  cat("=== cyCombine 批次效应检测和矫正 (修复版本) ===\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")

  # --- 7.1 检查数据 ---
  cat("\n--- 检查数据 ---\n")
  cat(sprintf("总细胞数：%d\n", nrow(ALL_CELL_DATA)))
  cat("\n批次分布：\n")
  print(table(ALL_CELL_DATA$batch))
  cat("\n组分布：\n")
  print(table(ALL_CELL_DATA$group))

  # --- 7.2 准备cyCombine格式的数据（修复版本）---
  cat("\n--- 准备数据（转换为cyCombine tibble格式）---\n")

  # **关键修复1：检查CD51数据是否可用**
  has_cd51 <- sum(!is.na(ALL_CELL_DATA$CD51_PE_intensity)) > 0
  cat(sprintf("CD51数据可用性：%s\n", ifelse(has_cd51, "是", "否")))

  # **关键修复2：添加原始行号作为ID（在过滤前）**
  ALL_CELL_DATA$original_row_id <- seq_len(nrow(ALL_CELL_DATA))

  # **关键修复3：根据CD51数据可用性决定使用的markers**
  if (has_cd51) {
    markers <- c("GFP_intensity", "CD51_PE_intensity")
    cat("使用的Markers:", paste(markers, collapse = ", "), "\n")

    # 清理数据：移除NA和无效值
    df <- ALL_CELL_DATA %>%
      filter(!is.na(GFP_intensity), !is.na(CD51_PE_intensity)) %>%
      filter(GFP_intensity > 0, CD51_PE_intensity > 0)
  } else {
    markers <- c("GFP_intensity")
    cat("使用的Markers:", paste(markers, collapse = ", "), "\n")
    cat("警告：仅使用GFP强度进行批次矫正\n")

    # 清理数据：仅检查GFP
    df <- ALL_CELL_DATA %>%
      filter(!is.na(GFP_intensity)) %>%
      filter(GFP_intensity > 0)
  }

  # **关键修复4：确保batch是字符串类型，不是factor**
  df <- df %>%
    mutate(
      id = original_row_id,  # 使用原始行号作为ID
      batch = as.character(batch),  # 字符串类型
      sample = as.character(sample),
      condition = as.character(group)
    ) %>%
    as_tibble() %>%
    select(id, sample, batch, condition, all_of(markers), everything())

  cat(sprintf("清理后数据：%d个细胞\n", nrow(df)))

  # **关键修复5：检查数据完整性**
  cat("\n数据完整性检查：\n")
  for (m in markers) {
    n_na <- sum(is.na(df[[m]]))
    n_inf <- sum(is.infinite(df[[m]]))
    n_neg <- sum(df[[m]] <= 0)
    cat(sprintf("  %s: NA=%d, Inf=%d, <=0=%d\n", m, n_na, n_inf, n_neg))
  }

  # 如果还有NA或Inf，移除这些行
  for (m in markers) {
    df <- df %>%
      filter(!is.na(.data[[m]]), !is.infinite(.data[[m]]), .data[[m]] > 0)
  }
  cat(sprintf("最终数据：%d个细胞\n", nrow(df)))

  # asinh变换
  cat("\n对数据进行asinh变换（cofactor=150）...\n")
  cofactor <- 150

  for (m in markers) {
    df[[m]] <- asinh(df[[m]] / cofactor)
  }

  cat("变换后数据范围：\n")
  for (m in markers) {
    cat(sprintf("  %s: [%.3f, %.3f]\n", m, min(df[[m]]), max(df[[m]])))
  }

  # 再次检查变换后的数据
  cat("\n变换后数据完整性检查：\n")
  for (m in markers) {
    n_na <- sum(is.na(df[[m]]))
    n_inf <- sum(is.infinite(df[[m]]))
    cat(sprintf("  %s: NA=%d, Inf=%d\n", m, n_na, n_inf))
  }

  # --- 7.3 批次效应检测 ---
  cat("\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  cat("=== 批次效应检测 ===\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")

  # 创建输出目录
  out_dir <- "batch_effect_analysis"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # 快速批次效应检测
  cat("\n--- 快速批次效应检测 ---\n")
  tryCatch({
    detect_batch_effect_express(
      df,
      downsample = min(10000, nrow(df)),
      out_dir = file.path(out_dir, "express"),
      markers = markers,
      batch_col = "batch"
    )
    cat("✓ 快速检测完成！\n")
    cat(sprintf("  结果保存在: %s/express/\n", out_dir))
  }, error = function(e) {
    cat(sprintf("✗ 快速检测失败: %s\n", e$message))
    cat("  继续执行后续步骤...\n")
  })

  # 完整批次效应检测
  cat("\n--- 完整批次效应检测 ---\n")
  tryCatch({
    detect_batch_effect(
      df,
      batch_col = "batch",
      out_dir = file.path(out_dir, "full"),
      seed = 473,
      name = "FACS_GFP_CD51",
      markers = markers
    )
    cat("✓ 完整检测完成！\n")
    cat(sprintf("  结果保存在: %s/full/\n", out_dir))
  }, error = function(e) {
    cat(sprintf("✗ 完整检测失败: %s\n", e$message))
    cat("  继续执行批次矫正...\n")
  })

  # 手动创建基本可视化
  cat("\n--- 手动创建批次效应可视化 ---\n")

  tryCatch({
    # 密度图
    p_gfp_density <- ggplot(df, aes(x = GFP_intensity, fill = batch)) +
      geom_density(alpha = 0.5) +
      labs(title = "GFP强度分布（按批次）",
           subtitle = "如果分布明显不同，说明存在批次效应",
           x = "GFP强度 (asinh变换)", y = "密度") +
      theme_bw() +
      theme(legend.position = "bottom")

    if (has_cd51) {
      p_cd51_density <- ggplot(df, aes(x = CD51_PE_intensity, fill = batch)) +
        geom_density(alpha = 0.5) +
        labs(title = "CD51强度分布（按批次）",
             x = "CD51强度 (asinh变换)", y = "密度") +
        theme_bw() +
        theme(legend.position = "bottom")
    }

    # 箱线图
    p_batch_boxplot <- ggplot(df, aes(x = batch, y = GFP_intensity, fill = batch)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.1, size = 0.3) +
      labs(title = "GFP强度（按批次）") +
      theme_bw() +
      theme(legend.position = "none")

    p_group_boxplot <- ggplot(df, aes(x = condition, y = GFP_intensity, fill = condition)) +
      geom_boxplot(alpha = 0.7) +
      labs(title = "GFP强度（按组）") +
      theme_bw() +
      theme(legend.position = "none")

    # 组合图
    if (has_cd51) {
      batch_effect_plots <- plot_grid(
        p_gfp_density, p_cd51_density,
        p_batch_boxplot, p_group_boxplot,
        ncol = 2, nrow = 2
      )
    } else {
      batch_effect_plots <- plot_grid(
        p_gfp_density, p_batch_boxplot, p_group_boxplot,
        ncol = 2, nrow = 2
      )
    }

    print(batch_effect_plots)
    ggsave(file.path(out_dir, "manual_batch_effect_check.png"),
           batch_effect_plots, width = 12, height = 10, dpi = 150)
    cat(sprintf("✓ 手动可视化已保存: %s/manual_batch_effect_check.png\n", out_dir))
  }, error = function(e) {
    cat(sprintf("✗ 手动可视化失败: %s\n", e$message))
  })

  # --- 7.4 批次矫正 ---
  cat("\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  cat("=== 批次矫正 ===\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")

  # 检查批次平衡性
  cat("\n批次细胞数：\n")
  print(table(df$batch))
  batch_counts <- as.numeric(table(df$batch))
  if (min(batch_counts) / max(batch_counts) < 0.1) {
    cat("\n⚠️ 警告：批次严重不平衡（最小批次 < 最大批次的10%）\n")
    cat("   这可能影响批次矫正的效果\n\n")
  }

  cat("\n运行batch_correct()...\n")

  tryCatch({
    corrected <- df %>%
      batch_correct(
        markers = markers,
        covar = "condition",
        xdim = 8,
        ydim = 8,
        norm_method = "scale",
        seed = 473
      )

    cat("\n✓ 批次矫正完成！\n")

    # --- 7.5 评估矫正效果（可选步骤）---
    cat("\n")
    cat(paste0(rep("=", 60), collapse = ""), "\n")
    cat("=== 评估矫正效果（可选步骤）===\n")
    cat(paste0(rep("=", 60), collapse = ""), "\n")

    # 将评估步骤放到单独的tryCatch中，失败不影响归一化值保存
    tryCatch({
      # 聚类分析
      cat("\n--- 聚类分析 ---\n")
      labels <- corrected %>%
        create_som(
          rlen = 10,
          xdim = 8,
          ydim = 8,
          markers = markers
        )

      corrected <- corrected %>%
        mutate(som = labels)

      df <- df %>%
        mutate(som = labels)

      cat("✓ 聚类完成\n")

      # EMD评估（需要emdist包）
      cat("\n--- EMD评估 ---\n")
      if (!requireNamespace("emdist", quietly = TRUE)) {
        cat("⚠️ emdist包未安装，跳过EMD评估\n")
        cat("   如需EMD评估，请运行: install.packages('emdist')\n")
      } else {
        emd_val <- df %>%
          evaluate_emd(
            corrected,
            binSize = 0.1,
            markers = markers,
            cell_col = "som"
          )

        cat(sprintf("EMD减少: %.2f\n", emd_val$reduction))
        cat("  (越接近1越好，表示批次效应消除得越彻底)\n")

        emd_plots <- plot_grid(emd_val$violin, emd_val$scatterplot)
        print(emd_plots)
        ggsave(file.path(out_dir, "EMD_evaluation.png"), emd_plots,
               width = 12, height = 5, dpi = 150)
        cat("✓ EMD评估完成\n")
      }

      # MAD评估
      cat("\n--- MAD评估 ---\n")
      mad_val <- df %>%
        evaluate_mad(
          corrected,
          filter_limit = NULL,
          markers = markers,
          cell_col = "som"
        )

      cat(sprintf("MAD分数: %.4f\n", mad_val))
      cat("  (越接近0越好，表示生物学信息保留得越完整)\n")
      cat("✓ MAD评估完成\n")

    }, error = function(e) {
      cat(sprintf("\n⚠️ 评估步骤失败: %s\n", e$message))
      cat("   这不影响归一化值的保存，继续执行...\n")
    })

    # 密度图对比（独立的tryCatch）
    tryCatch({
      cat("\n--- 矫正前后密度对比 ---\n")

      p_before <- ggplot(df, aes(x = GFP_intensity, fill = batch)) +
        geom_density(alpha = 0.5) +
        labs(title = "矫正前", x = "GFP强度") +
        theme_bw() +
        coord_cartesian(xlim = range(c(df$GFP_intensity, corrected$GFP_intensity)))

      p_after <- ggplot(corrected, aes(x = GFP_intensity, fill = batch)) +
        geom_density(alpha = 0.5) +
        labs(title = "矫正后", x = "GFP强度") +
        theme_bw() +
        coord_cartesian(xlim = range(c(df$GFP_intensity, corrected$GFP_intensity)))

      comparison_plot <- plot_grid(p_before, p_after, ncol = 2)
      print(comparison_plot)
      ggsave(file.path(out_dir, "GFP_before_after_comparison.png"),
             comparison_plot, width = 12, height = 5, dpi = 150)
      cat("✓ 密度对比图已保存\n")
    }, error = function(e) {
      cat(sprintf("⚠️ 密度图生成失败: %s\n", e$message))
    })

    # --- 7.6 保存矫正后的数据（修复版本）---
    cat("\n")
    cat(paste0(rep("=", 60), collapse = ""), "\n")
    cat("=== 保存矫正后的数据 ===\n")
    cat(paste0(rep("=", 60), collapse = ""), "\n")

    # **关键修复6：反向变换并重命名列**
    corrected_export <- corrected

    # 反向asinh变换
    corrected_export$GFP_intensity_corrected_raw <- sinh(corrected_export$GFP_intensity) * cofactor
    corrected_export$GFP_intensity_corrected_asinh <- corrected_export$GFP_intensity

    if (has_cd51) {
      corrected_export$CD51_PE_intensity_corrected_raw <- sinh(corrected_export$CD51_PE_intensity) * cofactor
      corrected_export$CD51_PE_intensity_corrected_asinh <- corrected_export$CD51_PE_intensity
    }

    # 保存数据
    saveRDS(corrected_export, file.path(out_dir, "corrected_data.rds"))
    write_csv(corrected_export, file.path(out_dir, "corrected_data.csv"))

    # **关键修复7：将矫正值添加回ALL_CELL_DATA（改进匹配逻辑）**
    cat("\n将矫正值添加回ALL_CELL_DATA...\n")

    ALL_CELL_DATA$GFP_intensity_normalized <- NA_real_
    if (has_cd51) {
      ALL_CELL_DATA$CD51_PE_intensity_normalized <- NA_real_
    }

    # 通过id列匹配（id对应original_row_id）
    cat(sprintf("  矫正数据行数：%d\n", nrow(corrected_export)))
    cat(sprintf("  原始数据行数：%d\n", nrow(ALL_CELL_DATA)))

    # 创建ID到归一化值的映射
    id_to_gfp <- setNames(corrected_export$GFP_intensity_corrected_raw, corrected_export$id)

    if (has_cd51) {
      id_to_cd51 <- setNames(corrected_export$CD51_PE_intensity_corrected_raw, corrected_export$id)
    }

    # 批量匹配
    matched_ids <- intersect(ALL_CELL_DATA$original_row_id, corrected_export$id)
    cat(sprintf("  成功匹配%d个细胞的ID\n", length(matched_ids)))

    # 通过向量化操作赋值
    for (cell_id in matched_ids) {
      row_idx <- which(ALL_CELL_DATA$original_row_id == cell_id)
      ALL_CELL_DATA$GFP_intensity_normalized[row_idx] <- id_to_gfp[[as.character(cell_id)]]

      if (has_cd51) {
        ALL_CELL_DATA$CD51_PE_intensity_normalized[row_idx] <- id_to_cd51[[as.character(cell_id)]]
      }
    }

    # 验证匹配结果
    n_normalized <- sum(!is.na(ALL_CELL_DATA$GFP_intensity_normalized))
    cat(sprintf("  最终有%d个细胞获得了归一化值\n", n_normalized))
    cat(sprintf("  归一化率：%.1f%%\n", 100 * n_normalized / nrow(ALL_CELL_DATA)))

    # 保存更新后的ALL_CELL_DATA
    saveRDS(ALL_CELL_DATA, file.path(out_dir, "ALL_CELL_DATA_with_normalization.rds"))
    write_csv(ALL_CELL_DATA, file.path(out_dir, "ALL_CELL_DATA_with_normalization.csv"))
    cat(sprintf("✓ 更新后的数据保存在: %s/\n", out_dir))

    cat("\n✓ 批次矫正和归一化完成！\n")

  }, error = function(e) {
    cat(sprintf("\n✗ 批次矫正失败: %s\n", e$message))
    cat("\n错误诊断信息：\n")
    cat(sprintf("  - 数据行数：%d\n", nrow(df)))
    cat(sprintf("  - 批次数：%d\n", length(unique(df$batch))))
    cat(sprintf("  - Markers数：%d\n", length(markers)))
    cat("\n将继续分析，但不使用归一化值\n")
    DO_BATCH_CORRECTION <- FALSE
  })

} else if (nrow(ALL_CELL_DATA) == 0) {
  cat("\n警告：没有细胞数据，跳过批次归一化\n")
  DO_BATCH_CORRECTION <- FALSE
} else {
  cat("\n批次归一化已关闭\n")
}

# --- 8. 可视化 ---
# ==================

if (nrow(ALL_CELL_DATA) > 0) {
  cat("\n=== 生成可视化图表 ===\n")

  # 准备绘图数据
  plot_data <- ALL_CELL_DATA
  if (nrow(plot_data) > 100000) {
    cat(sprintf("  数据量较大（%d个细胞），随机采样100,000个细胞用于可视化\n", nrow(plot_data)))
    plot_data <- plot_data[sample(nrow(plot_data), 100000), ]
  }

  plot_data <- plot_data %>%
    filter(GFP_intensity > 0)

  if (!all(is.na(plot_data$CD51_PE_intensity))) {
    plot_data <- plot_data %>% filter(CD51_PE_intensity > 0)
  }

  cat(sprintf("  绘图使用%d个细胞\n", nrow(plot_data)))

  # 图1：GFP强度分布
  cat("  生成图1：GFP荧光强度分布（原始值）...\n")
  p1 <- ggplot(plot_data, aes(x = GFP_intensity, fill = group)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~batch, ncol = 1) +
    scale_x_log10() +
    labs(title = "CD51+ GFP+ 细胞的GFP荧光强度分布（按批次分面）",
         subtitle = "原始值（批次矫正前）",
         x = "GFP荧光强度 (log10)",
         y = "密度") +
    theme_bw() +
    theme(legend.position = "bottom")

  print(p1)

  # 图2：各样本GFP强度分布
  cat("  生成图2：各样本GFP强度分布...\n")
  p2 <- ggplot(plot_data, aes(x = sample, y = GFP_intensity, fill = group)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, alpha = 0.5, outlier.size = 0.5) +
    facet_wrap(~batch, scales = "free_x", ncol = 1) +
    scale_y_log10() +
    labs(title = "各样本GFP荧光强度分布",
         x = "样本",
         y = "GFP荧光强度 (log10)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")

  print(p2)

  # 图3：批次效应
  if (length(unique(plot_data$batch)) > 1) {
    cat("  生成图3：批次效应（归一化前）...\n")
    p3 <- ggplot(plot_data, aes(x = batch, y = GFP_intensity, fill = batch)) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.2, alpha = 0.5, outlier.size = 0.5) +
      scale_y_log10() +
      labs(title = "批次效应（归一化前）",
           x = "批次",
           y = "GFP荧光强度 (log10)") +
      theme_bw() +
      theme(legend.position = "none")

    print(p3)
  }

  # 归一化后的图
  if (DO_BATCH_CORRECTION && "GFP_intensity_normalized" %in% colnames(plot_data)) {
    plot_data_norm <- plot_data %>%
      filter(!is.na(GFP_intensity_normalized), GFP_intensity_normalized > 0)

    if (nrow(plot_data_norm) > 0) {
      cat("  生成图4：GFP荧光强度分布（归一化后）...\n")
      p4 <- ggplot(plot_data_norm, aes(x = GFP_intensity_normalized, fill = group)) +
        geom_density(alpha = 0.5) +
        facet_wrap(~batch, ncol = 1) +
        scale_x_log10() +
        labs(title = "CD51+ GFP+ 细胞的GFP荧光强度分布（按批次分面）",
             subtitle = "归一化后（批次矫正后）",
             x = "GFP荧光强度（归一化）(log10)",
             y = "密度") +
        theme_bw() +
        theme(legend.position = "bottom")

      print(p4)

      if (length(unique(plot_data_norm$batch)) > 1) {
        cat("  生成图5：批次效应（归一化后）...\n")
        p5 <- ggplot(plot_data_norm, aes(x = batch, y = GFP_intensity_normalized, fill = batch)) +
          geom_violin(alpha = 0.7) +
          geom_boxplot(width = 0.2, alpha = 0.5, outlier.size = 0.5) +
          scale_y_log10() +
          labs(title = "批次效应（归一化后）",
               x = "批次",
               y = "GFP荧光强度（归一化）(log10)") +
          theme_bw() +
          theme(legend.position = "none")

        print(p5)
      }
    }
  }

  # GFP vs CD51散点图
  if (!all(is.na(plot_data$CD51_PE_intensity))) {
    cat("  生成图6：GFP vs CD51散点图...\n")
    p6 <- ggplot(plot_data, aes(x = CD51_PE_intensity, y = GFP_intensity, color = group)) +
      geom_point(alpha = 0.2, size = 0.5) +
      facet_wrap(~batch) +
      scale_x_log10() +
      scale_y_log10() +
      labs(title = "GFP vs CD51 荧光强度",
           x = "CD51 (PE)荧光强度 (log10)",
           y = "GFP荧光强度 (log10)") +
      theme_bw() +
      theme(legend.position = "bottom")

    print(p6)
  }

  cat("  ✓ 可视化完成\n")
}

# --- 9. 导出结果 ---
# ====================

if (EXPORT_DATA && nrow(ALL_CELL_DATA) > 0) {
  cat("\n=== 导出结果 ===\n")

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # 确定要导出的列
  cols_to_export <- c("group", "batch", "sample", "cell_id",
                      "GFP_intensity", "CD51_PE_intensity")

  if (DO_BATCH_CORRECTION && "GFP_intensity_normalized" %in% colnames(ALL_CELL_DATA)) {
    cols_to_export <- c(cols_to_export, "GFP_intensity_normalized")
    if ("CD51_PE_intensity_normalized" %in% colnames(ALL_CELL_DATA)) {
      cols_to_export <- c(cols_to_export, "CD51_PE_intensity_normalized")
    }
  }

  # CSV导出
  if (EXPORT_CSV) {
    cat("\n[ 1/3 ] 导出CSV文件...\n")

    sample_list <- split(ALL_CELL_DATA, ALL_CELL_DATA$sample)

    for (sample_name in names(sample_list)) {
      sample_data <- sample_list[[sample_name]]

      safe_name <- gsub("[^A-Za-z0-9_-]", "_", sample_name)
      csv_file <- file.path(OUTPUT_DIR, paste0(safe_name, "_SingleCell_GFP.csv"))

      export_data <- sample_data[, intersect(cols_to_export, colnames(sample_data))]

      write.csv(export_data, csv_file, row.names = FALSE)

      cat(sprintf("  ✓ %s: %d个细胞 → %s\n",
                  sample_name, nrow(sample_data), basename(csv_file)))
    }

    cat(sprintf("\nCSV文件已保存到：%s/\n", OUTPUT_DIR))
  }

  # Excel导出
  if (EXPORT_EXCEL) {
    cat("\n[ 2/3 ] 导出Excel文件...\n")
    excel_file <- file.path(OUTPUT_DIR, paste0("FACS_SingleCell_GFP_", timestamp, ".xlsx"))

    wb <- createWorkbook()

    # 工作表1：样本汇总
    addWorksheet(wb, "Sample_Summary")
    writeDataTable(wb, "Sample_Summary", SAMPLE_SUMMARY)
    cat("  ✓ 工作表1: Sample_Summary\n")

    # 工作表2：完整单细胞数据
    max_excel_rows <- 500000

    if (nrow(ALL_CELL_DATA) <= max_excel_rows) {
      addWorksheet(wb, "All_Cells_Data")
      export_data <- ALL_CELL_DATA[, intersect(cols_to_export, colnames(ALL_CELL_DATA))]
      writeDataTable(wb, "All_Cells_Data", export_data)
      cat(sprintf("  ✓ 工作表2: All_Cells_Data（%d个细胞）\n", nrow(ALL_CELL_DATA)))
    } else {
      cat(sprintf("  ⚠ 单细胞数据过大（%d行），仅在CSV中提供\n", nrow(ALL_CELL_DATA)))
    }

    # 工作表3：README
    addWorksheet(wb, "README")
    readme_rows <- list(
      c("group", "实验组"),
      c("batch", "批次标识"),
      c("sample", "样本名称"),
      c("cell_id", "细胞编号"),
      c("GFP_intensity", "GFP荧光强度原始值"),
      c("CD51_PE_intensity", "CD51/PE荧光强度原始值")
    )

    if (DO_BATCH_CORRECTION && "GFP_intensity_normalized" %in% colnames(ALL_CELL_DATA)) {
      readme_rows <- c(readme_rows, list(
        c("GFP_intensity_normalized", "GFP荧光强度归一化值（批次矫正后）"),
        c("CD51_PE_intensity_normalized", "CD51/PE荧光强度归一化值（批次矫正后）")
      ))
    }

    readme_text <- as.data.frame(do.call(rbind, readme_rows))
    colnames(readme_text) <- c("Column", "Description")

    writeDataTable(wb, "README", readme_text)
    cat("  ✓ 工作表3: README\n")

    saveWorkbook(wb, excel_file, overwrite = TRUE)
    cat(sprintf("\nExcel文件已保存：%s\n", basename(excel_file)))
  }

  # R工作空间
  cat("\n[ 3/3 ] 保存R工作空间...\n")
  rdata_file <- file.path(OUTPUT_DIR, paste0("FACS_Workspace_", timestamp, ".RData"))

  save(ALL_CELL_DATA,
       SAMPLE_SUMMARY,
       ALL_GATING_SETS,
       SAMPLE_INFO,
       file = rdata_file)

  cat(sprintf("  ✓ R工作空间已保存：%s\n", basename(rdata_file)))
}

# --- 10. 最终报告 ---
# =====================

cat("\n")
cat(paste0(rep("=", 80), collapse = ""), "\n")
cat("=== 分析完成 ===\n")
cat(paste0(rep("=", 80), collapse = ""), "\n\n")

cat(sprintf("分析时间：%s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("总样本数：%d\n", nrow(SAMPLE_SUMMARY)))
cat(sprintf("总细胞数：%d\n", nrow(ALL_CELL_DATA)))

if (nrow(SAMPLE_SUMMARY) > 0) {
  cat("\n各样本细胞数：\n")
  print(SAMPLE_SUMMARY)
}

cat("\n导出的数据列：\n")
cat("  原始数据：\n")
cat("    - GFP_intensity: 每个细胞的GFP荧光强度（原始值）\n")
cat("    - CD51_PE_intensity: 每个细胞的CD51/PE荧光强度（原始值）\n")

if (DO_BATCH_CORRECTION && "GFP_intensity_normalized" %in% colnames(ALL_CELL_DATA)) {
  cat("  归一化数据（批次矫正后）：\n")
  cat("    - GFP_intensity_normalized: GFP荧光强度（归一化值）\n")
  if ("CD51_PE_intensity_normalized" %in% colnames(ALL_CELL_DATA)) {
    cat("    - CD51_PE_intensity_normalized: CD51/PE荧光强度（归一化值）\n")
  }
  cat("\n✓ 批次矫正成功！建议使用归一化值进行组间比较。\n")
} else {
  cat("\n⚠ 批次矫正未执行或失败，仅提供原始值。\n")
}

cat("\n输出文件位置：", OUTPUT_DIR, "/\n")

cat(paste0(rep("=", 80), collapse = ""), "\n\n")

# ==============================================================================
# 脚本结束
# ==============================================================================
