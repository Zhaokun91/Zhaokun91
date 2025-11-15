# 批次矫正关键修复 - 快速参考

## 7个关键修复点

### ✅ 修复1：CD51数据智能处理

**问题**: 过滤时要求CD51不为NA，但有些样本没有CD51数据

**原代码**:
```r
df <- ALL_CELL_DATA %>%
  filter(!is.na(GFP_intensity), !is.na(CD51_PE_intensity)) %>%  # ❌ 如果CD51全是NA，数据被清空
  filter(GFP_intensity > 0, CD51_PE_intensity > 0)
```

**修复后**:
```r
# 先检查CD51是否可用
has_cd51 <- sum(!is.na(ALL_CELL_DATA$CD51_PE_intensity)) > 0

if (has_cd51) {
  markers <- c("GFP_intensity", "CD51_PE_intensity")
  df <- ALL_CELL_DATA %>%
    filter(!is.na(GFP_intensity), !is.na(CD51_PE_intensity)) %>%
    filter(GFP_intensity > 0, CD51_PE_intensity > 0)
} else {
  markers <- c("GFP_intensity")  # ✅ 只用GFP
  df <- ALL_CELL_DATA %>%
    filter(!is.na(GFP_intensity)) %>%
    filter(GFP_intensity > 0)
}
```

---

### ✅ 修复2：batch列类型

**问题**: cyCombine需要字符串类型的batch，不是factor

**原代码**:
```r
df <- df %>%
  mutate(
    batch = as.factor(batch)  # ❌ factor类型
  )
```

**修复后**:
```r
df <- df %>%
  mutate(
    batch = as.character(batch)  # ✅ 字符串类型
  )
```

---

### ✅ 修复3：添加原始行号ID（在过滤前）

**问题**: 过滤后无法匹配回原始数据

**原代码**:
```r
df <- ALL_CELL_DATA %>%
  filter(...) %>%
  mutate(id = row_number())  # ❌ 过滤后的行号，无法匹配原始数据
```

**修复后**:
```r
# ✅ 在过滤前添加原始行号
ALL_CELL_DATA$original_row_id <- seq_len(nrow(ALL_CELL_DATA))

df <- ALL_CELL_DATA %>%
  filter(...) %>%
  mutate(id = original_row_id)  # ✅ 使用原始行号
```

---

### ✅ 修复4：数据完整性验证

**问题**: NA/Inf值导致cyCombine函数失败

**新增代码**:
```r
# asinh变换前检查
cat("\n数据完整性检查：\n")
for (m in markers) {
  n_na <- sum(is.na(df[[m]]))
  n_inf <- sum(is.infinite(df[[m]]))
  n_neg <- sum(df[[m]] <= 0)
  cat(sprintf("  %s: NA=%d, Inf=%d, <=0=%d\n", m, n_na, n_inf, n_neg))
}

# 移除问题行
for (m in markers) {
  df <- df %>%
    filter(!is.na(.data[[m]]), !is.infinite(.data[[m]]), .data[[m]] > 0)
}

# 变换后再次检查
cat("\n变换后数据完整性检查：\n")
for (m in markers) {
  n_na <- sum(is.na(df[[m]]))
  n_inf <- sum(is.infinite(df[[m]]))
  cat(sprintf("  %s: NA=%d, Inf=%d\n", m, n_na, n_inf))
}
```

---

### ✅ 修复5：改进归一化值匹配

**问题**: 矫正值无法正确添加回ALL_CELL_DATA

**原代码**:
```r
# 逐行匹配（慢且容易出错）
for (i in 1:nrow(corrected_export)) {
  matched <- which(ALL_CELL_DATA$sample == corrected_export$sample[i] &
                   ALL_CELL_DATA$cell_id == corrected_export$cell_id[i])
  # 可能匹配失败
}
```

**修复后**:
```r
# ✅ 通过ID创建映射表（快速且准确）
id_to_gfp <- setNames(corrected_export$GFP_intensity_corrected_raw,
                      corrected_export$id)

# 批量匹配
matched_ids <- intersect(ALL_CELL_DATA$original_row_id, corrected_export$id)

for (cell_id in matched_ids) {
  row_idx <- which(ALL_CELL_DATA$original_row_id == cell_id)
  ALL_CELL_DATA$GFP_intensity_normalized[row_idx] <- id_to_gfp[[as.character(cell_id)]]
}

# ✅ 验证匹配率
n_normalized <- sum(!is.na(ALL_CELL_DATA$GFP_intensity_normalized))
cat(sprintf("  归一化率：%.1f%%\n", 100 * n_normalized / nrow(ALL_CELL_DATA)))
```

---

### ✅ 修复6：CD51归一化值处理

**问题**: CD51归一化值没有正确生成

**新增代码**:
```r
# 反向asinh变换
corrected_export$GFP_intensity_corrected_raw <- sinh(corrected_export$GFP_intensity) * cofactor

# ✅ 只在有CD51数据时才处理CD51归一化
if (has_cd51) {
  corrected_export$CD51_PE_intensity_corrected_raw <- sinh(corrected_export$CD51_PE_intensity) * cofactor

  # 添加回ALL_CELL_DATA
  ALL_CELL_DATA$CD51_PE_intensity_normalized <- NA_real_
  id_to_cd51 <- setNames(corrected_export$CD51_PE_intensity_corrected_raw, corrected_export$id)

  for (cell_id in matched_ids) {
    row_idx <- which(ALL_CELL_DATA$original_row_id == cell_id)
    ALL_CELL_DATA$CD51_PE_intensity_normalized[row_idx] <- id_to_cd51[[as.character(cell_id)]]
  }
}
```

---

### ✅ 修复7：增强错误处理

**问题**: 错误信息不清晰，难以调试

**新增代码**:
```r
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

}, error = function(e) {
  # ✅ 详细的错误诊断
  cat(sprintf("\n✗ 批次矫正失败: %s\n", e$message))
  cat("\n错误诊断信息：\n")
  cat(sprintf("  - 数据行数：%d\n", nrow(df)))
  cat(sprintf("  - 批次数：%d\n", length(unique(df$batch))))
  cat(sprintf("  - Markers数：%d\n", length(markers)))
  cat("\n建议：\n")
  cat("  - 检查批次平衡性\n")
  cat("  - 尝试减小SOM网格大小（xdim=4, ydim=4）\n")
  cat("  - 或者不进行批次矫正，在下游分析中将批次作为协变量控制\n")
  DO_BATCH_CORRECTION <- FALSE
})
```

---

## 快速测试检查清单

运行修复后的脚本，应该看到：

- [ ] ✅ 没有 "NA/NaN argument" 警告
- [ ] ✅ 没有 "invalid 'type' (character)" 警告
- [ ] ✅ 输出显示 "✓ 批次矫正完成！"
- [ ] ✅ 输出显示归一化率接近100%
- [ ] ✅ 生成了 `batch_effect_analysis/` 目录和图表
- [ ] ✅ CSV文件包含 `GFP_intensity_normalized` 列
- [ ] ✅ Excel文件包含归一化列
- [ ] ✅ `GFP_intensity_normalized` 列不全是NA

## 验证命令

```R
# 加载结果
data <- read.csv("batch_effect_analysis/ALL_CELL_DATA_with_normalization.csv")

# 检查1：归一化列是否存在
"GFP_intensity_normalized" %in% colnames(data)
# 期望: TRUE

# 检查2：归一化值分布
summary(data$GFP_intensity_normalized)
# 期望: 不全是NA，有Min/Max/Mean等统计值

# 检查3：归一化率
sum(!is.na(data$GFP_intensity_normalized)) / nrow(data) * 100
# 期望: 接近100%

# 检查4：原始值vs归一化值对比
library(ggplot2)
ggplot(data, aes(x = GFP_intensity, y = GFP_intensity_normalized)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "原始值 vs 归一化值")
# 期望: 看到点的分布，不是空白图
```

## 如果还是失败

1. 检查R包版本是否最新
2. 检查数据完整性（是否有足够的细胞数）
3. 查看详细的错误诊断信息
4. 参考 `BATCH_CORRECTION_FIX_NOTES.md` 进行深度排查
