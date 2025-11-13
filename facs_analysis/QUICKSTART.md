# 快速开始指南

## 5分钟上手FACS数据分析

### 步骤1: 安装依赖

```bash
cd facs_analysis
pip install -r requirements.txt
```

### 步骤2: 准备数据

将您的FCS文件放入 `data/fcs_files/` 目录：

```
data/fcs_files/
├── BM_1.fcs
├── BM_2.fcs
├── Culture_1.fcs
├── Sabq_1.fcs
├── WT_unstained_A.fcs
└── WT_unstained_B.fcs
```

### 步骤3: 配置参数

复制并编辑配置文件：

```bash
cp config_template.json config.json
```

修改 `config.json` 中的样本ID和通道名称以匹配您的数据。

### 步骤4: 运行分析

```bash
cd scripts
python facs_analyzer.py --config ../config.json
```

### 步骤5: 查看结果

结果将保存在 `output/` 目录：

- `all_samples_GFP+CD51+_corrected.csv` - 所有样本的合并数据
- `summary_statistics.csv` - 汇总统计
- `[样本ID]_GFP+CD51+.csv` - 单个样本数据

## 常见任务

### 查看FCS文件的通道名称

```python
from flowkit import Sample

sample = Sample('data/fcs_files/BM_1.fcs')
print(sample.pnn_labels)
```

### 只分析特定样本

修改 `config.json` 中的 `batch_mapping` -> `samples` 列表。

### 调整门控阈值

在配置文件中设置：

```json
"analysis_parameters": {
  "threshold_percentile": 99.5
}
```

### 更改批次矫正方法

在脚本中修改 `perform_batch_correction` 的 `method` 参数：
- `"quantile"` - 分位数归一化（推荐）
- `"mean"` - 均值矫正
- `"median"` - 中位数矫正

## 获取帮助

运行以下命令查看所有选项：

```bash
python facs_analyzer.py --help
```

查看完整文档: [README.md](README.md)
