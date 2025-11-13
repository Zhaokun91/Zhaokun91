# FACS流式细胞数据分析工具

## 项目简介

本工具用于分析FACS（流式细胞术）数据，特别针对具有批次效应的实验设计。支持：

- 读取 `.fcs` 原始数据文件
- 读取 FlowJo `.wsp` 工作空间文件
- 自动补偿矫正
- 基于WT对照的门控定义
- 批次效应矫正
- 提取GFP+CD51+细胞的荧光强度

## 实验设计

本工具专为以下实验设计开发：

### 样本组
- **BM组**: 新鲜骨髓
- **Culture组**: 培养后骨髓
- **Sabq组**: 异体骨骨髓

### 批次效应
- **批次A**: BM组 + Culture组（同一天运行）
- **批次B**: Sabq组（另一天运行）

### 对照
- **门控对照**: 每个批次都有WT（野生型）老鼠的不染骨髓作为阴性对照
- **补偿对照**: 每个批次都有当天的单染管

### 特殊情况
- 基因编辑老鼠样本，目标细胞自带荧光

## 目录结构

```
facs_analysis/
├── data/                      # 数据目录
│   ├── fcs_files/            # FCS原始文件
│   └── analysis.wsp          # FlowJo工作空间文件
├── output/                    # 输出目录
├── scripts/                   # 脚本目录
│   └── facs_analyzer.py      # 主分析脚本
├── utils/                     # 工具函数
├── config_template.json       # 配置文件模板
└── requirements.txt           # Python依赖包
```

## 安装

### 1. 克隆仓库

```bash
git clone https://github.com/Zhaokun91/Zhaokun91.git
cd Zhaokun91/facs_analysis
```

### 2. 创建虚拟环境（推荐）

```bash
python -m venv venv
source venv/bin/activate  # Linux/Mac
# 或
venv\Scripts\activate     # Windows
```

### 3. 安装依赖

```bash
pip install -r requirements.txt
```

## 配置

### 1. 准备配置文件

复制配置文件模板：

```bash
cp config_template.json config.json
```

### 2. 编辑配置文件

编辑 `config.json`，填入您的实验信息：

```json
{
  "data_paths": {
    "fcs_directory": "./data/fcs_files",
    "wsp_file": "./data/analysis.wsp",
    "output_directory": "./output"
  },

  "channels": {
    "gfp": "GFP-A",
    "cd51": "CD51-APC-A"
  },

  "batch_mapping": {
    "batch_A": {
      "samples": [
        "BM_1", "BM_2", "BM_3",
        "Culture_1", "Culture_2", "Culture_3"
      ],
      "wt_control": "WT_unstained_A",
      "compensation_controls": [
        "comp_GFP_A", "comp_CD51_A"
      ]
    },
    "batch_B": {
      "samples": [
        "Sabq_1", "Sabq_2", "Sabq_3"
      ],
      "wt_control": "WT_unstained_B",
      "compensation_controls": [
        "comp_GFP_B", "comp_CD51_B"
      ]
    }
  }
}
```

**重要提示**：

1. **样本ID**: `samples` 列表中的样本ID必须与FCS文件名（不含.fcs扩展名）完全匹配
2. **通道名称**: `gfp` 和 `cd51` 的通道名称必须与FCS文件中的通道名称完全一致
3. **对照样本**: 确保 `wt_control` 和 `compensation_controls` 中的样本ID也在FCS文件中

### 3. 准备数据

将您的数据文件放入相应目录：

```
data/
├── fcs_files/
│   ├── BM_1.fcs
│   ├── BM_2.fcs
│   ├── Culture_1.fcs
│   ├── Sabq_1.fcs
│   ├── WT_unstained_A.fcs
│   ├── WT_unstained_B.fcs
│   ├── comp_GFP_A.fcs
│   ├── comp_CD51_A.fcs
│   └── ...
└── analysis.wsp               # FlowJo工作空间文件（可选）
```

## 使用方法

### 方法1：使用配置文件（推荐）

```bash
cd scripts
python facs_analyzer.py --config ../config.json
```

### 方法2：命令行参数

```bash
cd scripts
python facs_analyzer.py \
    --fcs-dir ../data/fcs_files \
    --wsp ../data/analysis.wsp \
    --output-dir ../output \
    --gfp-channel GFP-A \
    --cd51-channel CD51-APC-A
```

### 方法3：在Python中使用

```python
from facs_analyzer import FACSAnalyzer

# 创建分析器
analyzer = FACSAnalyzer(config_path='config.json')

# 加载FCS文件
analyzer.load_fcs_files('./data/fcs_files')

# 加载FlowJo工作空间（可选）
analyzer.load_flowjo_workspace('./data/analysis.wsp', './data/fcs_files')

# 组织批次信息
batch_mapping = {
    'batch_A': {
        'samples': ['BM_1', 'BM_2', 'Culture_1'],
        'wt_control': 'WT_unstained_A',
        'compensation_controls': ['comp_GFP_A', 'comp_CD51_A']
    },
    'batch_B': {
        'samples': ['Sabq_1', 'Sabq_2'],
        'wt_control': 'WT_unstained_B',
        'compensation_controls': ['comp_GFP_B', 'comp_CD51_B']
    }
}
analyzer.organize_samples_by_batch(batch_mapping)

# 执行分析
results = analyzer.analyze_all_samples(
    gfp_channel='GFP-A',
    cd51_channel='CD51-APC-A',
    output_dir='./output'
)

# 生成汇总统计
analyzer.generate_summary_statistics(results, output_dir='./output')
```

## 输出结果

分析完成后，会在输出目录生成以下文件：

### 1. 单个样本结果

每个样本的GFP+CD51+细胞数据：

```
output/
├── BM_1_GFP+CD51+.csv
├── BM_2_GFP+CD51+.csv
├── Culture_1_GFP+CD51+.csv
├── Sabq_1_GFP+CD51+.csv
└── ...
```

每个CSV文件包含：
- `GFP-A`: GFP荧光强度（批次矫正后）
- `CD51-APC-A`: CD51荧光强度（批次矫正后）
- `sample_id`: 样本ID
- `batch`: 批次名称

### 2. 合并结果

所有样本的合并数据（批次矫正后）：

```
output/all_samples_GFP+CD51+_corrected.csv
```

### 3. 汇总统计

每个样本的统计信息：

```
output/summary_statistics.csv
```

包含每个样本的：
- 细胞数量
- 荧光强度均值、中位数、标准差
- 最小值、最大值

## 分析流程

本工具执行以下分析步骤：

1. **加载数据**
   - 读取FCS文件
   - 读取FlowJo工作空间（如有）

2. **补偿矫正**
   - 自动应用FCS文件中的补偿矩阵
   - 或使用补偿对照计算补偿矩阵

3. **门控定义**
   - 使用WT不染对照定义阴性/阳性阈值
   - 默认使用99%分位数作为阈值
   - 筛选GFP+CD51+双阳性细胞

4. **批次矫正**
   - 使用分位数归一化方法
   - 将批次B的数据映射到批次A的分布
   - 支持其他方法：均值矫正、中位数矫正

5. **数据提取**
   - 提取每个细胞的荧光强度
   - 保存单个样本和合并数据
   - 生成汇总统计

## 高级功能

### 自定义门控阈值

如果您想手动指定阈值而不是从WT对照计算：

```python
# 在analyze_all_samples之前
analyzer.batch_info['batch_A']['gfp_threshold'] = 1500
analyzer.batch_info['batch_A']['cd51_threshold'] = 2000
```

### 选择批次矫正方法

支持三种批次矫正方法：

1. **quantile**（分位数归一化，默认）- 推荐用于大多数情况
2. **mean**（均值矫正）- 简单快速
3. **median**（中位数矫正）- 对异常值更稳健

在 `perform_batch_correction` 方法中指定：

```python
corrected = analyzer.perform_batch_correction(
    batch_a_data,
    batch_b_data,
    method='quantile'  # 或 'mean', 'median'
)
```

### 使用FlowJo门控

如果您在FlowJo中已经定义了门控，可以直接使用：

```python
gate_indices = analyzer.get_gate_indices(
    sample,
    gate_name='GFP+CD51+',  # FlowJo中的门控名称
    # 或使用门控路径
    gate_path=['Lymphocytes', 'GFP+', 'CD51+']
)
```

## 常见问题

### Q1: 如何确定通道名称？

使用以下代码查看FCS文件中的所有通道：

```python
from flowkit import Sample
sample = Sample('path/to/file.fcs')
print(sample.pnn_labels)  # 打印所有通道名称
```

### Q2: 如果没有FlowJo工作空间文件怎么办？

没问题！本工具可以独立工作，不需要FlowJo工作空间文件。只需使用手动门控即可。

### Q3: 如何调整门控阈值的严格程度？

修改 `calculate_threshold_from_wt` 方法中的 `percentile` 参数：

- 99（默认）：较宽松，包含更多细胞
- 99.5：中等
- 99.9：较严格，只包含最强的阳性细胞

### Q4: 批次矫正是否必需？

如果您的所有样本都在同一天运行，可以跳过批次矫正。但对于多批次实验，强烈推荐进行批次矫正以消除技术变异。

### Q5: 如何处理多个荧光标记？

修改 `extract_gated_data` 调用，指定额外的通道：

```python
gated_data = analyzer.extract_gated_data(
    sample,
    gate_indices,
    channels=['GFP-A', 'CD51-APC-A', 'CD45-PE-A', 'CD11b-PerCP-A']
)
```

## 引用和贡献

如果您在研究中使用了本工具，请引用本仓库。

欢迎提交问题报告和改进建议！

## 许可证

本项目采用 MIT 许可证。

## 联系方式

如有问题，请通过GitHub Issues联系。

---

**作者**: Zhaokun91
**最后更新**: 2025-11-13
