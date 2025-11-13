#!/usr/bin/env python3
"""
辅助工具函数
用于FACS数据分析的各种实用函数
"""

import numpy as np
import pandas as pd
from pathlib import Path


def list_fcs_files(directory):
    """
    列出目录中的所有FCS文件

    参数:
        directory: FCS文件目录

    返回:
        FCS文件路径列表
    """
    path = Path(directory)
    fcs_files = list(path.glob('*.fcs'))
    return sorted([str(f) for f in fcs_files])


def get_sample_info_from_filename(filename):
    """
    从文件名提取样本信息

    参数:
        filename: FCS文件名

    返回:
        字典，包含样本信息
    """
    # 假设文件名格式为: SampleType_Number.fcs
    # 例如: BM_1.fcs, Culture_2.fcs

    stem = Path(filename).stem
    parts = stem.split('_')

    info = {
        'sample_id': stem,
        'sample_type': parts[0] if len(parts) > 0 else 'unknown',
        'replicate': parts[1] if len(parts) > 1 else '1'
    }

    return info


def check_channel_exists(sample, channel_name):
    """
    检查样本中是否存在指定通道

    参数:
        sample: Sample对象
        channel_name: 通道名称

    返回:
        布尔值
    """
    try:
        sample.get_channel_index(channel_name)
        return True
    except:
        return False


def print_channel_info(sample):
    """
    打印样本的通道信息

    参数:
        sample: Sample对象
    """
    print(f"\n样本: {sample.id}")
    print(f"事件数: {sample.event_count}")
    print(f"\n通道列表:")
    print("-" * 60)

    for i, (pnn, pns) in enumerate(zip(sample.pnn_labels, sample.pns_labels)):
        print(f"{i+1:3d}. {pnn:20s} | {pns}")


def calculate_cv(data):
    """
    计算变异系数 (Coefficient of Variation)

    参数:
        data: 数据数组

    返回:
        CV值（百分比）
    """
    mean = np.mean(data)
    std = np.std(data)
    cv = (std / mean) * 100 if mean != 0 else 0
    return cv


def calculate_signal_to_noise(signal, background):
    """
    计算信噪比

    参数:
        signal: 信号强度（如阳性细胞的荧光强度）
        background: 背景强度（如阴性细胞的荧光强度）

    返回:
        信噪比
    """
    signal_mean = np.mean(signal)
    background_mean = np.mean(background)
    snr = signal_mean / background_mean if background_mean != 0 else 0
    return snr


def export_to_csv(data, output_path, include_index=False):
    """
    导出数据到CSV文件

    参数:
        data: DataFrame
        output_path: 输出文件路径
        include_index: 是否包含索引
    """
    data.to_csv(output_path, index=include_index)
    print(f"数据已导出到: {output_path}")


def combine_csv_files(file_list, output_path):
    """
    合并多个CSV文件

    参数:
        file_list: CSV文件路径列表
        output_path: 输出文件路径

    返回:
        合并后的DataFrame
    """
    dfs = []
    for file_path in file_list:
        df = pd.read_csv(file_path)
        dfs.append(df)

    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv(output_path, index=False)

    print(f"合并了 {len(file_list)} 个文件")
    print(f"总行数: {len(combined)}")
    print(f"输出到: {output_path}")

    return combined


def calculate_percentile_threshold(data, percentile=99):
    """
    计算百分位数阈值

    参数:
        data: 数据数组
        percentile: 百分位数

    返回:
        阈值
    """
    threshold = np.percentile(data, percentile)
    return threshold


def apply_log_transform(data, base=10):
    """
    应用对数变换

    参数:
        data: 数据数组
        base: 对数底数

    返回:
        变换后的数据
    """
    # 避免log(0)
    data_shifted = data + 1
    transformed = np.log(data_shifted) / np.log(base)
    return transformed


def apply_logicle_transform(data, T=262144, M=4.5, W=0.5, A=0):
    """
    应用Logicle变换（用于流式细胞数据）

    参数:
        data: 数据数组
        T: 顶部值
        M: 对数十进制数
        W: 线性范围宽度
        A: 额外负值范围

    返回:
        变换后的数据

    注意: 这是简化版本，完整的Logicle变换需要专门的库
    """
    # 这里使用简化的双曲正弦反函数变换
    # 实际应用中建议使用flowkit或FlowUtils库的Logicle实现
    r = np.sinh(data * np.log(10) * M / (2 * T)) / np.sinh(np.log(10) * M / (2 * T))
    return r


def detect_outliers_iqr(data, multiplier=1.5):
    """
    使用IQR方法检测异常值

    参数:
        data: 数据数组
        multiplier: IQR乘数（默认1.5）

    返回:
        异常值的布尔掩码
    """
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1

    lower_bound = q1 - multiplier * iqr
    upper_bound = q3 + multiplier * iqr

    outliers = (data < lower_bound) | (data > upper_bound)
    return outliers


def normalize_to_control(data, control_data, method='median'):
    """
    相对于对照组归一化数据

    参数:
        data: 待归一化的数据
        control_data: 对照组数据
        method: 归一化方法（'median', 'mean', 'max'）

    返回:
        归一化后的数据
    """
    if method == 'median':
        control_value = np.median(control_data)
    elif method == 'mean':
        control_value = np.mean(control_data)
    elif method == 'max':
        control_value = np.max(control_data)
    else:
        raise ValueError(f"未知的归一化方法: {method}")

    if control_value == 0:
        print("警告: 对照值为0，无法归一化")
        return data

    normalized = data / control_value
    return normalized


def calculate_fold_change(treatment, control):
    """
    计算倍数变化

    参数:
        treatment: 处理组数据
        control: 对照组数据

    返回:
        倍数变化值
    """
    treatment_mean = np.mean(treatment)
    control_mean = np.mean(control)

    if control_mean == 0:
        return np.inf

    fold_change = treatment_mean / control_mean
    return fold_change


def print_summary_stats(data, name="数据"):
    """
    打印数据的汇总统计

    参数:
        data: 数据数组或DataFrame列
        name: 数据名称
    """
    print(f"\n{name} 统计:")
    print(f"  数量:    {len(data)}")
    print(f"  均值:    {np.mean(data):.2f}")
    print(f"  中位数:  {np.median(data):.2f}")
    print(f"  标准差:  {np.std(data):.2f}")
    print(f"  最小值:  {np.min(data):.2f}")
    print(f"  最大值:  {np.max(data):.2f}")
    print(f"  CV:      {calculate_cv(data):.2f}%")


# 批次矫正相关函数

def quantile_normalize(data_list):
    """
    分位数归一化（用于多个样本）

    参数:
        data_list: 数据列表（每个元素是一个样本的数据）

    返回:
        归一化后的数据列表
    """
    # 将数据转换为矩阵
    n_samples = len(data_list)
    max_len = max([len(d) for d in data_list])

    # 创建矩阵
    matrix = np.zeros((max_len, n_samples))
    for i, data in enumerate(data_list):
        matrix[:len(data), i] = data

    # 排序每列
    sorted_matrix = np.sort(matrix, axis=0)

    # 计算行均值
    row_means = np.mean(sorted_matrix, axis=1)

    # 重新排列
    normalized_list = []
    for i in range(n_samples):
        ranks = data_list[i].argsort().argsort()
        normalized = row_means[ranks]
        normalized_list.append(normalized)

    return normalized_list


if __name__ == '__main__':
    # 测试函数
    print("FACS分析辅助工具函数")
    print("可用函数:")
    print("  - list_fcs_files: 列出FCS文件")
    print("  - get_sample_info_from_filename: 提取样本信息")
    print("  - print_channel_info: 打印通道信息")
    print("  - calculate_cv: 计算变异系数")
    print("  - calculate_signal_to_noise: 计算信噪比")
    print("  - 以及更多...")
