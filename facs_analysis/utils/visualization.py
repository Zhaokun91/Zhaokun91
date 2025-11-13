#!/usr/bin/env python3
"""
FACS数据可视化工具
用于生成流式细胞数据的各种图表
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# 设置绘图样式
sns.set_style("whitegrid")
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def plot_scatter_2d(data_x, data_y, xlabel='GFP-A', ylabel='CD51-APC-A',
                    title='2D Scatter Plot', threshold_x=None, threshold_y=None,
                    output_path=None, figsize=(8, 8)):
    """
    绘制2D散点图（类似FlowJo的dot plot）

    参数:
        data_x: X轴数据
        data_y: Y轴数据
        xlabel: X轴标签
        ylabel: Y轴标签
        title: 图标题
        threshold_x: X轴阈值线
        threshold_y: Y轴阈值线
        output_path: 输出文件路径
        figsize: 图形大小
    """
    fig, ax = plt.subplots(figsize=figsize)

    # 绘制散点图（使用半透明点）
    ax.scatter(data_x, data_y, alpha=0.3, s=1, c='blue')

    # 添加阈值线
    if threshold_x is not None:
        ax.axvline(threshold_x, color='red', linestyle='--', linewidth=2,
                   label=f'{xlabel} threshold')

    if threshold_y is not None:
        ax.axhline(threshold_y, color='red', linestyle='--', linewidth=2,
                   label=f'{ylabel} threshold')

    # 计算每个象限的细胞数
    if threshold_x is not None and threshold_y is not None:
        q1 = np.sum((data_x < threshold_x) & (data_y < threshold_y))  # 左下
        q2 = np.sum((data_x >= threshold_x) & (data_y < threshold_y))  # 右下
        q3 = np.sum((data_x < threshold_x) & (data_y >= threshold_y))  # 左上
        q4 = np.sum((data_x >= threshold_x) & (data_y >= threshold_y))  # 右上（双阳性）

        total = len(data_x)

        # 在图上标注百分比
        ax.text(0.02, 0.02, f'{q1}\n{q1/total*100:.1f}%',
                transform=ax.transAxes, fontsize=10)
        ax.text(0.98, 0.02, f'{q2}\n{q2/total*100:.1f}%',
                transform=ax.transAxes, fontsize=10, ha='right')
        ax.text(0.02, 0.98, f'{q3}\n{q3/total*100:.1f}%',
                transform=ax.transAxes, fontsize=10, va='top')
        ax.text(0.98, 0.98, f'{q4}\n{q4/total*100:.1f}%',
                transform=ax.transAxes, fontsize=10, ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend()

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"图表保存到: {output_path}")

    plt.show()
    return fig


def plot_histogram(data, channel_name='GFP-A', threshold=None,
                   title='Histogram', bins=256, output_path=None,
                   figsize=(10, 6)):
    """
    绘制单参数直方图

    参数:
        data: 数据数组
        channel_name: 通道名称
        threshold: 阈值线
        title: 图标题
        bins: 直方图bin数
        output_path: 输出文件路径
        figsize: 图形大小
    """
    fig, ax = plt.subplots(figsize=figsize)

    # 绘制直方图
    n, bins_edges, patches = ax.hist(data, bins=bins, alpha=0.7,
                                      color='steelblue', edgecolor='black')

    # 添加阈值线
    if threshold is not None:
        ax.axvline(threshold, color='red', linestyle='--', linewidth=2,
                   label=f'Threshold: {threshold:.0f}')

        # 计算阳性百分比
        positive_pct = np.sum(data > threshold) / len(data) * 100
        ax.text(threshold, max(n) * 0.9,
                f'Positive: {positive_pct:.1f}%',
                fontsize=10, color='red')

    ax.set_xlabel(channel_name, fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend()

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"图表保存到: {output_path}")

    plt.show()
    return fig


def plot_batch_comparison(batch_a_data, batch_b_data, channel_name='GFP-A',
                          output_path=None, figsize=(12, 6)):
    """
    比较两个批次的数据分布

    参数:
        batch_a_data: 批次A数据
        batch_b_data: 批次B数据
        channel_name: 通道名称
        output_path: 输出文件路径
        figsize: 图形大小
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    # 直方图比较
    ax1.hist(batch_a_data, bins=100, alpha=0.5, label='Batch A', color='blue')
    ax1.hist(batch_b_data, bins=100, alpha=0.5, label='Batch B', color='red')
    ax1.set_xlabel(channel_name, fontsize=12)
    ax1.set_ylabel('Count', fontsize=12)
    ax1.set_title('Distribution Comparison', fontsize=14)
    ax1.legend()

    # 箱线图比较
    data_combined = pd.DataFrame({
        'Batch A': batch_a_data,
        'Batch B': batch_b_data
    })
    data_combined.boxplot(ax=ax2)
    ax2.set_ylabel(channel_name, fontsize=12)
    ax2.set_title('Boxplot Comparison', fontsize=14)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"图表保存到: {output_path}")

    plt.show()
    return fig


def plot_sample_comparison(data_dict, channel_name='GFP-A',
                          plot_type='violin', output_path=None,
                          figsize=(12, 8)):
    """
    比较多个样本的数据分布

    参数:
        data_dict: 字典，键为样本ID，值为数据数组
        channel_name: 通道名称
        plot_type: 图表类型（'violin', 'box', 'strip'）
        output_path: 输出文件路径
        figsize: 图形大小
    """
    # 准备数据
    df_list = []
    for sample_id, data in data_dict.items():
        df = pd.DataFrame({
            'value': data,
            'sample': sample_id
        })
        df_list.append(df)

    df_combined = pd.concat(df_list, ignore_index=True)

    # 绘图
    fig, ax = plt.subplots(figsize=figsize)

    if plot_type == 'violin':
        sns.violinplot(data=df_combined, x='sample', y='value', ax=ax)
    elif plot_type == 'box':
        sns.boxplot(data=df_combined, x='sample', y='value', ax=ax)
    elif plot_type == 'strip':
        sns.stripplot(data=df_combined, x='sample', y='value', ax=ax, alpha=0.3)
    else:
        raise ValueError(f"未知的图表类型: {plot_type}")

    ax.set_xlabel('Sample', fontsize=12)
    ax.set_ylabel(channel_name, fontsize=12)
    ax.set_title(f'{channel_name} Comparison Across Samples', fontsize=14)
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"图表保存到: {output_path}")

    plt.show()
    return fig


def plot_correlation_heatmap(df, output_path=None, figsize=(10, 8)):
    """
    绘制相关性热图

    参数:
        df: DataFrame，包含多个通道的数据
        output_path: 输出文件路径
        figsize: 图形大小
    """
    # 只保留数值列
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    df_numeric = df[numeric_cols]

    # 计算相关性矩阵
    corr_matrix = df_numeric.corr()

    # 绘制热图
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0,
                square=True, linewidths=1, ax=ax, fmt='.2f')

    ax.set_title('Correlation Heatmap', fontsize=14, fontweight='bold')

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"图表保存到: {output_path}")

    plt.show()
    return fig


def plot_batch_correction_effect(before_a, before_b, after_a, after_b,
                                 channel_name='GFP-A', output_path=None,
                                 figsize=(14, 10)):
    """
    展示批次矫正的效果

    参数:
        before_a: 矫正前批次A数据
        before_b: 矫正前批次B数据
        after_a: 矫正后批次A数据
        after_b: 矫正后批次B数据
        channel_name: 通道名称
        output_path: 输出文件路径
        figsize: 图形大小
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)

    # 矫正前直方图
    axes[0, 0].hist(before_a, bins=100, alpha=0.5, label='Batch A', color='blue')
    axes[0, 0].hist(before_b, bins=100, alpha=0.5, label='Batch B', color='red')
    axes[0, 0].set_title('Before Correction - Histogram', fontsize=12)
    axes[0, 0].set_xlabel(channel_name)
    axes[0, 0].legend()

    # 矫正后直方图
    axes[0, 1].hist(after_a, bins=100, alpha=0.5, label='Batch A', color='blue')
    axes[0, 1].hist(after_b, bins=100, alpha=0.5, label='Batch B', color='red')
    axes[0, 1].set_title('After Correction - Histogram', fontsize=12)
    axes[0, 1].set_xlabel(channel_name)
    axes[0, 1].legend()

    # 矫正前箱线图
    df_before = pd.DataFrame({
        'Batch A': before_a,
        'Batch B': before_b
    })
    df_before.boxplot(ax=axes[1, 0])
    axes[1, 0].set_title('Before Correction - Boxplot', fontsize=12)
    axes[1, 0].set_ylabel(channel_name)

    # 矫正后箱线图
    df_after = pd.DataFrame({
        'Batch A': after_a,
        'Batch B': after_b
    })
    df_after.boxplot(ax=axes[1, 1])
    axes[1, 1].set_title('After Correction - Boxplot', fontsize=12)
    axes[1, 1].set_ylabel(channel_name)

    plt.suptitle(f'Batch Correction Effect - {channel_name}',
                 fontsize=16, fontweight='bold', y=1.00)
    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"图表保存到: {output_path}")

    plt.show()
    return fig


def plot_qc_report(summary_df, output_path=None, figsize=(14, 10)):
    """
    生成质量控制报告图表

    参数:
        summary_df: 汇总统计DataFrame
        output_path: 输出文件路径
        figsize: 图形大小
    """
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

    # 1. 细胞数量条形图
    ax1 = fig.add_subplot(gs[0, :])
    ax1.bar(summary_df['sample_id'], summary_df['n_cells'])
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Cell Count')
    ax1.set_title('Cell Count per Sample', fontweight='bold')
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # 2. GFP均值比较
    if 'GFP-A_mean' in summary_df.columns:
        ax2 = fig.add_subplot(gs[1, 0])
        colors = ['blue' if b == 'batch_A' else 'red'
                 for b in summary_df['batch']]
        ax2.bar(summary_df['sample_id'], summary_df['GFP-A_mean'], color=colors)
        ax2.set_xlabel('Sample')
        ax2.set_ylabel('GFP Mean Intensity')
        ax2.set_title('GFP Mean Intensity', fontweight='bold')
        plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # 3. CD51均值比较
    if 'CD51-APC-A_mean' in summary_df.columns:
        ax3 = fig.add_subplot(gs[1, 1])
        colors = ['blue' if b == 'batch_A' else 'red'
                 for b in summary_df['batch']]
        ax3.bar(summary_df['sample_id'], summary_df['CD51-APC-A_mean'],
                color=colors)
        ax3.set_xlabel('Sample')
        ax3.set_ylabel('CD51 Mean Intensity')
        ax3.set_title('CD51 Mean Intensity', fontweight='bold')
        plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # 4. CV比较
    if 'GFP-A_std' in summary_df.columns and 'GFP-A_mean' in summary_df.columns:
        ax4 = fig.add_subplot(gs[2, 0])
        cv = (summary_df['GFP-A_std'] / summary_df['GFP-A_mean']) * 100
        ax4.bar(summary_df['sample_id'], cv)
        ax4.set_xlabel('Sample')
        ax4.set_ylabel('CV (%)')
        ax4.set_title('GFP Coefficient of Variation', fontweight='bold')
        plt.setp(ax4.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # 5. 批次比较
    ax5 = fig.add_subplot(gs[2, 1])
    batch_counts = summary_df.groupby('batch')['n_cells'].sum()
    ax5.pie(batch_counts, labels=batch_counts.index, autopct='%1.1f%%',
            startangle=90)
    ax5.set_title('Cells per Batch', fontweight='bold')

    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"QC报告保存到: {output_path}")

    plt.show()
    return fig


if __name__ == '__main__':
    print("FACS数据可视化工具")
    print("可用函数:")
    print("  - plot_scatter_2d: 2D散点图")
    print("  - plot_histogram: 直方图")
    print("  - plot_batch_comparison: 批次比较")
    print("  - plot_sample_comparison: 样本比较")
    print("  - plot_correlation_heatmap: 相关性热图")
    print("  - plot_batch_correction_effect: 批次矫正效果")
    print("  - plot_qc_report: 质量控制报告")
