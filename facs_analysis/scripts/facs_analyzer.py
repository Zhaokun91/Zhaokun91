#!/usr/bin/env python3
"""
FACS数据分析脚本
用于处理FlowJo工作空间和FCS文件，提取GFP+CD51+门控下的细胞荧光强度

作者: Zhaokun91
日期: 2025-11-13
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    import flowkit as fk
    from flowkit import Sample, Session
except ImportError:
    print("警告: flowkit未安装，请运行: pip install flowkit")
    sys.exit(1)

try:
    import flowutils
except ImportError:
    print("警告: flowutils未安装，请运行: pip install flowutils")
    sys.exit(1)


class FACSAnalyzer:
    """FACS数据分析器"""

    def __init__(self, config_path=None):
        """
        初始化分析器

        参数:
            config_path: 配置文件路径（JSON格式）
        """
        self.config = self._load_config(config_path) if config_path else {}
        self.session = None
        self.samples = {}
        self.batch_info = {
            'batch_A': {'samples': [], 'controls': {}},
            'batch_B': {'samples': [], 'controls': {}}
        }

    def _load_config(self, config_path):
        """加载配置文件"""
        with open(config_path, 'r', encoding='utf-8') as f:
            return json.load(f)

    def load_flowjo_workspace(self, wsp_path, fcs_dir):
        """
        加载FlowJo工作空间文件

        参数:
            wsp_path: .wsp文件路径
            fcs_dir: FCS文件目录
        """
        print(f"正在加载FlowJo工作空间: {wsp_path}")
        try:
            self.session = Session(
                fcs_samples=fcs_dir,
                wsp_file=wsp_path
            )
            print(f"成功加载 {len(self.session.sample_ids)} 个样本")
            return True
        except Exception as e:
            print(f"加载工作空间时出错: {str(e)}")
            return False

    def load_fcs_files(self, fcs_dir, file_pattern='*.fcs'):
        """
        加载FCS文件

        参数:
            fcs_dir: FCS文件目录
            file_pattern: 文件匹配模式
        """
        print(f"正在从 {fcs_dir} 加载FCS文件...")
        fcs_path = Path(fcs_dir)
        fcs_files = list(fcs_path.glob(file_pattern))

        print(f"找到 {len(fcs_files)} 个FCS文件")

        for fcs_file in fcs_files:
            try:
                sample_id = fcs_file.stem
                sample = Sample(str(fcs_file))
                self.samples[sample_id] = sample
                print(f"  ✓ 加载: {sample_id}")
            except Exception as e:
                print(f"  ✗ 加载失败 {fcs_file.name}: {str(e)}")

        return len(self.samples)

    def organize_samples_by_batch(self, batch_mapping):
        """
        按批次组织样本

        参数:
            batch_mapping: 批次映射字典
            格式: {
                'batch_A': {
                    'samples': ['BM_1', 'BM_2', 'Culture_1', ...],
                    'wt_control': 'WT_unstained_A',
                    'compensation_controls': ['comp_GFP_A', 'comp_CD51_A', ...]
                },
                'batch_B': {...}
            }
        """
        self.batch_info = batch_mapping
        print("\n批次组织:")
        for batch_name, batch_data in batch_mapping.items():
            print(f"\n{batch_name}:")
            print(f"  样本数: {len(batch_data.get('samples', []))}")
            print(f"  WT对照: {batch_data.get('wt_control', 'N/A')}")
            print(f"  补偿对照数: {len(batch_data.get('compensation_controls', []))}")

    def apply_compensation(self, sample, compensation_matrix=None):
        """
        应用补偿矫正

        参数:
            sample: Sample对象
            compensation_matrix: 补偿矩阵（如果为None，使用FCS文件中的矩阵）

        返回:
            补偿后的Sample对象
        """
        try:
            if compensation_matrix is not None:
                sample.apply_compensation(compensation_matrix)
            else:
                # 使用FCS文件中嵌入的补偿矩阵
                if sample.has_compensation:
                    sample.apply_compensation()
                else:
                    print(f"警告: 样本无补偿矩阵")
            return sample
        except Exception as e:
            print(f"应用补偿时出错: {str(e)}")
            return sample

    def get_gate_indices(self, sample, gate_name=None, gate_path=None):
        """
        获取门控索引

        参数:
            sample: Sample对象
            gate_name: 门控名称（在FlowJo中定义）
            gate_path: 门控路径（用于层级门控）

        返回:
            布尔数组，指示哪些事件通过门控
        """
        if self.session is None:
            print("警告: 未加载FlowJo工作空间，无法应用门控")
            return None

        try:
            # 从FlowJo工作空间获取门控
            if gate_name:
                gate_indices = self.session.get_gate_indices(
                    sample_id=sample.id,
                    gate_name=gate_name
                )
                return gate_indices
            elif gate_path:
                gate_indices = self.session.get_gate_indices(
                    sample_id=sample.id,
                    gate_path=gate_path
                )
                return gate_indices
        except Exception as e:
            print(f"获取门控索引时出错: {str(e)}")
            return None

    def manual_gating(self, sample, channel_x, channel_y,
                      threshold_x=None, threshold_y=None,
                      positive_x=True, positive_y=True):
        """
        手动门控（用于定义阳性细胞）

        参数:
            sample: Sample对象
            channel_x: X轴通道名称（如 'GFP-A'）
            channel_y: Y轴通道名称（如 'CD51-APC-A'）
            threshold_x: X轴阈值
            threshold_y: Y轴阈值
            positive_x: True表示 > threshold_x，False表示 < threshold_x
            positive_y: True表示 > threshold_y，False表示 < threshold_y

        返回:
            布尔数组，指示哪些事件通过门控
        """
        try:
            # 获取事件数据
            events = sample.get_events()

            # 获取通道索引
            channel_idx_x = sample.get_channel_index(channel_x)
            channel_idx_y = sample.get_channel_index(channel_y)

            # 获取数据
            data_x = events[:, channel_idx_x]
            data_y = events[:, channel_idx_y]

            # 应用门控
            if positive_x:
                gate_x = data_x > threshold_x
            else:
                gate_x = data_x < threshold_x

            if positive_y:
                gate_y = data_y > threshold_y
            else:
                gate_y = data_y < threshold_y

            # 组合门控
            gate_indices = gate_x & gate_y

            return gate_indices

        except Exception as e:
            print(f"手动门控时出错: {str(e)}")
            return None

    def calculate_threshold_from_wt(self, wt_sample, channel, percentile=99):
        """
        从WT不染样本计算阈值

        参数:
            wt_sample: WT样本对象
            channel: 通道名称
            percentile: 百分位数（默认99%）

        返回:
            阈值
        """
        try:
            events = wt_sample.get_events()
            channel_idx = wt_sample.get_channel_index(channel)
            data = events[:, channel_idx]
            threshold = np.percentile(data, percentile)
            print(f"  通道 {channel} 的 {percentile}% 阈值: {threshold:.2f}")
            return threshold
        except Exception as e:
            print(f"计算阈值时出错: {str(e)}")
            return None

    def extract_gated_data(self, sample, gate_indices, channels=None):
        """
        提取门控后的数据

        参数:
            sample: Sample对象
            gate_indices: 门控索引（布尔数组）
            channels: 要提取的通道列表（如果为None，提取所有通道）

        返回:
            DataFrame，包含门控后的数据
        """
        try:
            events = sample.get_events()

            # 应用门控
            if gate_indices is not None:
                gated_events = events[gate_indices]
            else:
                gated_events = events

            # 获取通道名称
            if channels is None:
                channels = sample.pnn_labels

            # 创建DataFrame
            df = pd.DataFrame(gated_events, columns=sample.pnn_labels)

            # 只保留指定的通道
            if channels:
                available_channels = [ch for ch in channels if ch in df.columns]
                df = df[available_channels]

            return df

        except Exception as e:
            print(f"提取门控数据时出错: {str(e)}")
            return None

    def perform_batch_correction(self, batch_a_data, batch_b_data,
                                  method='quantile'):
        """
        执行批次矫正

        参数:
            batch_a_data: 批次A的数据（DataFrame）
            batch_b_data: 批次B的数据（DataFrame）
            method: 矫正方法（'quantile', 'mean', 'median'）

        返回:
            矫正后的数据字典 {'batch_A': df_a, 'batch_B': df_b}
        """
        print(f"\n正在执行批次矫正（方法: {method}）...")

        corrected_data = {}

        if method == 'quantile':
            # 分位数归一化
            # 将批次B的数据映射到批次A的分布
            for col in batch_a_data.columns:
                if col in batch_b_data.columns:
                    # 计算批次A的分位数
                    quantiles_a = np.linspace(0, 1, 100)
                    values_a = np.quantile(batch_a_data[col], quantiles_a)

                    # 将批次B映射到批次A的分布
                    batch_b_data[col] = np.interp(
                        batch_b_data[col].rank(pct=True),
                        quantiles_a,
                        values_a
                    )

            corrected_data['batch_A'] = batch_a_data
            corrected_data['batch_B'] = batch_b_data

        elif method == 'mean':
            # 均值矫正
            for col in batch_a_data.columns:
                if col in batch_b_data.columns:
                    mean_a = batch_a_data[col].mean()
                    mean_b = batch_b_data[col].mean()
                    batch_b_data[col] = batch_b_data[col] * (mean_a / mean_b)

            corrected_data['batch_A'] = batch_a_data
            corrected_data['batch_B'] = batch_b_data

        elif method == 'median':
            # 中位数矫正
            for col in batch_a_data.columns:
                if col in batch_b_data.columns:
                    median_a = batch_a_data[col].median()
                    median_b = batch_b_data[col].median()
                    batch_b_data[col] = batch_b_data[col] * (median_a / median_b)

            corrected_data['batch_A'] = batch_a_data
            corrected_data['batch_B'] = batch_b_data

        print("  批次矫正完成")
        return corrected_data

    def analyze_all_samples(self, gfp_channel='GFP-A', cd51_channel='CD51-APC-A',
                           output_dir='./output'):
        """
        分析所有样本并提取GFP+CD51+细胞的荧光强度

        参数:
            gfp_channel: GFP通道名称
            cd51_channel: CD51通道名称
            output_dir: 输出目录

        返回:
            结果字典
        """
        print("\n" + "="*60)
        print("开始分析所有样本")
        print("="*60)

        results = {
            'batch_A': {},
            'batch_B': {}
        }

        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)

        # 处理每个批次
        for batch_name, batch_data in self.batch_info.items():
            print(f"\n处理 {batch_name}...")

            # 1. 从WT对照计算阈值
            wt_control_id = batch_data.get('wt_control')
            if wt_control_id and wt_control_id in self.samples:
                wt_sample = self.samples[wt_control_id]
                print(f"\n使用WT对照: {wt_control_id}")

                threshold_gfp = self.calculate_threshold_from_wt(
                    wt_sample, gfp_channel, percentile=99
                )
                threshold_cd51 = self.calculate_threshold_from_wt(
                    wt_sample, cd51_channel, percentile=99
                )
            else:
                print(f"警告: 未找到WT对照，使用默认阈值")
                threshold_gfp = 1000
                threshold_cd51 = 1000

            # 2. 处理该批次的所有样本
            batch_results = []

            for sample_id in batch_data.get('samples', []):
                if sample_id not in self.samples:
                    print(f"  警告: 样本 {sample_id} 未找到")
                    continue

                print(f"\n  分析样本: {sample_id}")
                sample = self.samples[sample_id]

                # 应用补偿
                sample = self.apply_compensation(sample)

                # 应用门控 (GFP+ AND CD51+)
                gate_indices = self.manual_gating(
                    sample,
                    channel_x=gfp_channel,
                    channel_y=cd51_channel,
                    threshold_x=threshold_gfp,
                    threshold_y=threshold_cd51,
                    positive_x=True,
                    positive_y=True
                )

                if gate_indices is not None:
                    n_gated = np.sum(gate_indices)
                    total = len(gate_indices)
                    percentage = (n_gated / total) * 100
                    print(f"    门控细胞数: {n_gated} / {total} ({percentage:.2f}%)")

                    # 提取门控后的数据
                    gated_data = self.extract_gated_data(
                        sample,
                        gate_indices,
                        channels=[gfp_channel, cd51_channel]
                    )

                    if gated_data is not None and len(gated_data) > 0:
                        # 添加样本信息
                        gated_data['sample_id'] = sample_id
                        gated_data['batch'] = batch_name
                        batch_results.append(gated_data)

                        # 保存单个样本的结果
                        sample_output_path = os.path.join(
                            output_dir,
                            f"{sample_id}_GFP+CD51+.csv"
                        )
                        gated_data.to_csv(sample_output_path, index=False)
                        print(f"    保存到: {sample_output_path}")

            # 合并该批次的所有结果
            if batch_results:
                batch_combined = pd.concat(batch_results, ignore_index=True)
                results[batch_name] = batch_combined
                print(f"\n  {batch_name} 合并数据: {len(batch_combined)} 个细胞")

        # 3. 执行批次矫正
        if results['batch_A'] is not None and results['batch_B'] is not None:
            if len(results['batch_A']) > 0 and len(results['batch_B']) > 0:
                print("\n执行批次矫正...")

                # 只对荧光强度列进行矫正
                batch_a_fi = results['batch_A'][[gfp_channel, cd51_channel]]
                batch_b_fi = results['batch_B'][[gfp_channel, cd51_channel]]

                corrected = self.perform_batch_correction(
                    batch_a_fi,
                    batch_b_fi,
                    method='quantile'
                )

                # 更新结果
                results['batch_A'][[gfp_channel, cd51_channel]] = corrected['batch_A']
                results['batch_B'][[gfp_channel, cd51_channel]] = corrected['batch_B']

        # 4. 保存最终结果
        final_combined = pd.concat([
            results['batch_A'] if len(results['batch_A']) > 0 else pd.DataFrame(),
            results['batch_B'] if len(results['batch_B']) > 0 else pd.DataFrame()
        ], ignore_index=True)

        final_output_path = os.path.join(output_dir, 'all_samples_GFP+CD51+_corrected.csv')
        final_combined.to_csv(final_output_path, index=False)

        print("\n" + "="*60)
        print(f"分析完成！")
        print(f"总细胞数: {len(final_combined)}")
        print(f"最终结果保存到: {final_output_path}")
        print("="*60)

        return results

    def generate_summary_statistics(self, results, output_dir='./output'):
        """
        生成汇总统计

        参数:
            results: 分析结果字典
            output_dir: 输出目录
        """
        print("\n生成汇总统计...")

        summary_data = []

        for batch_name, batch_df in results.items():
            if len(batch_df) == 0:
                continue

            # 按样本分组统计
            for sample_id in batch_df['sample_id'].unique():
                sample_df = batch_df[batch_df['sample_id'] == sample_id]

                # 提取荧光通道
                fi_channels = [col for col in sample_df.columns
                              if col not in ['sample_id', 'batch']]

                summary_row = {
                    'sample_id': sample_id,
                    'batch': batch_name,
                    'n_cells': len(sample_df)
                }

                # 计算每个通道的统计量
                for channel in fi_channels:
                    summary_row[f'{channel}_mean'] = sample_df[channel].mean()
                    summary_row[f'{channel}_median'] = sample_df[channel].median()
                    summary_row[f'{channel}_std'] = sample_df[channel].std()
                    summary_row[f'{channel}_min'] = sample_df[channel].min()
                    summary_row[f'{channel}_max'] = sample_df[channel].max()

                summary_data.append(summary_row)

        # 创建汇总DataFrame
        summary_df = pd.DataFrame(summary_data)

        # 保存汇总统计
        summary_path = os.path.join(output_dir, 'summary_statistics.csv')
        summary_df.to_csv(summary_path, index=False)

        print(f"汇总统计保存到: {summary_path}")
        print("\n样本汇总:")
        print(summary_df[['sample_id', 'batch', 'n_cells']].to_string(index=False))

        return summary_df


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(
        description='FACS数据分析工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 使用配置文件
  python facs_analyzer.py --config config.json

  # 手动指定参数
  python facs_analyzer.py --fcs-dir ./data --wsp ./data/analysis.wsp
        """
    )

    parser.add_argument('--config', type=str, help='配置文件路径（JSON格式）')
    parser.add_argument('--fcs-dir', type=str, help='FCS文件目录')
    parser.add_argument('--wsp', type=str, help='FlowJo工作空间文件（.wsp）')
    parser.add_argument('--output-dir', type=str, default='./output',
                       help='输出目录（默认: ./output）')
    parser.add_argument('--gfp-channel', type=str, default='GFP-A',
                       help='GFP通道名称（默认: GFP-A）')
    parser.add_argument('--cd51-channel', type=str, default='CD51-APC-A',
                       help='CD51通道名称（默认: CD51-APC-A）')

    args = parser.parse_args()

    # 创建分析器
    analyzer = FACSAnalyzer(config_path=args.config)

    # 加载数据
    if args.fcs_dir:
        analyzer.load_fcs_files(args.fcs_dir)

    if args.wsp and args.fcs_dir:
        analyzer.load_flowjo_workspace(args.wsp, args.fcs_dir)

    # 如果使用配置文件，加载批次信息
    if args.config:
        with open(args.config, 'r', encoding='utf-8') as f:
            config = json.load(f)
            if 'batch_mapping' in config:
                analyzer.organize_samples_by_batch(config['batch_mapping'])

    # 执行分析
    if len(analyzer.samples) > 0:
        results = analyzer.analyze_all_samples(
            gfp_channel=args.gfp_channel,
            cd51_channel=args.cd51_channel,
            output_dir=args.output_dir
        )

        # 生成汇总统计
        analyzer.generate_summary_statistics(results, output_dir=args.output_dir)
    else:
        print("错误: 未加载任何样本")


if __name__ == '__main__':
    main()
