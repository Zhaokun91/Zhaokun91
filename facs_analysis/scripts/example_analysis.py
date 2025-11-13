#!/usr/bin/env python3
"""
FACS数据分析示例脚本
展示如何使用FACSAnalyzer类进行数据分析
"""

from facs_analyzer import FACSAnalyzer

def main():
    """示例分析流程"""

    # 1. 创建分析器实例
    print("="*60)
    print("FACS数据分析示例")
    print("="*60)

    analyzer = FACSAnalyzer()

    # 2. 加载FCS文件
    print("\n步骤1: 加载FCS文件")
    fcs_dir = '../data/fcs_files'  # 根据实际路径修改
    n_samples = analyzer.load_fcs_files(fcs_dir)
    print(f"成功加载 {n_samples} 个样本")

    # 3. （可选）加载FlowJo工作空间
    # 如果您有.wsp文件，可以使用以下代码：
    # wsp_file = '../data/analysis.wsp'
    # analyzer.load_flowjo_workspace(wsp_file, fcs_dir)

    # 4. 组织批次信息
    print("\n步骤2: 组织批次信息")
    batch_mapping = {
        'batch_A': {
            'description': '批次A：BM和Culture组',
            'samples': [
                # 根据您的实际样本名称修改
                'BM_1', 'BM_2', 'BM_3',
                'Culture_1', 'Culture_2', 'Culture_3'
            ],
            'wt_control': 'WT_unstained_A',
            'compensation_controls': [
                'comp_GFP_A',
                'comp_CD51_A'
            ]
        },
        'batch_B': {
            'description': '批次B：Sabq组',
            'samples': [
                # 根据您的实际样本名称修改
                'Sabq_1', 'Sabq_2', 'Sabq_3'
            ],
            'wt_control': 'WT_unstained_B',
            'compensation_controls': [
                'comp_GFP_B',
                'comp_CD51_B'
            ]
        }
    }

    analyzer.organize_samples_by_batch(batch_mapping)

    # 5. 执行分析
    print("\n步骤3: 执行分析")
    results = analyzer.analyze_all_samples(
        gfp_channel='GFP-A',      # 根据实际通道名称修改
        cd51_channel='CD51-APC-A',  # 根据实际通道名称修改
        output_dir='../output'
    )

    # 6. 生成汇总统计
    print("\n步骤4: 生成汇总统计")
    summary = analyzer.generate_summary_statistics(
        results,
        output_dir='../output'
    )

    print("\n分析完成！")
    print("="*60)


if __name__ == '__main__':
    main()
