# output_optimizer.py - 优化输出显示的工具模块

import sys
import warnings

# 全局静默设置
def setup_silent_mode():
    """设置静默模式，抑制冗余输出"""
    
    # 静默RDKit警告
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
    except:
        pass
    
    # 静默Python警告
    warnings.filterwarnings('ignore', category=DeprecationWarning)
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=UserWarning)
    
    # 抑制OpenMP警告
    import os
    os.environ['OMP_NUM_THREADS'] = '1'
    
    # 抑制sklearn警告
    try:
        from sklearn.exceptions import ConvergenceWarning
        warnings.filterwarnings('ignore', category=ConvergenceWarning)
    except:
        pass
    
    # 重定向stderr中的RDKit输出
    class StderrFilter:
        """过滤stderr输出的包装器"""
        def __init__(self):
            self.original_stderr = sys.stderr
            self.buffer = []
            
        def write(self, text):
            # 过滤掉特定的警告信息
            if any(keyword in text for keyword in [
                'Explicit valence',
                'DEPRECATION',
                'RuntimeWarning',
                'threadpoolctl',
                'OpenMP',
                'libiomp',
                'libomp',
                'KMeans is known to have a memory leak',
                'Omitted undefined stereo'
            ]):
                return
            self.original_stderr.write(text)
            
        def flush(self):
            self.original_stderr.flush()
    
    # 应用过滤器（可选，如果需要完全静默）
    # sys.stderr = StderrFilter()


def get_clean_summary(success_count, failed_count, total_count, duplicates_count, unique_count):
    """生成简洁的统计摘要"""
    
    success_rate = (success_count / total_count * 100) if total_count > 0 else 0
    
    summary = f"""
{'='*60}
C6H6 异构体枚举统计
{'='*60}
  总候选结构: {total_count}
  成功处理: {success_count} ({success_rate:.1f}%)
  处理失败: {failed_count}
  重复去重: {duplicates_count}
{'-'*60}
  最终唯一结构: {unique_count} 个有效C6H6异构体
{'='*60}
"""
    return summary


def show_detailed_errors(show=False):
    """控制是否显示详细错误信息"""
    global SHOW_DETAILED_ERRORS
    SHOW_DETAILED_ERRORS = show


# 初始化时自动调用
setup_silent_mode()
