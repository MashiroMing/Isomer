# 异构体生成模块
# 提供各碳数烷烃异构体的生成功能

# 饱和链烃生成器（独立模块）
from .alkane_tree_generator import generate_alkane_isomers, get_theoretical_count

# 合并的非饱和烃/含氧有机物生成器
from .isomer_combined import (
    generate_isomers_c1_c5,
    generate_isomers_c1_c5_hetero,
    generate_isomers_c6,
    generate_c6h6_general,
    generate_isomers_c6_hetero,
    generate_isomers_general,
    generate_isomers_combined  # 统一入口
)

__all__ = [
    # 饱和链烃
    'generate_alkane_isomers',
    'get_theoretical_count',
    # 非饱和烃/含氧有机物
    'generate_isomers_c1_c5',
    'generate_isomers_c1_c5_hetero',
    'generate_isomers_c6',
    'generate_c6h6_general',
    'generate_isomers_c6_hetero',
    'generate_isomers_general',
    'generate_isomers_combined',  # 统一入口
]
