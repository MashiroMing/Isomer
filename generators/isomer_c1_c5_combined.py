# isomer_c1_c5_combined.py
# 合并C1-C5的所有异构体生成模块

import os
import sys
import re

# 添加项目根目录到 Python 路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.mol_utils import process_smiles_to_data
# 避免循环导入：在函数内部延迟导入
# from .alkane_tree_generator import generate_alkane_isomers

# ==============================================================================
# 核心结构数据库 (C1-C5) - 纯碳氢化合物
# ==============================================================================

STRUCTURE_DATABASE_C1_C5 = {
    # --------------------------------------------------------------------------
    # 说明：饱和烷烃 (CnH2n+2) 已迁移至 alkane_tree_generator.py
    # 以下为非饱和烃（烯烃/炔烃/环烃）
    # --------------------------------------------------------------------------

    # C2 烯烃/炔烃
    'C2H4': ['C=C'],
    'C2H2': ['C#C'],

    # C3 烯烃/炔烃/环烃
    'C3H6': [
        'C=CC',
        'C1CC1'
    ],
    'C3H4': [
        'C#CC',
        'C=C=C'
    ],

    # C4 烯烃/炔烃/环烃
    'C4H8': [
        'C=CCC',
        'CC=CC',
        'C=C(C)C',
        'C1CCC1',
        'C1CC(C)1'
    ],
    'C4H6': [
        'C=CC=C',
        'C=C=CC',
        'C#CCC',
        'CC#CC',
        'C1=CCC1'
    ],

    # C5 烯烃/炔烃/环烃
    'C5H10': [             # 戊烯/环戊烷 (6种常见结构)
        'C=CCCC',
        'CC=CCC',
        'C=C(C)CC',
        'C=CC(C)C',
        'C1CCCC1',
        'C1(C)CCC1',
    ],
    'C5H8': [              # 戊二烯/戊炔/环状结构 (20 种结构)
        # --- 炔烃 (3种) ---
        'C#CCCC',          # 1-戊炔
        'CC#CCC',          # 2-戊炔
        'CC(C)C#C',        # 3-甲基-1-丁炔

        # --- 二烯烃 (4种结构异构体 + E/Z 两种) ---
        'C=CCCC=C',         # 1,4-戊二烯
        'C=C(C)C=C',        # 2-甲基-1,3-丁二烯 (异戊二烯)
        r'C/C=C\C=C',        # (3Z)-1,3-戊二烯
        r'C/C=C/C=C',        # (3E)-1,3-戊二烯

        # --- 环状结构 (13 种) ---
        'C1=CCCC1',         # 环戊烯
        'CC1=CCC1',         # 1-甲基环丁烯
        'CC1C=CC1',         # 3-甲基环丁烯
        'C=C1CCC1',         # 亚甲基环丁烷
        'CC1=C(C)C1',       # 1,2-二甲基环丙烯
        'CC1=CC(C)C1',      # 1,3-二甲基环丙烯
        'CC(C)C1=CC1',      # 3,3-二甲基环丙烯
        'C1(CC)=C C1',      # 1-乙基环丙烯
        'CCC1=CC1',         # 3-乙基环丙烯
        'C=C1C(C)C1',       # 1-甲基-2-亚甲基环丙烷
        'C=C=C1CC1',        # 亚乙烯基环丙烷
        'C=CC1CC1',         # 乙烯基环丙烷
        'C12CC1CC2'         # 螺[2,2]戊烷
    ],
}


def generate_isomers_c1_c5(formula):
    """
    针对 C1 到 C5 的纯碳氢分子式，生成并处理同分异构体。

    策略：
    - 饱和烷烃 (CnH2n+2): 使用树结构生成算法 (alkane_tree_generator)
    - 其他 (烯烃/炔烃/环烃): 使用硬编码数据库
    
    Args:
        formula: 分子式（如 'C3H8'）

    Returns:
        处理后的异构体数据列表
    """
    # 解析分子式
    match = re.match(r'C(\d+)H(\d+)', formula)
    if not match:
        return {"error": f"无效的分子式: {formula}"}

    n_c = int(match.group(1))
    n_h = int(match.group(2))

    # 判断是否为饱和烷烃: H = 2C + 2
    if n_h == 2 * n_c + 2:
        # 使用树结构生成算法（通用烷烃生成器）
        # 延迟导入避免循环依赖
        from .alkane_tree_generator import generate_alkane_isomers
        smiles_list = generate_alkane_isomers(n_c)
        return process_smiles_to_data(smiles_list)

    # 非饱和烷烃（烯烃、炔烃、环烃等），使用硬编码数据
    if formula not in STRUCTURE_DATABASE_C1_C5:
        # 饱和烷烃使用树算法，非饱和烃列出可用的
        return {"error": f"分子式 {formula} 不在 C1-C5 非饱和烃预设列表中。\n饱和烷烃 (CnH2n+2) 已自动使用算法生成。\n当前支持的非饱和烃: {', '.join(sorted(STRUCTURE_DATABASE_C1_C5.keys()))}"}

    smiles_list = STRUCTURE_DATABASE_C1_C5[formula]

    # 调用核心处理函数
    return process_smiles_to_data(smiles_list)


# ==============================================================================
# 核心结构数据库 (C1-C5 含氧结构)
# ==============================================================================

STRUCTURE_DATABASE_C1_C5_HETERO = {
    # --------------------------------------------------------------------------
    # C1
    # --------------------------------------------------------------------------
    'CH2O': ['C=O'],              # 甲醛 (Aldehyde)
    'CH2O2': ['C(=O)O'],         # 甲酸 (Carboxylic Acid)
    'C1H2O': ['C=O'],             # 甲醛 (Aldehyde) - 兼容CH2O输入
    'C1H2O2': ['C(=O)O'],        # 甲酸 (Carboxylic Acid) - 兼容CH2O2输入
    'C1H4O': ['CO'],               # 甲醇 (Alcohol)
    
    # --------------------------------------------------------------------------
    # C2
    # --------------------------------------------------------------------------
    'C2H6O': [
        'CCO',         # 乙醇 (Alcohol)
        'COC'          # 二甲醚 (Ether)
    ],
    'C2H4O': [
        'CC=O',        # 乙醛 (Aldehyde)
        'C1CO1'        # 环氧乙烷 (Ether/Epoxide)
    ],
    'C2H4O2': [
        'CC(=O)O'      # 乙酸 (Carboxylic Acid)
    ],
    
    # --------------------------------------------------------------------------
    # C3
    # --------------------------------------------------------------------------
    'C3H8O': [
        'CCCO',        # 丙-1-醇
        'CC(C)O',      # 丙-2-醇
        'CCOC'         # 甲基乙基醚
    ],
    'C3H6O': [
        'CCC=O',       # 丙醛 (Aldehyde)
        'CC(=O)C',     # 丙酮 (Ketone)
        'C1CCCO1'      # 氧杂环丁烷 (Ether/cyclic)
    ],
    'C3H6O2': [        
        'CCC(=O)O',    # 丙酸 (Carboxylic Acid)
        'COC(=O)C',    # 乙酸甲酯 (Ester)
        'CCOC=O',      # 甲酸乙酯 (Ester)
        'C=C(O)CO'     # 羟基丙烯醛（enol/diol）
    ],
    
    # --------------------------------------------------------------------------
    # C4
    # --------------------------------------------------------------------------
    'C4H10O': [
        'CCCCO',        # 丁-1-醇
        'CCC(O)C',      # 丁-2-醇
        'CC(C)CO',      # 2-甲基丙-1-醇
        'CC(C)(O)C',    # 叔丁醇
        'CCOCC',        # 二乙醚
        'CC(C)OC'       # 甲基异丙醚
    ],
    'C4H8O2': [
        'CCCC(=O)O',    # 丁酸 (Carboxylic Acid)
        'CCC(=O)OC',    # 丙酸甲酯 (Ester)
        'CCOC(=O)C',    # 乙酸乙酯 (Ester)
        'COC(=O)CCC',   # 甲酸丙酯 (Ester)
        'C1CCC(=O)O1'   # 丁内酯 (Lactone/Cyclic Ester)
    ],
}


def generate_isomers_c1_c5_hetero(formula):
    """
    针对 C1 到 C5 的含氧分子式，生成并处理同分异构体。
    
    Args:
        formula: 分子式（如 'C2H6O'）
        
    Returns:
        处理后的异构体数据列表
    """
    if formula not in STRUCTURE_DATABASE_C1_C5_HETERO:
        return {"error": f"分子式 {formula} 不在 C1-C5 含氧结构预设列表中。目前支持的分子式类型有: {', '.join(STRUCTURE_DATABASE_C1_C5_HETERO.keys())}"}

    smiles_list = STRUCTURE_DATABASE_C1_C5_HETERO[formula]
    
    return process_smiles_to_data(smiles_list)


def generate_isomers_c1_to_c5(formula):
    """
    统一的C1-C5异构体生成入口
    
    Args:
        formula: 分子式
        
    Returns:
        处理后的异构体数据列表或错误信息
    """
    # 先尝试纯碳氢化合物
    result = generate_isomers_c1_c5(formula)
    
    # 如果没找到，尝试含氧化合物
    if isinstance(result, dict) and 'error' in result:
        result = generate_isomers_c1_c5_hetero(formula)
    
    return result
