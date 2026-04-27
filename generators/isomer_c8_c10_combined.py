# isomer_c8_c10_combined.py
# 合并C8-C10的所有异构体生成模块

import os
import sys

# 添加项目根目录到 Python 路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.mol_utils import process_smiles_to_data
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from typing import List, Dict, Any

# 导入通用烷烃生成算法
from generators.alkane_tree_generator import generate_alkane_isomers, get_theoretical_count

# ==============================================================================
# C8H18 异构体生成（18种异构体）
# ==============================================================================

class C8IsomerGenerator:
    """C8H18异构体生成器"""
    
    def __init__(self):
        self.expected_count = 18
    
    def generate_all_c8h18_isomers(self) -> List[str]:
        """生成所有C8H18异构体 - 使用树结构算法"""
        return generate_alkane_isomers(8)


# ==============================================================================
# C9H20 异构体生成（35种异构体）
# ==============================================================================

class C9IsomerGenerator:
    """C9H20异构体生成器"""
    
    def __init__(self):
        self.expected_count = 35
    
    def generate_all_c9h20_isomers(self) -> List[str]:
        """生成所有C9H20异构体 - 使用树结构算法"""
        return generate_alkane_isomers(9)


# ==============================================================================
# C10H22 异构体生成（75种异构体）
# ==============================================================================

class C10IsomerGenerator:
    """C10H22异构体生成器"""
    
    def __init__(self):
        self.expected_count = 75
    
    def generate_all_c10h22_isomers(self) -> List[str]:
        """生成所有C10H22异构体 - 使用树结构算法"""
        return generate_alkane_isomers(10)


# ==============================================================================
# 统一入口函数
# ==============================================================================

def generate_isomers_c8_c10(formula: str) -> List[Dict[str, Any]]:
    """
    统一的C8-C10异构体生成入口
    
    Args:
        formula: 分子式（如 'C8H18', 'C9H20', 'C10H22'）
        
    Returns:
        处理后的异构体数据列表或错误信息
    """
    if formula == "C8H18":
        generator = C8IsomerGenerator()
        smiles_list = generator.generate_all_c8h18_isomers()
        result = process_smiles_to_data(smiles_list)
        theory = get_theoretical_count(8)
        print(f"C8H18生成: {len(smiles_list)}/{theory} 种异构体")
        return result
    
    elif formula == "C9H20":
        generator = C9IsomerGenerator()
        smiles_list = generator.generate_all_c9h20_isomers()
        result = process_smiles_to_data(smiles_list)
        theory = get_theoretical_count(9)
        print(f"C9H20生成: {len(smiles_list)}/{theory} 种异构体")
        return result
    
    elif formula == "C10H22":
        generator = C10IsomerGenerator()
        smiles_list = generator.generate_all_c10h22_isomers()
        result = process_smiles_to_data(smiles_list)
        theory = get_theoretical_count(10)
        print(f"C10H22生成: {len(smiles_list)}/{theory} 种异构体")
        return result
    
    else:
        return [{"error": f"此模块支持 C8H18, C9H20, C10H22，输入为 {formula}"}]
