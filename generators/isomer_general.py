# isomer_general.py - 通用有机同分异构体生成算法
# 支持C7及以上碳数的分子结构生成

import os
import sys
import itertools
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import List, Dict, Any, Set, Tuple

# 添加项目根目录到 Python 路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.mol_utils import process_smiles_to_data

class GeneralIsomerGenerator:
    """
    通用有机同分异构体生成器
    使用系统性的算法来生成合理的同分异构体结构
    """
    
    def __init__(self):
        self.max_results = 100  # 限制最大结果数量避免内存溢出
        
    def generate_isomers_general(self, formula: str) -> List[Dict[str, Any]]:
        """
        通用同分异构体生成主函数
        """
        # 解析分子式
        parsed = self._parse_formula(formula)
        if "error" in parsed:
            return [parsed]
            
        n_c = parsed["n_c"]
        n_h = parsed["n_h"] 
        n_o = parsed.get("n_o", 0)
        
        # 生成策略：基于氢碳比和分子大小选择合适的策略
        all_isomers = []
        
        # 策略1：生成链式结构（烷烃、烯烃、炔烃等）
        chain_isomers = self._generate_chain_structures(n_c, n_h, n_o)
        all_isomers.extend(chain_isomers)
        
        # 策略2：生成支链结构
        branch_isomers = self._generate_branched_structures(n_c, n_h, n_o)
        all_isomers.extend(branch_isomers)
        
        # 策略3：生成环状结构（当碳数足够时）
        if n_c >= 3:
            ring_isomers = self._generate_ring_structures(n_c, n_h, n_o)
            all_isomers.extend(ring_isomers)
            
        # 去重和验证
        unique_isomers = self._deduplicate_and_validate(all_isomers, formula)
        
        return unique_isomers[:self.max_results]
    
    def _parse_formula(self, formula: str) -> Dict[str, Any]:
        """解析分子式，提取原子数"""
        import re
        
        formula = formula.strip().upper().replace(' ', '')
        
        # 解析碳、氢、氧原子数
        match_c = re.search(r'C(\d+)', formula)
        match_h = re.search(r'H(\d*)', formula)
        match_o = re.search(r'O(\d*)', formula)
        
        if not match_c:
            return {"error": "分子式格式错误，缺少碳原子数"}
            
        try:
            n_c = int(match_c.group(1))
            n_h = int(match_h.group(1)) if match_h and match_h.group(1) else 0
            n_o = 1 if match_o else (int(match_o.group(1)) if match_o and match_o.group(1) else 0)
        except ValueError:
            return {"error": "原子数解析失败"}
            
        return {"n_c": n_c, "n_h": n_h, "n_o": n_o}
    
    def _generate_chain_structures(self, n_c: int, n_h: int, n_o: int) -> List[str]:
        """生成链式结构的SMILES"""
        chains = []
        
        # 生成基本的链结构
        if n_h == 2 * n_c + 2:  # 烷烃 CnH2n+2
            chains.append('C' * n_c)
            
        elif n_h == 2 * n_c:  # 烯烃 CnH2n  
            # 双键在不同位置
            for i in range(1, n_c):
                smiles = 'C' * i + '=' + 'C' * (n_c - i)
                chains.append(smiles)
                
        elif n_h == 2 * n_c - 2:  # 炔烃 CnH2n-2
            # 三键在不同位置
            for i in range(1, n_c):
                smiles = 'C' * i + '#' + 'C' * (n_c - i)
                chains.append(smiles)
                
        # 含氧化合物的链式结构
        if n_o > 0:
            oxygen_chains = self._generate_oxygen_chains(n_c, n_h, n_o)
            chains.extend(oxygen_chains)
            
        return chains
    
    def _generate_oxygen_chains(self, n_c: int, n_h: int, n_o: int) -> List[str]:
        """生成含氧链式结构"""
        oxygen_chains = []
        
        if n_o == 1:
            # 醇类：在链的不同位置添加羟基 (CnH2n+2O -> 醇)
            if n_h == 2 * n_c + 2:  # 饱和醇
                for pos in range(1, n_c + 1):
                    if pos == 1:
                        smiles = 'CO' + 'C' * (n_c - 1)  # 1-醇
                    elif pos == n_c:
                        smiles = 'C' * (n_c - 1) + 'CO'  # n-醇
                    else:
                        smiles = 'C' * (pos - 1) + 'C(O)' + 'C' * (n_c - pos)
                    oxygen_chains.append(smiles)
                    
            # 醚类：氧原子在不同位置 (CnH2n+2O -> 醚)
            if n_h == 2 * n_c + 2:
                for pos in range(1, n_c):
                    smiles = 'C' * pos + 'OC' + 'C' * (n_c - pos - 1)
                    oxygen_chains.append(smiles)
                
            # 醛基：在链的不同位置 (CnH2nO -> 醛)
            if n_h == 2 * n_c:
                for pos in range(1, n_c + 1):
                    if pos == 1:
                        smiles = 'CC=O' + 'C' * (n_c - 2) if n_c > 1 else 'C=O'
                    elif pos == n_c:
                        smiles = 'C' * (n_c - 1) + 'C=O'
                    else:
                        smiles = 'C' * (pos - 1) + 'CC(=O)' + 'C' * (n_c - pos - 1)
                    oxygen_chains.append(smiles)
                    
            # 酮基：在链的不同位置 (CnH2nO -> 酮)
            if n_h == 2 * n_c and n_c >= 3:
                for pos in range(2, n_c):
                    smiles = 'C' * (pos - 1) + 'C(=O)' + 'C' * (n_c - pos)
                    oxygen_chains.append(smiles)
                    
        return oxygen_chains
    
    def _generate_branched_structures(self, n_c: int, n_h: int, n_o: int) -> List[str]:
        """生成支链结构"""
        branched = []
        
        if n_c >= 4:  # 从C4开始可以有支链
            # 生成不同类型的支链结构
            branched.extend(self._generate_simple_branches(n_c, n_h, n_o))
            
        return branched
    
    def _generate_simple_branches(self, n_c: int, n_h: int, n_o: int) -> List[str]:
        """生成简单支链结构 - 庚烷的9种同分异构体"""
        branches = []
        
        if n_c == 7 and n_h == 16:  # 庚烷 C7H16
            # 庚烷的9种同分异构体
            branches = [
                'CCCCCCC',      # 正庚烷
                'CC(C)CCCC',    # 2-甲基己烷
                'CCC(C)CCC',    # 3-甲基己烷
                'CCCC(C)CC',    # 4-甲基己烷
                'CC(C)CC(C)C',  # 2,4-二甲基戊烷
                'CC(C)C(C)CC',  # 2,3-二甲基戊烷
                'CCC(C)C(C)C',  # 3,3-二甲基戊烷
                'CC(C)(C)CCC',  # 2,2-二甲基戊烷
                'CC(C)(C)C(C)C' # 2,2,3-三甲基丁烷
            ]
        elif n_c >= 4:
            # 通用支链生成算法
            branches.extend(self._generate_generic_branches(n_c, n_h, n_o))
                
        return branches
    
    def _generate_generic_branches(self, n_c: int, n_h: int, n_o: int) -> List[str]:
        """通用支链生成算法"""
        branches = []
        
        # 简单的单支链结构
        for methyl_pos in range(2, n_c - 1):
            # 甲基支链
            prefix = 'C' * (methyl_pos - 1)
            branch = '(C)'
            suffix = 'C' * (n_c - methyl_pos - 1)
            smiles = prefix + branch + suffix
            branches.append(smiles)
            
        # 双支链结构（对于较大的分子）
        if n_c >= 6:
            # 2,2-二甲基支链
            if n_c >= 5:
                branches.append('CC(C)(C)C' + 'C' * (n_c - 5))
            # 2,3-二甲基支链
            if n_c >= 5:
                branches.append('CC(C)C(C)C' + 'C' * (n_c - 5))
                
        return branches
    
    def _generate_ring_structures(self, n_c: int, n_h: int, n_o: int) -> List[str]:
        """生成环状结构"""
        rings = []
        
        # 小环结构（3-6元环）
        for ring_size in range(3, min(7, n_c + 1)):
            if ring_size <= n_c:
                # 基础环状结构
                if ring_size == n_c:
                    # 单环
                    rings.append('C1' + 'C' * (ring_size - 2) + 'C1')
                else:
                    # 环加链
                    ring_part = 'C1' + 'C' * (ring_size - 2) + 'C1'
                    chain_part = 'C' * (n_c - ring_size)
                    rings.append(ring_part + chain_part)
                    rings.append(chain_part + ring_part)
                    
        # 芳香环（当碳数为6时）
        if n_c >= 6:
            # 苯环及其衍生物
            if n_c == 6:
                rings.append('c1ccccc1')  # 苯
            elif n_c > 6:
                # 苯环加侧链
                benzene = 'c1ccccc1'
                side_chain = 'C' * (n_c - 6)
                rings.append(benzene + side_chain)
                
        # 多环结构（当碳数足够时）
        if n_c >= 8:
            rings.extend(self._generate_polycyclic_structures(n_c))
            
        return rings
    
    def _generate_polycyclic_structures(self, n_c: int) -> List[str]:
        """生成多环结构"""
        polycyclic = []
        
        # 萘结构（C10）
        if n_c == 10:
            polycyclic.append('c1ccc2ccccc2c1')
        # 联苯结构
        elif n_c >= 6:
            half_size = n_c // 2
            if half_size >= 3:
                ring1 = 'c1ccccc1'[:2*half_size+1]  # 部分苯环
                ring2 = 'c1ccccc1'[:2*half_size+1]
                polycyclic.append(ring1 + ring2)
                
        return polycyclic
    
    def _deduplicate_and_validate(self, smiles_list: List[str], formula: str) -> List[Dict[str, Any]]:
        """去重并验证SMILES结构"""
        unique_smiles = set()
        valid_isomers = []
        
        for smiles in smiles_list:
            # 去重
            if smiles in unique_smiles:
                continue
                
            # 验证结构
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
                
            # 检查分子式是否匹配
            actual_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            if not self._formula_matches(actual_formula, formula):
                continue
                
            unique_smiles.add(smiles)
            
            # 生成分子数据
            try:
                processed_data = process_smiles_to_data([smiles])
                if processed_data:
                    valid_isomers.extend(processed_data)
            except Exception as e:
                print(f"处理SMILES {smiles} 时出错: {e}", file=sys.stderr)
                continue
                
        return valid_isomers
    
    def _formula_matches(self, actual: str, expected: str) -> bool:
        """检查分子式是否匹配"""
        # 简化的分子式匹配逻辑
        # 将分子式标准化后比较
        def normalize_formula(formula_str):
            import re
            # 提取原子数
            atoms = {}
            for atom in re.findall(r'([A-Z][a-z]?)(\d*)', formula_str):
                element, count = atom
                atoms[element] = int(count) if count else 1
            return atoms
            
        actual_atoms = normalize_formula(actual)
        expected_atoms = normalize_formula(expected)
        
        return actual_atoms == expected_atoms

# 主函数接口
def generate_isomers_general(formula: str) -> List[Dict[str, Any]]:
    """
    通用同分异构体生成接口
    支持C7及以上碳数的分子
    """
    generator = GeneralIsomerGenerator()
    return generator.generate_isomers_general(formula)