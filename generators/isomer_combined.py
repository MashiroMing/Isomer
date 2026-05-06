# isomer_combined.py
# 合并的有机同分异构体生成器
# 处理非饱和链烃：烯烃、炔烃、环烃、含氧化合物等
# 饱和链烃请使用 alkane_tree_generator.py

import os
import sys
import re
from typing import List, Dict, Any

# 添加项目根目录到 Python 路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils.mol_utils import process_smiles_to_data
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# ==============================================================================
# 辅助函数：图转SMILES
# ==============================================================================

def graph_to_smiles(G) -> str:
    """将networkx图转换为SMILES字符串"""
    try:
        mol = graph_to_rdkit_mol(G)
        if mol:
            return Chem.MolToSmiles(mol, isomericSmiles=False)
    except:
        pass
    return None

def graph_to_rdkit_mol(G) -> Any:
    """将networkx图转换为RDKit分子对象"""
    from rdkit.Chem import AddHs, RemoveHs
    n = G.number_of_nodes()
    if n == 0:
        return None
    
    # 创建空分子
    mol = Chem.RWMol(Chem.MolFromSmiles(''))
    for _ in range(n):
        mol.AddAtom(Chem.Atom(6))  # 碳原子
    
    # 添加键
    for u, v, d in G.edges(data=True):
        bond_type = d.get('bond_type', 'single')
        if bond_type == 'double':
            mol.AddBond(int(u), int(v), Chem.BondType.DOUBLE)
        elif bond_type == 'triple':
            mol.AddBond(int(u), int(v), Chem.BondType.TRIPLE)
        else:
            mol.AddBond(int(u), int(v), Chem.BondType.SINGLE)
    
    return RemoveHs(mol)

def generate_cyclic_isomers(n_c: int, n_h: int) -> List[str]:
    """使用内置逻辑生成环状异构体"""
    # 计算不饱和度，确定双键数
    saturated_h = 2 * n_c + 2
    unsaturation = (saturated_h - n_h) // 2
    
    if unsaturation == 0:
        # 环烷烃
        return [f'C{n_c-1}1C{n_c-1}1'] if n_c >= 3 else []
    
    # 环烯烃：简单生成（数据库覆盖主要情况）
    return []

# ==============================================================================
# C1-C5 非饱和烃数据库（炔烃等无法用算法生成的）
# ==============================================================================

STRUCTURE_DATABASE_C1_C5 = {
    # C2 烯烃/炔烃
    'C2H4': ['C=C'],
    'C2H2': ['C#C'],

    # C3 烯烃/炔烃/环烃
    'C3H6': ['C=CC', 'C1CC1'],
    'C3H4': ['C#CC', 'C=C=C'],

    # C4 烯烃/炔烃/环烃
    'C4H8': [
        'C=CCC', 'CC=CC', 'C=C(C)C',
        'C1CCC1', 'C1CC(C)1'
    ],
    'C4H6': [
        'C=CC=C', 'C=C=CC', 'C#CCC', 'CC#CC', 'C1=CCC1'
    ],

    # C5 烯烃/炔烃/环烃
    'C5H10': ['C=CCCC', 'CC=CCC', 'C=C(C)CC', 'C=CC(C)C', 'C1CCCC1', 'C1(C)CCC1'],
    'C5H8': [
        'C#CCCC', 'CC#CCC', 'CC(C)C#C',       # 炔烃
        'C=CCCC=C', 'C=C(C)C=C', r'C/C=C\C=C', r'C/C=C/C=C',  # 二烯烃
        'C1=CCCC1', 'CC1=CCC1', 'CC1C=CC1', 'C=C1CCC1',       # 环状
        'CC1=C(C)C1', 'CC1=CC(C)C1', 'CC(C)C1=CC1', 'C1(CC)=C C1',
        'CCC1=CC1', 'C=C1C(C)C1', 'C=C=C1CC1', 'C=CC1CC1',
        'C12CC1CC2'  # 螺[2,2]戊烷
    ],
}

# ==============================================================================
# C1-C5 含氧有机物数据库
# ==============================================================================

STRUCTURE_DATABASE_C1_C5_HETERO = {
    # C1
    'CH2O': ['C=O'],              # 甲醛
    'CH2O2': ['C(=O)O'],          # 甲酸
    'C1H2O': ['C=O'],             # 甲醛（兼容格式）
    'C1H2O2': ['C(=O)O'],         # 甲酸（兼容格式）
    'C1H4O': ['CO'],              # 甲醇

    # C2
    'C2H6O': ['CCO', 'COC'],      # 乙醇、二甲醚
    'C2H4O': ['CC=O', 'C1CO1'],  # 乙醛、环氧乙烷
    'C2H4O2': ['CC(=O)O'],       # 乙酸

    # C3
    'C3H8O': ['CCCO', 'CC(C)O', 'CCOC'],      # 丙醇、醚
    'C3H6O': ['CCC=O', 'CC(=O)C', 'C1CCCO1'], # 丙醛、丙酮、氧杂环丁烷
    'C3H6O2': ['CCC(=O)O', 'COC(=O)C', 'CCOC=O', 'C=C(O)CO'],

    # C4
    'C4H10O': [
        'CCCCO', 'CCC(O)C', 'CC(C)CO', 'CC(C)(O)C',  # 醇
        'CCOCC', 'CC(C)OC'                            # 醚
    ],
    'C4H8O2': [
        'CCCC(=O)O', 'CCC(=O)OC', 'CCOC(=O)C',
        'COC(=O)CCC', 'C1CCC(=O)O1'   # 丁酸、酯、内酯
    ],
}

# ==============================================================================
# C6 含氧有机物数据库
# ==============================================================================

STRUCTURE_DATABASE_C6_HETERO = {
    'C6H14O': [
        'CCCCCC(O)', 'CCCCC(O)C', 'CCCC(O)CC', 'CC(C)CCC(O)C',
        'CC(C)(C)CCCO', 'CCCCCOCC', 'CCC(C)OC(C)C'
    ],
    'C6H12O': [
        'CCCCCC=O', 'CCCCC(=O)C', 'CCCC(=O)CC',
        'C1CCCCC1O', 'C1CCOC(C)C1', 'COC(C)CCCC'
    ],
    'C6H10O2': [
        'CCCCCC(=O)O', 'CC(=O)OC(C)CC', 'C1CC(C)C(=O)OC1',
    ],
    'C6H6O': ['c1ccccc1O'],  # 苯酚
}

# ==============================================================================
# C7+ 预定义异构体数据库（包含主要不饱和化合物）
# ==============================================================================

STRUCTURE_DATABASE_C7_PLUS = {
    # C7H16 - 庚烷 9种
    'C7H16': [
        'CCCCCCC', 'CC(C)CCCC', 'CCC(C)CCC', 'CCCC(C)CC',
        'CC(C)CC(C)C', 'CC(C)C(C)CC', 'CCC(C)C(C)C',
        'CC(C)(C)CCC', 'CC(C)(C)C(C)C'
    ],
    
    # C7H8 - 甲苯系化合物 5种（来自分子库）
    'C7H8': [
        'c1ccc(cc1)C',      # 甲苯 (苯环+甲基)
        'C=C1C=CC=CC1',     # 1,3,5-庚三烯
        'C=C1C=CC=C1C',     # 4-乙烯基-1-甲基-1,3-环己二烯
        'C=CC=CC#CC',       # 1-庚烯-3,5-二炔
        'C=CC1=CC=CC1',     # 1-乙烯基-1,3-环己二烯
    ],
    
    # C7H12 - 庚烯/庚炔
    'C7H12': [
        'C=CCCCCC', 'CC=CCCCC', 'CCC=CCCC', 'CCCC=CCC',
        'CC(C)=CCCC', 'CCC(C)=CCC', 'CCCC(C)=CC',
        'C#CCCCCC', 'CC#CCCCC', 'CCC#CCCC',
        'C1CCCCCC1', 'C1CCC(C)CC1', 'C1CCCCC(C)1',
    ],
    'C8H18': [  # 辛烷 18种
        'CCCCCCCC', 'CC(C)CCCCC', 'CCC(C)CCCC', 'CCCC(C)CCC',
        'CC(C)CC(C)CC', 'CC(C)C(C)CCC', 'CCC(C)C(C)CC', 'CCCC(C)C(C)C',
        'CC(C)(C)CCCC', 'CC(C)(C)C(C)CC', 'CC(C)(C)CC(C)C', 'CC(C)C(C)(C)CC',
        'CCC(C)(C)CCC', 'C(CC(C)C)CC(C)C', 'CC(C)CC(C)(C)C', 'CC(C)(C)CC(C)C',
        'CC(C)(C)C(C)(C)C', 'CC(C)(C)(C)CCC'
    ],
    'C9H20': [  # 壬烷 35种
        'CCCCCCCCC', 'CC(C)CCCCCC', 'CCC(C)CCCCC', 'CCCC(C)CCCC',
        'CC(C)CC(C)CCC', 'CC(C)C(C)CCCC', 'CCC(C)C(C)CCC', 'CCCC(C)C(C)CC',
        'CC(C)(C)CCCCC', 'CC(C)(C)C(C)CCC', 'CC(C)(C)CC(C)CC', 'CC(C)C(C)(C)CCC',
        'CCC(C)(C)CCCC', 'C(CC(C)C)CCC(C)C', 'CC(C)CC(C)(C)CC', 'CC(C)(C)CC(C)(C)C',
        'CC(C)(C)C(C)(C)CC', 'CCC(C)(C)C(C)CC', 'CCCC(C)(C)(C)CC', 'CC(C)(C)(C)CCCC',
        'CC(C)C(CC)CCC', 'CC(C)(CC)CCCC', 'C(C(C)C)CCCC', 'C(C)CCCCCCC',
        'CCCCCCCC(C)C', 'CCCCC(C)CCC', 'CCCC(C)CCCC(C)C', 'CC(C)CCCCC(C)C',
        'CCCC(C)(C)CCC', 'CC(C)(C)CCCC(C)C', 'CC(C)CCCC(C)(C)C', 'CCC(C)CCC(C)(C)C',
        'CCCC(C)CC(C)(C)C', 'CC(C)CC(C)CC(C)C'
    ],
}

# ==============================================================================
# 绘图配置（C6H6 特殊处理）
# ==============================================================================

def setup_drawing_options():
    """设置绘图选项"""
    from rdkit.Chem.Draw.MolDrawing import DrawingOptions
    opts = DrawingOptions()
    opts.atomLabelFontSize = 400
    opts.dblBondLengthFrac = 0.8
    opts.includeAtomNumbers = True
    opts.dotsPerAngstrom = 1000
    return opts

def enhance_benzene_structure(mol, kekulize=False):
    """增强苯结构的芳香性表达"""
    if mol is None:
        return None
    try:
        if kekulize:
            Chem.Kekulize(mol, clearAromaticFlags=True)
        else:
            for bond in mol.GetBonds():
                if bond.GetIsAromatic():
                    bond.SetBondType(Chem.BondType.AROMATIC)
            for atom in mol.GetAtoms():
                if atom.GetIsAromatic():
                    atom.SetIsAromatic(True)
                    atom.SetHybridization(Chem.HybridizationType.SP2)
    except Exception as e:
        print(f"芳香性处理警告: {e}")
    return mol

def classify_c6h6_structure(smiles: str) -> dict:
    """智能分类C6H6结构类型"""
    classification = {}
    if smiles == 'c1ccccc1':
        classification = {'StructureType': '芳香烃 (苯)', 'Aromatic': True, 'RingSystem': '6元芳香环', 'Category': 'benzene'}
    elif '#' in smiles:
        bond_count = smiles.count('#')
        classification = {'StructureType': '多炔烃' if bond_count >= 2 else '炔烃', 'Aromatic': False, 'Category': 'polyalkyne' if bond_count >= 2 else 'alkyne'}
    elif '=' in smiles and '1' in smiles:
        classification = {'StructureType': '环状烯烃', 'Aromatic': False, 'Category': 'cycloalkene'}
    elif '=' in smiles:
        classification = {'StructureType': '共轭烯烃', 'Aromatic': False, 'Category': 'conjugated_alkene'}
    elif '1' in smiles:
        classification = {'StructureType': '环状烷烃', 'Aromatic': False, 'Category': 'cycloalkane'}
    else:
        classification = {'StructureType': '烷烃/饱和环', 'Aromatic': False, 'Category': 'alkane'}
    return classification

def is_pure_c6h6_hydrocarbon(mol) -> bool:
    """检查分子是否为纯C6H6碳氢化合物"""
    try:
        mol_with_h = Chem.AddHs(mol)
        carbon_count = sum(1 for atom in mol_with_h.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count != 6:
            return False
        formula = rdMolDescriptors.CalcMolFormula(mol_with_h)
        return formula == "C6H6"
    except:
        return False

def generate_c6h6_general(formula: str) -> list:
    """C6H6纯碳氢化合物算法"""
    if formula != "C6H6":
        return [{"error": "此模块专门用于C6H6"}]
    
    # 预定义C6H6异构体列表
    verified_c6h6_isomers = [
        'c1ccccc1',   # 苯
        'C1=CC=CC=C1', # 苯（Kekule结构，等价于c1ccccc1）
        'C#CCCC#C',   # 1,5-己二炔
    ]
    
    all_results = []
    seen_smiles = set()  # 去重
    for i, smiles in enumerate(verified_c6h6_isomers):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None or not is_pure_c6h6_hydrocarbon(mol):
                continue
            processed_data = process_smiles_to_data([smiles])
            if processed_data:
                for item in processed_data:
                    canon_smiles = item.get('SMILES', '')
                    if canon_smiles in seen_smiles:
                        continue
                    seen_smiles.add(canon_smiles)
                    item.update(classify_c6h6_structure(smiles))
                    item['IsomerNumber'] = i + 1
                    item['TotalCount'] = len(verified_c6h6_isomers)
                    all_results.append(item)
        except:
            continue
    
    return all_results

# ==============================================================================
# 链烯烃/炔烃生成器
# ==============================================================================

class ChainAlkeneGenerator:
    """链状烯烃/炔烃生成器"""
    
    @staticmethod
    def generate_chain_alkenes(n_c: int, n_h: int) -> List[str]:
        """生成链烯烃异构体"""
        if n_h != 2 * n_c:  # 非烯烃
            return []
        
        if n_c < 2:
            return []
        
        alkenes = []
        # 双键在不同位置
        for pos in range(1, n_c):
            # 生成线性链烯烃
            smiles = 'C' * pos + '=' + 'C' * (n_c - pos)
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                canon_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
                alkenes.append(canon_smiles)
        
        return list(set(alkenes))
    
    @staticmethod
    def generate_chain_alkynes(n_c: int, n_h: int) -> List[str]:
        """生成链炔烃异构体"""
        if n_h != 2 * n_c - 2:  # 非炔烃
            return []
        
        if n_c < 2:
            return []
        
        alkynes = []
        # 三键在不同位置
        for pos in range(1, n_c):
            smiles = 'C' * pos + '#' + 'C' * (n_c - pos)
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                canon_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
                alkynes.append(canon_smiles)
        
        return list(set(alkynes))

# ==============================================================================
# 通用异构体生成器（C7+）
# ==============================================================================

class GeneralIsomerGenerator:
    """通用有机同分异构体生成器"""
    
    def __init__(self):
        self.max_results = 100
        self.chain_gen = ChainAlkeneGenerator()
    
    def generate_isomers_general(self, formula: str) -> List[Dict[str, Any]]:
        """通用同分异构体生成主函数"""
        # 先检查预定义数据库（饱和烷烃）
        if formula in STRUCTURE_DATABASE_C7_PLUS:
            return process_smiles_to_data(STRUCTURE_DATABASE_C7_PLUS[formula])
        
        # 解析分子式
        parsed = self._parse_formula(formula)
        if "error" in parsed:
            return [parsed]
        
        n_c, n_h, n_o = parsed["n_c"], parsed["n_h"], parsed.get("n_o", 0)
        
        all_isomers = []
        
        # 生成链烯烃/炔烃（使用算法）
        if n_o == 0:
            chain_alkenes = self.chain_gen.generate_chain_alkenes(n_c, n_h)
            chain_alkynes = self.chain_gen.generate_chain_alkynes(n_c, n_h)
            all_isomers.extend(chain_alkenes)
            all_isomers.extend(chain_alkynes)
        
        # 生成环状结构（使用环烯烃生成器）
        if n_c >= 3 and n_o == 0:
            cyclic_smiles = generate_cyclic_isomers(n_c, n_h)
            all_isomers.extend(cyclic_smiles)
        
        # 生成含氧有机物
        if n_o > 0:
            all_isomers.extend(self._generate_oxygen_chains(n_c, n_h, n_o))
        
        return self._deduplicate_and_validate(all_isomers, formula)[:self.max_results]
    
    def _parse_formula(self, formula: str) -> Dict[str, Any]:
        import re
        formula = formula.strip().upper().replace(' ', '')
        match_c = re.search(r'C(\d+)', formula)
        match_h = re.search(r'H(\d*)', formula)
        match_o = re.search(r'O(\d*)', formula)
        if not match_c:
            return {"error": "分子式格式错误，缺少碳原子数"}
        try:
            n_c = int(match_c.group(1))
            n_h = int(match_h.group(1)) if match_h and match_h.group(1) else 0
            n_o = int(match_o.group(1)) if match_o and match_o.group(1) else 0
        except ValueError:
            return {"error": "原子数解析失败"}
        return {"n_c": n_c, "n_h": n_h, "n_o": n_o}
    
    def _generate_oxygen_chains(self, n_c: int, n_h: int, n_o: int) -> List[str]:
        oxygen_chains = []
        if n_o == 1:
            if n_h == 2 * n_c + 2:  # 饱和醇
                for pos in range(1, n_c + 1):
                    if pos == 1:
                        oxygen_chains.append('CO' + 'C' * (n_c - 1))
                    elif pos == n_c:
                        oxygen_chains.append('C' * (n_c - 1) + 'CO')
                    else:
                        oxygen_chains.append('C' * (pos - 1) + 'C(O)' + 'C' * (n_c - pos))
            if n_h == 2 * n_c + 2:  # 醚
                for pos in range(1, n_c):
                    oxygen_chains.append('C' * pos + 'OC' + 'C' * (n_c - pos - 1))
            if n_h == 2 * n_c:  # 醛/酮
                if n_c >= 2:
                    oxygen_chains.append('CC=O' + 'C' * max(0, n_c - 2))
                if n_c >= 3:
                    for pos in range(2, n_c):
                        oxygen_chains.append('C' * (pos - 1) + 'C(=O)' + 'C' * (n_c - pos))
        return oxygen_chains
    
    def _deduplicate_and_validate(self, smiles_list: List[str], formula: str) -> List[Dict[str, Any]]:
        unique_smiles = set()
        valid_isomers = []
        for smiles in smiles_list:
            if smiles in unique_smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            try:
                processed_data = process_smiles_to_data([smiles])
                if processed_data:
                    valid_isomers.extend(processed_data)
                    unique_smiles.add(smiles)
            except:
                continue
        return valid_isomers

def generate_isomers_general(formula: str) -> List[Dict[str, Any]]:
    """通用同分异构体生成接口"""
    generator = GeneralIsomerGenerator()
    return generator.generate_isomers_general(formula)

# ==============================================================================
# 统一入口函数
# ==============================================================================

def generate_isomers_c1_c5(formula: str):
    """C1-C5 非饱和烃入口"""
    match = re.match(r'C(\d+)H(\d+)', formula)
    if not match:
        return {"error": f"无效的分子式: {formula}"}
    
    n_c = int(match.group(1))
    n_h = int(match.group(2))
    
    # 饱和烷烃使用树算法
    if n_h == 2 * n_c + 2 and n_c >= 1:
        from generators.alkane_tree_generator import generate_alkane_isomers
        smiles_list = generate_alkane_isomers(n_c)
        return process_smiles_to_data(smiles_list)
    
    # 非饱和烃：先尝试使用算法生成
    if n_c >= 3:
        all_smiles = []
        
        # 链烯烃
        if n_h == 2 * n_c:
            chain_gen = ChainAlkeneGenerator()
            all_smiles.extend(chain_gen.generate_chain_alkenes(n_c, n_h))
        
        # 链炔烃
        if n_h == 2 * n_c - 2:
            chain_gen = ChainAlkeneGenerator()
            all_smiles.extend(chain_gen.generate_chain_alkynes(n_c, n_h))
        
        # 环状结构（使用环烯烃生成器）
        cyclic_smiles = generate_cyclic_isomers(n_c, n_h)
        all_smiles.extend(cyclic_smiles)
        
        if all_smiles:
            return process_smiles_to_data(all_smiles)
    
    # 回退到数据库
    if formula not in STRUCTURE_DATABASE_C1_C5:
        return {"error": f"分子式 {formula} 不在 C1-C5 非饱和烃预设列表中。"}
    
    return process_smiles_to_data(STRUCTURE_DATABASE_C1_C5[formula])

def generate_isomers_c1_c5_hetero(formula: str):
    """C1-C5 含氧有机物入口"""
    if formula not in STRUCTURE_DATABASE_C1_C5_HETERO:
        return {"error": f"分子式 {formula} 不在 C1-C5 含氧结构预设列表中。"}
    return process_smiles_to_data(STRUCTURE_DATABASE_C1_C5_HETERO[formula])

def generate_isomers_c6(formula: str):
    """C6 异构体统一入口"""
    if formula == "C6H6":
        return generate_c6h6_general(formula)
    if formula == "C6H14":
        from generators.alkane_tree_generator import generate_alkane_isomers
        return process_smiles_to_data(generate_alkane_isomers(6))
    
    # C6H12 环己烷或链烯烃
    if formula == "C6H12":
        smiles = ['C1CCCCC1', 'C=CCCCC', 'CC=CCCC', 'CCC=CCC', 'CCC(C)=CC']
        return process_smiles_to_data(smiles)
    
    # C6H10 炔烃或环烯烃
    if formula == "C6H10":
        all_smiles = []
        # 链炔烃
        chain_gen = ChainAlkeneGenerator()
        all_smiles.extend(chain_gen.generate_chain_alkynes(6, 10))
        # 环烯烃
        cyclo_smiles = ['C1=CCCCC1']
        for sm in cyclo_smiles:
            mol = Chem.MolFromSmiles(sm)
            if mol:
                all_smiles.append(sm)
        if all_smiles:
            return process_smiles_to_data(all_smiles)
        return {"error": f"无法生成 {formula} 异构体"}
    
    # C6H8 单环二烯烃
    if formula == "C6H8":
        cyclo_smiles = ['C1=CC=CCC1', 'C1=CCC=CC1']
        return process_smiles_to_data(cyclo_smiles)
    
    if formula in STRUCTURE_DATABASE_C6_HETERO:
        return process_smiles_to_data(STRUCTURE_DATABASE_C6_HETERO[formula])
    return {"error": f"分子式 {formula} 不在支持列表中。"}

def generate_isomers_c6_hetero(formula: str):
    """C6 含氧有机物入口"""
    if formula not in STRUCTURE_DATABASE_C6_HETERO:
        return {"error": f"分子式 {formula} 不在 C6 含氧结构预设列表中。"}
    return process_smiles_to_data(STRUCTURE_DATABASE_C6_HETERO[formula])

def generate_isomers_combined(formula: str):
    """
    统一的异构体生成入口
    自动根据分子式选择合适的生成策略
    """
    formula = formula.strip().upper().replace(' ', '')
    
    # 规范化分子式
    if formula == 'C':
        formula = 'C1H4'
    elif re.match(r'^C[H,O]', formula):
        formula = 'C1' + formula[1:]
    
    # 解析碳数
    match = re.match(r'C(\d+)', formula)
    if not match:
        return {"error": "分子式格式错误"}
    
    n_c = int(match.group(1))
    contains_o = 'O' in formula
    
    # 分发到对应函数
    if 1 <= n_c <= 5:
        if contains_o:
            return generate_isomers_c1_c5_hetero(formula)
        return generate_isomers_c1_c5(formula)
    elif n_c == 6:
        return generate_isomers_c6(formula)
    else:  # C7+
        return generate_isomers_general(formula)