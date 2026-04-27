# structure_filter.py - 结构筛选模块
# 根据用户指定的结构特征筛选同分异构体

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from typing import List, Dict, Any, Set
import re
import sys

class StructureFilter:
    """
    结构筛选器
    根据用户指定的结构特征筛选同分异构体
    """
    
    def __init__(self):
        pass
    
    def filter_isomers(self, isomers: List[Dict[str, Any]], 
                     structure_prefs: Dict[str, str]) -> List[Dict[str, Any]]:
        """
        根据结构偏好筛选异构体
        
        Args:
            isomers: 原始异构体列表
            structure_prefs: 结构偏好字典
                - ring_type: 'all', 'acyclic', 'monocyclic', 'polycyclic'
                - bond_type: 'all', 'single', 'double', 'triple', 'mixed'
                - aromatic: 'all', 'none', 'benzene', 'aromatic'
                - functional_groups: 'all', 'alcohol', 'carbonyl', 'carboxyl', 'ester'
        
        Returns:
            筛选后的异构体列表
        """
        if not isomers:
            return isomers
            
        # 如果所有选项都是'all'，直接返回
        if all(pref == 'all' for pref in structure_prefs.values()):
            return isomers
        
        filtered_isomers = []
        
        for isomer in isomers:
            smiles = isomer.get('SMILES', '')
            if not smiles:
                continue
                
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                    
                if self._matches_preferences(mol, structure_prefs):
                    # 更新结构类型描述
                    enhanced_type = self._get_enhanced_structure_type(mol, isomer.get('StructureType', ''))
                    isomer['StructureType'] = enhanced_type
                    filtered_isomers.append(isomer)
                    
            except Exception as e:
                print(f"筛选分子时出错 {smiles}: {e}", file=sys.stderr)
                continue
        
        return filtered_isomers
    
    def _matches_preferences(self, mol: Chem.Mol, prefs: Dict[str, str]) -> bool:
        """检查分子是否符合用户偏好"""
        
        # 检查环结构
        if 'ring_type' in prefs and not self._matches_ring_type(mol, prefs['ring_type']):
            return False
            
        # 检查键类型
        if 'bond_type' in prefs and not self._matches_bond_type(mol, prefs['bond_type']):
            return False
            
        # 检查芳香性
        if 'aromatic' in prefs and not self._matches_aromatic_type(mol, prefs['aromatic']):
            return False
            
        # 检查官能团
        if 'functional_groups' in prefs and not self._matches_functional_groups(mol, prefs['functional_groups']):
            return False
            
        return True
    
    def _matches_ring_type(self, mol: Chem.Mol, ring_type: str) -> bool:
        """检查环类型匹配"""
        if ring_type == 'all':
            return True
            
        ring_info = mol.GetRingInfo()
        ring_count = ring_info.NumRings()
        
        if ring_type == 'acyclic':
            return ring_count == 0
        elif ring_type == 'monocyclic':
            return ring_count == 1
        elif ring_type == 'polycyclic':
            return ring_count > 1
            
        return True
    
    def _matches_bond_type(self, mol: Chem.Mol, bond_type: str) -> bool:
        """检查键类型匹配"""
        if bond_type == 'all':
            return True
            
        single_bonds = double_bonds = triple_bonds = 0
        
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                single_bonds += 1
            elif bond.GetBondType() == Chem.BondType.DOUBLE:
                double_bonds += 1
            elif bond.GetBondType() == Chem.BondType.TRIPLE:
                triple_bonds += 1
        
        if bond_type == 'single':
            return double_bonds == 0 and triple_bonds == 0 and single_bonds > 0
        elif bond_type == 'double':
            return double_bonds > 0 and triple_bonds == 0
        elif bond_type == 'triple':
            return triple_bonds > 0
        elif bond_type == 'mixed':
            return (double_bonds > 0 and single_bonds > 0) or \
                   (triple_bonds > 0 and single_bonds > 0) or \
                   (double_bonds > 0 and triple_bonds > 0)
            
        return True
    
    def _matches_aromatic_type(self, mol: Chem.Mol, aromatic_type: str) -> bool:
        """检查芳香性匹配"""
        if aromatic_type == 'all':
            return True
            
        aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        aromatic_bonds = sum(1 for bond in mol.GetBonds() if bond.GetIsAromatic())
        
        if aromatic_type == 'none':
            return aromatic_atoms == 0 and aromatic_bonds == 0
        elif aromatic_type == 'benzene':
            # 检查是否为纯碳的苯环结构（整个分子就是C6H6）
            return self._is_pure_carbon_benzene(mol)
        elif aromatic_type == 'aromatic':
            return aromatic_atoms > 0 or aromatic_bonds > 0
            
        return True
    
    def _has_benzene_ring(self, mol: Chem.Mol) -> bool:
        """检查是否含有苯环（包括杂环芳香环）"""
        try:
            # 寻找6元环
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if len(ring) == 6:
                    # 检查是否为芳香环
                    is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
                    if is_aromatic:
                        return True
            return False
        except:
            return False
    
    def _has_pure_carbon_benzene_ring(self, mol: Chem.Mol) -> bool:
        """检查是否含有纯碳苯环（只有碳原子的苯环）"""
        try:
            # 寻找6元环
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if len(ring) == 6:
                    # 检查是否为芳香环
                    is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
                    if is_aromatic:
                        # 检查环中的原子是否全部为碳原子
                        all_carbon = all(mol.GetAtomWithIdx(idx).GetSymbol() == 'C' for idx in ring)
                        if all_carbon:
                            return True
            return False
        except:
            return False
    
    def _is_pure_carbon_benzene(self, mol: Chem.Mol) -> bool:
        """检查是否为纯碳的苯环结构（整个分子就是C6H6）"""
        try:
            # 首先检查分子式
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            if formula != 'C6H6':
                return False
            
            # 检查是否为单个6元芳香环且全为碳原子
            ring_info = mol.GetRingInfo()
            rings = ring_info.AtomRings()
            
            # 必须只有一个环，且为6元环
            if len(rings) != 1 or len(rings[0]) != 6:
                return False
            
            # 检查环中所有原子是否为碳且芳香
            ring = rings[0]
            is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
            all_carbon = all(mol.GetAtomWithIdx(idx).GetSymbol() == 'C' for idx in ring)
            
            # 检查所有原子是否都在环中
            atom_count = mol.GetNumAtoms()
            if atom_count != 6:
                return False
            
            return is_aromatic and all_carbon
        except:
            return False
    
    def _matches_functional_groups(self, mol: Chem.Mol, fg_type: str) -> bool:
        """检查官能团匹配"""
        if fg_type == 'all':
            return True
            
        # 定义官能团的SMARTS模式
        functional_patterns = {
            'alcohol': '[OH]',                    # 羟基
            'carbonyl': '[CX3]=[O]',             # 羰基 (醛酮)
            'carboxyl': 'C(=O)[OX1H0-,OX2H+]', # 羧基
            'ester': 'C(=O)[OX2H0]',            # 酯基
        }
        
        if fg_type not in functional_patterns:
            return True
            
        pattern = functional_patterns[fg_type]
        substructure = Chem.MolFromSmarts(pattern)
        
        if substructure is None:
            return True
            
        return mol.HasSubstructMatch(substructure)
    
    def _get_enhanced_structure_type(self, mol: Chem.Mol, original_type: str) -> str:
        """获取增强的结构类型描述"""
        features = []
        
        # 环结构信息
        ring_info = mol.GetRingInfo()
        ring_count = ring_info.NumRings()
        if ring_count == 0:
            features.append("无环")
        elif ring_count == 1:
            features.append("单环")
        else:
            features.append("多环")
        
        # 芳香性信息
        if self._is_pure_carbon_benzene(mol):
            features.append("含苯环")
        elif self._has_pure_carbon_benzene_ring(mol):
            features.append("含苯环取代基")
        elif self._has_benzene_ring(mol):
            features.append("杂环芳香")
        elif sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()) > 0:
            features.append("芳香")
        
        # 键类型信息
        bond_types = set()
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                bond_types.add("双键")
            elif bond.GetBondType() == Chem.BondType.TRIPLE:
                bond_types.add("三键")
        
        if bond_types:
            features.extend(bond_types)
        
        # 官能团信息
        if mol.HasSubstructMatch(Chem.MolFromSmarts('[OH]')):
            features.append("羟基")
        if mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3]=[O]')):
            features.append("羰基")
        if mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)[OX1H0-,OX2H+]')):
            features.append("羧基")
        if mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)[OX2H0]')):
            features.append("酯基")
        
        # 组合特征
        if features:
            enhanced = f"{original_type} | {', '.join(features)}"
        else:
            enhanced = original_type
            
        return enhanced
    
    def get_structure_statistics(self, isomers: List[Dict[str, Any]]) -> Dict[str, Any]:
        """获取异构体的结构统计信息"""
        if not isomers:
            return {}
            
        stats = {
            'total_count': len(isomers),
            'ring_types': {'无环': 0, '单环': 0, '多环': 0},
            'bond_types': {'单键': 0, '双键': 0, '三键': 0, '混合': 0},
            'aromatic_types': {'无芳香': 0, '含苯环': 0, '其他芳香': 0},
            'functional_groups': {'羟基': 0, '羰基': 0, '羧基': 0, '酯基': 0}
        }
        
        for isomer in isomers:
            smiles = isomer.get('SMILES', '')
            if not smiles:
                continue
                
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                    
                # 统计环类型
                ring_info = mol.GetRingInfo()
                ring_count = ring_info.NumRings()
                if ring_count == 0:
                    stats['ring_types']['无环'] += 1
                elif ring_count == 1:
                    stats['ring_types']['单环'] += 1
                else:
                    stats['ring_types']['多环'] += 1
                
                # 统计键类型
                has_single = has_double = has_triple = False
                for bond in mol.GetBonds():
                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        has_single = True
                    elif bond.GetBondType() == Chem.BondType.DOUBLE:
                        has_double = True
                    elif bond.GetBondType() == Chem.BondType.TRIPLE:
                        has_triple = True
                
                if has_double and has_single:
                    stats['bond_types']['混合'] += 1
                elif has_triple and (has_single or has_double):
                    stats['bond_types']['混合'] += 1
                elif has_double:
                    stats['bond_types']['双键'] += 1
                elif has_triple:
                    stats['bond_types']['三键'] += 1
                else:
                    stats['bond_types']['单键'] += 1
                
                # 统计芳香性
                if self._is_pure_carbon_benzene(mol):
                    stats['aromatic_types']['含苯环'] += 1
                elif self._has_pure_carbon_benzene_ring(mol):
                    stats['aromatic_types']['含苯环'] += 1
                elif self._has_benzene_ring(mol):
                    stats['aromatic_types']['其他芳香'] += 1
                elif sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()) > 0:
                    stats['aromatic_types']['其他芳香'] += 1
                else:
                    stats['aromatic_types']['无芳香'] += 1
                
                # 统计官能团
                if mol.HasSubstructMatch(Chem.MolFromSmarts('[OH]')):
                    stats['functional_groups']['羟基'] += 1
                if mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3]=[O]')):
                    stats['functional_groups']['羰基'] += 1
                if mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)[OX1H0-,OX2H+]')):
                    stats['functional_groups']['羧基'] += 1
                if mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)[OX2H0]')):
                    stats['functional_groups']['酯基'] += 1
                    
            except Exception as e:
                print(f"统计分子结构时出错 {smiles}: {e}", file=sys.stderr)
                continue
        
        return stats