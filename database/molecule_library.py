# molecule_library.py - 分子库管理系统
# 支持用户手动添加、编辑、删除和搜索分子

import os
import json
from typing import Dict, List, Any, Optional
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import re

class MoleculeLibrary:
    """分子库管理类"""

    def __init__(self, library_file: str = None):
        if library_file is None:
            # 默认使用 data 目录中的文件
            library_file = os.path.join(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                'data',
                'molecule_library.json'
            )
        self.library_file = library_file
        self.library_data = self._load_library()
    
    def _load_library(self) -> Dict[str, Any]:
        """加载分子库数据"""
        if os.path.exists(self.library_file):
            try:
                with open(self.library_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                print(f"加载分子库失败: {e}")
                return self._get_empty_library()
        else:
            return self._get_empty_library()
    
    def _get_empty_library(self) -> Dict[str, Any]:
        """返回空的分子库结构"""
        return {
            "version": "1.0",
            "created_date": datetime.now().isoformat(),
            "last_updated": datetime.now().isoformat(),
            "categories": {
                "user_added": [],
                "predefined": [],
                "favorites": []
            },
            "total_count": 0
        }
    
    def _save_library(self) -> bool:
        """保存分子库数据"""
        try:
            self.library_data["last_updated"] = datetime.now().isoformat()
            self.library_data["total_count"] = (
                len(self.library_data["categories"]["user_added"]) +
                len(self.library_data["categories"]["predefined"]) +
                len(self.library_data["categories"]["favorites"])
            )
            
            with open(self.library_file, 'w', encoding='utf-8') as f:
                json.dump(self.library_data, f, ensure_ascii=False, indent=2)
            return True
        except IOError as e:
            print(f"保存分子库失败: {e}")
            return False
    
    def validate_molecule(self, smiles: str) -> Dict[str, Any]:
        """验证分子SMILES字符串"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"valid": False, "error": "SMILES字符串无效"}
            
            # 获取分子信息
            formula = rdMolDescriptors.CalcMolFormula(mol)
            mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
            
            # 解析分子式
            carbon_count = self._count_atoms_in_formula(formula, 'C')
            hydrogen_count = self._count_atoms_in_formula(formula, 'H')
            oxygen_count = self._count_atoms_in_formula(formula, 'O')
            
            return {
                "valid": True,
                "formula": formula,
                "molecular_weight": mol_weight,
                "carbon_count": carbon_count,
                "hydrogen_count": hydrogen_count,
                "oxygen_count": oxygen_count,
                "atom_count": mol.GetNumAtoms(),
                "bond_count": mol.GetNumBonds()
            }
        except Exception as e:
            return {"valid": False, "error": f"验证失败: {str(e)}"}
    
    def _count_atoms_in_formula(self, formula: str, element: str) -> int:
        """从分子式中统计指定元素的数量"""
        pattern = f"{element}(\\d*)"
        matches = re.findall(pattern, formula)
        total = 0
        for match in matches:
            if match == "":
                total += 1
            else:
                total += int(match)
        return total
    
    def add_molecule(self, 
                    smiles: str, 
                    name: str = "", 
                    description: str = "", 
                    category: str = "user_added",
                    tags: List[str] = None) -> Dict[str, Any]:
        """添加分子到库中"""
        
        # 验证SMILES
        validation = self.validate_molecule(smiles)
        if not validation["valid"]:
            return {"success": False, "error": validation["error"]}
        
        # 生成唯一ID
        molecule_id = f"mol_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{len(self.library_data['categories'][category])}"
        
        # 创建分子记录
        molecule_record = {
            "id": molecule_id,
            "name": name or f"Molecule_{molecule_id[-6:]}",
            "smiles": smiles,
            "formula": validation["formula"],
            "molecular_weight": validation["molecular_weight"],
            "carbon_count": validation["carbon_count"],
            "hydrogen_count": validation["hydrogen_count"],
            "oxygen_count": validation["oxygen_count"],
            "description": description,
            "category": category,
            "tags": tags or [],
            "created_date": datetime.now().isoformat(),
            "custom_properties": {}
        }
        
        # 添加到对应类别
        if category not in self.library_data["categories"]:
            self.library_data["categories"][category] = []
        
        self.library_data["categories"][category].append(molecule_record)
        
        # 保存库
        if self._save_library():
            return {"success": True, "molecule_id": molecule_id, "record": molecule_record}
        else:
            return {"success": False, "error": "保存失败"}
    
    def remove_molecule(self, molecule_id: str) -> Dict[str, Any]:
        """从库中删除分子"""
        for category, molecules in self.library_data["categories"].items():
            for i, molecule in enumerate(molecules):
                if molecule["id"] == molecule_id:
                    removed_molecule = molecules.pop(i)
                    if self._save_library():
                        return {"success": True, "removed": removed_molecule}
                    else:
                        # 如果保存失败，恢复分子
                        molecules.insert(i, removed_molecule)
                        return {"success": False, "error": "保存失败"}
        
        return {"success": False, "error": "分子不存在"}
    
    def update_molecule(self, molecule_id: str, updates: Dict[str, Any]) -> Dict[str, Any]:
        """更新分子信息"""
        for category, molecules in self.library_data["categories"].items():
            for molecule in molecules:
                if molecule["id"] == molecule_id:
                    # 保存原始信息
                    original = molecule.copy()
                    
                    # 更新允许的字段
                    updatable_fields = ["name", "description", "tags", "custom_properties"]
                    for field, value in updates.items():
                        if field in updatable_fields:
                            molecule[field] = value
                    
                    molecule["last_updated"] = datetime.now().isoformat()
                    
                    if self._save_library():
                        return {"success": True, "updated": molecule}
                    else:
                        # 如果保存失败，恢复原始信息
                        molecules[molecules.index(molecule)] = original
                        return {"success": False, "error": "保存失败"}
        
        return {"success": False, "error": "分子不存在"}
    
    def search_molecules(self, 
                        query: str = "", 
                        formula: str = "", 
                        category: str = "",
                        tags: List[str] = None,
                        carbon_range: tuple = None,
                        hydrogen_range: tuple = None) -> List[Dict[str, Any]]:
        """搜索分子"""
        results = []
        
        for cat_name, molecules in self.library_data["categories"].items():
            if category and cat_name != category:
                continue
            
            for molecule in molecules:
                # 文本搜索
                if query:
                    query_lower = query.lower()
                    if not (query_lower in molecule["name"].lower() or 
                           query_lower in molecule["description"].lower() or
                           query_lower in molecule["smiles"].lower()):
                        continue
                
                # 分子式搜索
                if formula and not molecule["formula"].startswith(formula.upper()):
                    continue
                
                # 标签搜索
                if tags:
                    if not any(tag in molecule["tags"] for tag in tags):
                        continue
                
                # 碳原子数范围
                if carbon_range:
                    c_count = molecule["carbon_count"]
                    if not (carbon_range[0] <= c_count <= carbon_range[1]):
                        continue
                
                # 氢原子数范围
                if hydrogen_range:
                    h_count = molecule["hydrogen_count"]
                    if not (hydrogen_range[0] <= h_count <= hydrogen_range[1]):
                        continue
                
                results.append(molecule)
        
        return results
    
    def get_molecule_by_id(self, molecule_id: str) -> Optional[Dict[str, Any]]:
        """根据ID获取分子"""
        for molecules in self.library_data["categories"].values():
            for molecule in molecules:
                if molecule["id"] == molecule_id:
                    return molecule
        return None
    
    def get_all_categories(self) -> List[str]:
        """获取所有类别"""
        return list(self.library_data["categories"].keys())
    
    def get_statistics(self) -> Dict[str, Any]:
        """获取库统计信息"""
        stats = {
            "total_molecules": 0,
            "by_category": {},
            "by_carbon_count": {},
            "recent_added": []
        }
        
        for category, molecules in self.library_data["categories"].items():
            stats["by_category"][category] = len(molecules)
            stats["total_molecules"] += len(molecules)
            
            # 统计碳原子数分布
            for molecule in molecules:
                c_count = molecule["carbon_count"]
                stats["by_carbon_count"][c_count] = stats["by_carbon_count"].get(c_count, 0) + 1
        
        # 获取最近添加的分子（按创建时间排序）
        all_molecules = []
        for molecules in self.library_data["categories"].values():
            all_molecules.extend(molecules)
        
        all_molecules.sort(key=lambda x: x.get("created_date", ""), reverse=True)
        stats["recent_added"] = all_molecules[:5]  # 最近5个
        
        return stats
    
    def export_molecules(self, molecule_ids: List[str], format_type: str = "json") -> Dict[str, Any]:
        """导出分子数据"""
        molecules = []
        for mol_id in molecule_ids:
            molecule = self.get_molecule_by_id(mol_id)
            if molecule:
                molecules.append(molecule)
        
        if format_type == "json":
            return {"success": True, "data": molecules, "format": "json"}
        elif format_type == "smiles":
            smiles_list = [mol["smiles"] for mol in molecules]
            return {"success": True, "data": smiles_list, "format": "smiles"}
        else:
            return {"success": False, "error": "不支持的导出格式"}
    
    def import_molecules(self, molecules_data: List[Dict[str, Any]]) -> Dict[str, Any]:
        """导入分子数据"""
        success_count = 0
        failed_imports = []
        
        for mol_data in molecules_data:
            if "smiles" not in mol_data:
                failed_imports.append({"error": "缺少SMILES", "data": mol_data})
                continue
            
            result = self.add_molecule(
                smiles=mol_data["smiles"],
                name=mol_data.get("name", ""),
                description=mol_data.get("description", ""),
                category=mol_data.get("category", "user_added"),
                tags=mol_data.get("tags", [])
            )
            
            if result["success"]:
                success_count += 1
            else:
                failed_imports.append({"error": result["error"], "data": mol_data})
        
        return {
            "success": True,
            "imported_count": success_count,
            "failed_count": len(failed_imports),
            "failed_imports": failed_imports
        }

# 全局分子库实例
molecule_library = MoleculeLibrary()

# 初始化预定义分子
if __name__ == "__main__":
    import sys
    from database.predefined_molecules import check_and_initialize_predefined
    check_and_initialize_predefined()