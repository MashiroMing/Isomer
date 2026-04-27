# predefined_molecules.py - 预定义分子库
# 包含常用的有机分子，可以初始化到分子库中

import os
import sys
from typing import List, Dict, Any

# 支持作为主程序运行
if __name__ == "__main__":
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, project_root)
    from database.molecule_library import molecule_library
else:
    from .molecule_library import molecule_library

def get_predefined_molecules() -> List[Dict[str, Any]]:
    """获取预定义分子列表"""
    return [
        # === 烷烃类 ===
        {
            "smiles": "CC",
            "name": "乙烷",
            "description": "最简单的烷烃分子",
            "category": "predefined",
            "tags": ["烷烃", "有机化学"]
        },
        {
            "smiles": "CCC",
            "name": "丙烷",
            "description": "三碳烷烃",
            "category": "predefined", 
            "tags": ["烷烃"]
        },
        {
            "smiles": "CCCC",
            "name": "丁烷",
            "description": "四碳烷烃，有正丁烷和异丁烷两种异构体",
            "category": "predefined",
            "tags": ["烷烃", "异构体"]
        },
        {
            "smiles": "CCCCCC",
            "name": "己烷",
            "description": "六碳烷烃",
            "category": "predefined",
            "tags": ["烷烃"]
        },
        {
            "smiles": "CCCCCCCC",
            "name": "辛烷",
            "description": "八碳烷烃，汽油的参考标准",
            "category": "predefined",
            "tags": ["烷烃", "燃料"]
        },
        
        # === 烯烃类 ===
        {
            "smiles": "C=C",
            "name": "乙烯",
            "description": "最简单的烯烃，植物激素",
            "category": "predefined",
            "tags": ["烯烃", "植物激素"]
        },
        {
            "smiles": "CC=C",
            "name": "丙烯",
            "description": "三碳烯烃",
            "category": "predefined",
            "tags": ["烯烃"]
        },
        {
            "smiles": "C=CC",
            "name": "1-丙烯",
            "description": "丙烯的1号位置双键",
            "category": "predefined",
            "tags": ["烯烃"]
        },
        {
            "smiles": "C=CCCC",
            "name": "1-戊烯",
            "description": "五碳单烯烃",
            "category": "predefined",
            "tags": ["烯烃"]
        },
        
        # === 炔烃类 ===
        {
            "smiles": "C#C",
            "name": "乙炔",
            "description": "最简单的炔烃",
            "category": "predefined",
            "tags": ["炔烃"]
        },
        {
            "smiles": "CC#C",
            "name": "丙炔",
            "description": "三碳炔烃",
            "category": "predefined",
            "tags": ["炔烃"]
        },
        {
            "smiles": "C#CCCC",
            "name": "1-戊炔",
            "description": "五碳炔烃",
            "category": "predefined",
            "tags": ["炔烃"]
        },
        
        # === 芳香烃类 ===
        {
            "smiles": "c1ccccc1",
            "name": "苯",
            "description": "最简单的芳香烃，平面六元环",
            "category": "predefined",
            "tags": ["芳香烃", "苯环"]
        },
        {
            "smiles": "c1ccc(cc1)C",
            "name": "甲苯",
            "description": "苯环上连接甲基",
            "category": "predefined",
            "tags": ["芳香烃", "取代苯"]
        },
        {
            "smiles": "c1ccc(cc1)CC",
            "name": "乙苯",
            "description": "苯环上连接乙基",
            "category": "predefined",
            "tags": ["芳香烃", "取代苯"]
        },
        {
            "smiles": "c1ccc(cc1)C(C)C",
            "name": "异丙苯",
            "description": "苯环上连接异丙基",
            "category": "predefined",
            "tags": ["芳香烃", "取代苯"]
        },
        {
            "smiles": "c1ccc(cc1)C(C)(C)C",
            "name": "叔丁苯",
            "description": "苯环上连接叔丁基",
            "category": "predefined",
            "tags": ["芳香烃", "取代苯"]
        },
        {
            "smiles": "c1ccc2c(c1)cccc2",
            "name": "萘",
            "description": "双环芳香烃",
            "category": "predefined",
            "tags": ["芳香烃", "多环"]
        },
        {
            "smiles": "c1ccc2c(c1)ccc3c2cccc3",
            "name": "蒽",
            "description": "三环芳香烃",
            "category": "predefined",
            "tags": ["芳香烃", "多环"]
        },
        
        # === 醇类 ===
        {
            "smiles": "CO",
            "name": "甲醇",
            "description": "最简单的醇类",
            "category": "predefined",
            "tags": ["醇", "含氧有机物"]
        },
        {
            "smiles": "CCO",
            "name": "乙醇",
            "description": "酒精的主要成分",
            "category": "predefined",
            "tags": ["醇", "含氧有机物"]
        },
        {
            "smiles": "CC(C)O",
            "name": "异丙醇",
            "description": "消毒酒精",
            "category": "predefined",
            "tags": ["醇", "含氧有机物"]
        },
        {
            "smiles": "CCCCO",
            "name": "正丁醇",
            "description": "四碳醇",
            "category": "predefined",
            "tags": ["醇", "含氧有机物"]
        },
        {
            "smiles": "c1ccc(cc1)O",
            "name": "苯酚",
            "description": "芳香醇",
            "category": "predefined",
            "tags": ["醇", "芳香烃", "含氧有机物"]
        },
        
        # === 醛类 ===
        {
            "smiles": "C=O",
            "name": "甲醛",
            "description": "最简单的醛类",
            "category": "predefined",
            "tags": ["醛", "含氧有机物"]
        },
        {
            "smiles": "CC=O",
            "name": "乙醛",
            "description": "二碳醛",
            "category": "predefined",
            "tags": ["醛", "含氧有机物"]
        },
        {
            "smiles": "CCC=O",
            "name": "丙醛",
            "description": "三碳醛",
            "category": "predefined",
            "tags": ["醛", "含氧有机物"]
        },
        {
            "smiles": "c1ccc(cc1)C=O",
            "name": "苯甲醛",
            "description": "芳香醛",
            "category": "predefined",
            "tags": ["醛", "芳香烃", "含氧有机物"]
        },
        
        # === 酮类 ===
        {
            "smiles": "CC(=O)C",
            "name": "丙酮",
            "description": "最简单的酮类，常用溶剂",
            "category": "predefined",
            "tags": ["酮", "含氧有机物", "溶剂"]
        },
        {
            "smiles": "CC(=O)CC",
            "name": "丁酮",
            "description": "四碳酮",
            "category": "predefined",
            "tags": ["酮", "含氧有机物"]
        },
        {
            "smiles": "CC(=O)CCC",
            "name": "戊酮",
            "description": "五碳酮",
            "category": "predefined",
            "tags": ["酮", "含氧有机物"]
        },
        {
            "smiles": "c1ccc(cc1)C(=O)C",
            "name": "苯乙酮",
            "description": "芳香酮",
            "category": "predefined",
            "tags": ["酮", "芳香烃", "含氧有机物"]
        },
        
        # === 羧酸类 ===
        {
            "smiles": "C(=O)O",
            "name": "甲酸",
            "description": "最简单的羧酸",
            "category": "predefined",
            "tags": ["羧酸", "含氧有机物"]
        },
        {
            "smiles": "CC(=O)O",
            "name": "乙酸",
            "description": "醋的主要成分",
            "category": "predefined",
            "tags": ["羧酸", "含氧有机物"]
        },
        {
            "smiles": "CCC(=O)O",
            "name": "丙酸",
            "description": "三碳羧酸",
            "category": "predefined",
            "tags": ["羧酸", "含氧有机物"]
        },
        {
            "smiles": "c1ccc(cc1)C(=O)O",
            "name": "苯甲酸",
            "description": "芳香羧酸",
            "category": "predefined",
            "tags": ["羧酸", "芳香烃", "含氧有机物"]
        },
        
        # === 酯类 ===
        {
            "smiles": "CC(=O)OC",
            "name": "乙酸甲酯",
            "description": "水果香味物质",
            "category": "predefined",
            "tags": ["酯", "含氧有机物", "香料"]
        },
        {
            "smiles": "CC(=O)OCC",
            "name": "乙酸乙酯",
            "description": "常用溶剂",
            "category": "predefined",
            "tags": ["酯", "含氧有机物", "溶剂"]
        },
        {
            "smiles": "c1ccc(cc1)C(=O)OC",
            "name": "苯甲酸甲酯",
            "description": "芳香酯",
            "category": "predefined",
            "tags": ["酯", "芳香烃", "含氧有机物"]
        },
        
        # === 卤代烃类 ===
        {
            "smiles": "CCl",
            "name": "氯甲烷",
            "description": "最简单的卤代烷",
            "category": "predefined",
            "tags": ["卤代烃", "含氯"]
        },
        {
            "smiles": "C(Cl)(Cl)(Cl)Cl",
            "name": "四氯化碳",
            "description": "曾经常用的溶剂",
            "category": "predefined",
            "tags": ["卤代烃", "含氯", "溶剂"]
        },
        {
            "smiles": "c1ccc(cc1)Cl",
            "name": "氯苯",
            "description": "芳香卤代烃",
            "category": "predefined",
            "tags": ["卤代烃", "芳香烃", "含氯"]
        },
        
        # === 含氮化合物 ===
        {
            "smiles": "CN",
            "name": "甲胺",
            "description": "最简单的胺类",
            "category": "predefined",
            "tags": ["胺", "含氮有机物"]
        },
        {
            "smiles": "CCN",
            "name": "乙胺",
            "description": "二碳胺",
            "category": "predefined",
            "tags": ["胺", "含氮有机物"]
        },
        {
            "smiles": "c1ccc(cc1)N",
            "name": "苯胺",
            "description": "芳香胺",
            "category": "predefined",
            "tags": ["胺", "芳香烃", "含氮有机物"]
        },
        {
            "smiles": "c1c[nH]cn1",
            "name": "咪唑",
            "description": "五元杂环",
            "category": "predefined",
            "tags": ["杂环", "含氮有机物"]
        },
        {
            "smiles": "c1ncccc1",
            "name": "吡啶",
            "description": "六元杂环",
            "category": "predefined",
            "tags": ["杂环", "含氮有机物"]
        },
        
        # === 特殊结构 ===
        {
            "smiles": "C1=CC2C1C2",
            "name": "杜瓦苯",
            "description": "苯的价键异构体",
            "category": "predefined",
            "tags": ["苯异构体", "价键异构体"]
        },
        {
            "smiles": "C=C1C=CC=C1",
            "name": "富烯",
            "description": "苯的价键异构体",
            "category": "predefined",
            "tags": ["苯异构体", "价键异构体"]
        },
        {
            "smiles": "C12C3C1C1C2C31",
            "name": "棱晶烷",
            "description": "苯的笼状异构体",
            "category": "predefined",
            "tags": ["苯异构体", "笼状结构"]
        }
    ]

def initialize_predefined_molecules() -> Dict[str, Any]:
    """初始化预定义分子到分子库中"""
    predefined = get_predefined_molecules()
    success_count = 0
    failed_imports = []
    
    for mol_data in predefined:
        result = molecule_library.add_molecule(
            smiles=mol_data["smiles"],
            name=mol_data["name"],
            description=mol_data["description"],
            category=mol_data["category"],
            tags=mol_data["tags"]
        )
        
        if result["success"]:
            success_count += 1
        else:
            failed_imports.append({
                "name": mol_data["name"],
                "error": result["error"]
            })
    
    return {
        "success": True,
        "imported_count": success_count,
        "failed_count": len(failed_imports),
        "failed_imports": failed_imports
    }

def check_and_initialize_predefined():
    """检查是否需要初始化预定义分子"""
    # 检查分子库是否为空
    stats = molecule_library.get_statistics()
    if stats["total_molecules"] == 0:
        print("分子库为空，正在初始化预定义分子...")
        result = initialize_predefined_molecules()
        
        print(f"初始化完成:")
        print(f"  成功导入: {result['imported_count']} 个分子")
        print(f"  导入失败: {result['failed_count']} 个分子")
        
        if result['failed_count'] > 0:
            print("失败的分子:")
            for failed in result['failed_imports']:
                print(f"  {failed['name']}: {failed['error']}")
        
        return result
    else:
        print(f"分子库已包含 {stats['total_molecules']} 个分子，跳过初始化")
        return {"success": True, "imported_count": 0, "failed_count": 0, "failed_imports": []}