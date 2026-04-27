# mol_utils.py

from rdkit import Chem
from rdkit.Chem import AllChem
import sys

# 静默RDKit警告
try:
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
except:
    pass

# --- RDKit Standardizer 兼容性设置 ---
try:
    from rdkit.Chem.MolStandardize import Standardizer as RDKitStandardizer
    standardizer_instance = RDKitStandardizer()
except ImportError:
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
        # 创建一个正确的清理函数包装器
        def cleanup_wrapper(mol):
            if mol is None:
                return None
            return rdMolStandardize.Cleanup(mol)
        standardizer_instance = type('DummyStandardizer', (), {
            'standardize': lambda self, mol: cleanup_wrapper(mol),
            'cleanup': lambda self, mol: cleanup_wrapper(mol)
        })()
    except Exception as e:
        print(f"RDKit 规范化工具初始化失败: {e}", file=sys.stderr)
        standardizer_instance = None
except Exception as e:
    print(f"RDKit 规范化工具初始化失败: {e}", file=sys.stderr)
    standardizer_instance = None

# --- NEW: SMARTS patterns for functional groups ---
FUNCTIONAL_GROUP_SMARTS = {
    # 含氧官能团 (按优先级，复杂/特异性高的优先)
    "羧基 (Carboxylic Acid)": Chem.MolFromSmarts('[CX3](=O)[OH]'),
    "酯基 (Ester)": Chem.MolFromSmarts('[CX3](=O)[O][#6]'),
    "醛基 (Aldehyde)": Chem.MolFromSmarts('[CX3H1](=O)[#6]'),
    "酮基 (Ketone)": Chem.MolFromSmarts('[CX3](=O)[#6][#6]'), # 确保是酮，而非醛/酸/酯
    "醚键 (Ether)": Chem.MolFromSmarts('[CX4]-[O]-[CX4]'),
    "醇/酚羟基 (Hydroxyl)": Chem.MolFromSmarts('[OX2H]'),
    
    # 烃类结构
    "芳香烃 (苯)": Chem.MolFromSmarts('c1ccccc1'),
    "炔烃": Chem.MolFromSmarts('C#C'),
    "烯烃": Chem.MolFromSmarts('C=C'),
}
# 优先级列表
FUNCTIONAL_GROUP_PRIORITY = list(FUNCTIONAL_GROUP_SMARTS.keys())





def process_smiles_to_data(smiles_list):
    """
    接收 SMILES 列表，处理每个分子并返回包含坐标和结构类型的数据列表。
    新增：结构识别逻辑使用 SMARTS 优先级列表进行更精确的官能团识别。
    """
    # 如果标准化工具不可用，仍然可以处理基本功能
    results = []
    seen_smiles = set()
    
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            # 1. 结构规范化（如果可用）
            canonical_smiles = smiles  # 默认使用原始SMILES
            if standardizer_instance is not None:
                try:
                    if hasattr(standardizer_instance, 'standardize'):
                        standardized_mol = standardizer_instance.standardize(mol)
                        if standardized_mol is not None:
                            mol = standardized_mol
                    elif hasattr(standardizer_instance, 'cleanup'):
                        cleaned_mol = standardizer_instance.cleanup(mol)
                        if cleaned_mol is not None:
                            mol = cleaned_mol
                    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
                except Exception:
                    # 即使规范化失败也继续处理
                    canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
            else:
                canonical_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
                
            if canonical_smiles in seen_smiles:
                continue
            seen_smiles.add(canonical_smiles)

            # 2. 坐标生成
            mol_with_h = Chem.AddHs(mol)
            dim_type = "3D"
            try:
                if AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG()) == -1 or mol_with_h.GetConformer() is None:
                    raise ValueError("3D Embedding failed, falling back to 2D.")
            except:
                 AllChem.Compute2DCoords(mol_with_h)
                 dim_type = "2D"
            
            conf = mol_with_h.GetConformer()
            coords = []
            for i in range(mol_with_h.GetNumAtoms()):
                atom = mol_with_h.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                z_coord = f", z={pos.z:.4f}" if dim_type == "3D" else ""
                coords.append(
                    f"Atom {i+1} ({atom.GetSymbol()}): x={pos.x:.4f}, y={pos.y:.4f}{z_coord}"
                )
            coords_str = "\n".join(coords)
            
            # 3. 结构识别
            structure_type = "未知结构"
            if standardizer_instance is not None:
                try:
                    detected_groups = []
                    for name in FUNCTIONAL_GROUP_PRIORITY:
                        pat = FUNCTIONAL_GROUP_SMARTS[name]
                        if mol.HasSubstructMatch(pat):
                            # 特殊处理：防止醇/酚羟基被重复检测
                            if name == "醇/酚羟基 (Hydroxyl)":
                                # 如果已检测到酸或酯，则不报告醇/酚羟基
                                if any(g in detected_groups for g in ["羧基 (Carboxylic Acid)", "酯基 (Ester)"]):
                                    continue
                                # 如果是酚 (直接连接到芳香环)，则单独标注
                                if mol.HasSubstructMatch(Chem.MolFromSmarts('c[OH]')):
                                    detected_groups.append("酚羟基 (Phenol)")
                                else:
                                    detected_groups.append("醇羟基 (Alcohol)")
                            else:
                                 detected_groups.append(name)
                    
                    if detected_groups:
                        # 移除重复项并生成结构描述
                        structure_type = " & ".join(sorted(list(set(detected_groups))))
                    else:
                        structure_type = "烷烃/饱和环"
                except Exception:
                    structure_type = "烷烃/饱和环"  # 默认类型
            else:
                structure_type = "未识别(缺少RDKit标准化工具)"
                
            inchikey = Chem.MolToInchiKey(mol)
            
            results.append({
                "SMILES": canonical_smiles,
                "InChIKey": inchikey,
                "StructureType": structure_type,
                "Coords": coords_str,
                "Dimension": dim_type,
                "Mol": mol 
            })

        except Exception as e:
            print(f"处理SMILES {smiles} 失败: {e}", file=sys.stderr)
            continue
    
    return results