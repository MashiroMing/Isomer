# deepseek_database_integration.py
# DeepSeek AI 安全接入数据库模块
# 实现AI动态更新分子库，扩展项目能力

import os
import json
import hashlib
from typing import List, Dict, Any, Optional, Tuple
from datetime import datetime
import logging

# 检查依赖
try:
    from openai import OpenAI
    HAS_OPENAI = True
except ImportError:
    HAS_OPENAI = False
    print("提示: DeepSeek AI 模块不可用，请安装 openai: pip install openai")

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, AllChem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("提示: RDKit 模块不可用，请安装 rdkit: conda install -c conda-forge rdkit")

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DeepSeekDatabaseIntegration:
    """DeepSeek AI 数据库集成类
    
    功能：
    1. AI验证和去重生成的分子
    2. AI补充缺失的异构体
    3. AI生成新的分子变种
    4. 安全的数据库更新机制
    5. 操作日志和回滚功能
    """
    
    def __init__(self, api_key: str = None, library_file: str = None):
        """初始化DeepSeek数据库集成

        Args:
            api_key: DeepSeek API密钥（可选，默认从环境变量读取）
            library_file: 分子库文件路径（默认使用 data/molecule_library.json）
        """
        # 设置默认分子库路径
        if library_file is None:
            # 使用 data 目录中的文件
            library_file = os.path.join(
                os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                'data',
                'molecule_library.json'
            )
        # 检查依赖
        if not HAS_OPENAI:
            raise ImportError(
                "DeepSeek AI 模块不可用。请安装 openai 包:\n"
                "  pip install openai\n"
                "或运行安装脚本:\n"
                "  install_dependencies.bat"
            )

        if not HAS_RDKIT:
            raise ImportError(
                "RDKit 模块不可用。请安装 rdkit 包:\n"
                "  conda install -c conda-forge rdkit\n"
                "或运行安装脚本:\n"
                "  install_dependencies.bat"
            )

        if api_key is None:
            api_key = os.environ.get('DEEPSEEK_API_KEY')

        if not api_key:
            raise ValueError("DeepSeek API密钥未提供。请设置环境变量DEEPSEEK_API_KEY或传入api_key参数")

        self.client = OpenAI(api_key=api_key, base_url="https://api.deepseek.com")
        self.model = "deepseek-chat"
        self.library_file = library_file
        self.backup_dir = "database_backups"
        self.operation_log = "database_operations.log"

        # 创建备份目录
        os.makedirs(self.backup_dir, exist_ok=True)

        # 加载分子库
        self.library_data = self._load_library()
    
    def _load_library(self) -> Dict[str, Any]:
        """加载分子库"""
        if os.path.exists(self.library_file):
            try:
                with open(self.library_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except Exception as e:
                logger.error(f"加载分子库失败: {e}")
                return self._get_empty_library()
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
                "favorites": [],
                "ai_generated": [],
                "ai_verified": []
            },
            "total_count": 0
        }
    
    def _create_backup(self, operation_type: str) -> str:
        """创建数据库备份
        
        Args:
            operation_type: 操作类型（如'add', 'update', 'delete'）
            
        Returns:
            备份文件路径
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_filename = f"backup_{operation_type}_{timestamp}.json"
        backup_path = os.path.join(self.backup_dir, backup_filename)
        
        try:
            with open(backup_path, 'w', encoding='utf-8') as f:
                json.dump(self.library_data, f, ensure_ascii=False, indent=2)
            logger.info(f"备份创建成功: {backup_path}")
            return backup_path
        except Exception as e:
            logger.error(f"创建备份失败: {e}")
            raise
    
    def _log_operation(self, operation: str, details: Dict[str, Any]):
        """记录数据库操作
        
        Args:
            operation: 操作类型
            details: 操作详情
        """
        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "operation": operation,
            "details": details
        }
        
        try:
            # 追加到操作日志文件
            if os.path.exists(self.operation_log):
                with open(self.operation_log, 'r', encoding='utf-8') as f:
                    logs = json.load(f)
            else:
                logs = []
            
            logs.append(log_entry)
            
            # 只保留最近100条记录
            if len(logs) > 100:
                logs = logs[-100:]
            
            with open(self.operation_log, 'w', encoding='utf-8') as f:
                json.dump(logs, f, ensure_ascii=False, indent=2)
        except Exception as e:
            logger.error(f"记录操作日志失败: {e}")
    
    def _validate_smiles_rdkit(self, smiles: str) -> Dict[str, Any]:
        """使用RDKit验证SMILES
        
        Args:
            smiles: SMILES字符串
            
        Returns:
            验证结果字典
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"valid": False, "error": "SMILES字符串无效"}
            
            formula = rdMolDescriptors.CalcMolFormula(mol)
            mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
            
            return {
                "valid": True,
                "formula": formula,
                "molecular_weight": mol_weight,
                "carbon_count": self._count_atoms_in_formula(formula, 'C'),
                "hydrogen_count": self._count_atoms_in_formula(formula, 'H'),
                "oxygen_count": self._count_atoms_in_formula(formula, 'O'),
                "atom_count": mol.GetNumAtoms(),
                "bond_count": mol.GetNumBonds()
            }
        except Exception as e:
            return {"valid": False, "error": f"验证失败: {str(e)}"}
    
    def _count_atoms_in_formula(self, formula: str, element: str) -> int:
        """从分子式中统计指定元素的数量"""
        import re
        pattern = f"{element}(\\d*)"
        matches = re.findall(pattern, formula)
        total = 0
        for match in matches:
            if match == "":
                total += 1
            else:
                total += int(match)
        return total
    
    def _normalize_smiles(self, smiles: str) -> Optional[str]:
        """标准化SMILES字符串"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return Chem.MolToSmiles(mol, canonical=True)
        except:
            return None
    
    def ai_verify_molecules(self, smiles_list: List[str], formula: str) -> Dict[str, Any]:
        """使用DeepSeek AI验证和筛选分子
        
        Args:
            smiles_list: 待验证的SMILES列表
            formula: 目标分子式
            
        Returns:
            验证结果字典
        """
        logger.info(f"开始AI验证 {len(smiles_list)} 个分子，目标分子式: {formula}")
        
        # 先用RDKit进行基础验证
        valid_smiles = []
        invalid_count = 0
        
        for smiles in smiles_list:
            validation = self._validate_smiles_rdkit(smiles)
            if validation["valid"] and validation["formula"] == formula:
                valid_smiles.append(smiles)
            else:
                invalid_count += 1
        
        logger.info(f"RDKit验证完成: 有效 {len(valid_smiles)}, 无效 {invalid_count}")
        
        # 使用AI进行深度验证
        smiles_str = "\n".join([f"{i+1}. {smi}" for i, smi in enumerate(valid_smiles[:50])])
        
        prompt = f"""你是一个专业的有机化学家，请验证以下SMILES字符串是否对应分子式 {formula}。

待验证的SMILES列表：
{smiles_str}

任务：
1. 检查每个SMILES的分子式是否确实是 {formula}
2. 标记出所有化学上不合理或错误的SMILES
3. 返回应该保留的有效SMILES的索引（从1开始）

输出格式：
只输出应该保留的SMILES索引，用逗号分隔，例如：1,3,5,7,9
不要输出任何解释或其他文字。
"""
        
        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": "你是一个专业的有机化学专家，专门验证分子结构的正确性。"},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.1,
                max_tokens=500
            )
            
            ai_result = response.choices[0].message.content.strip()
            
            # 解析AI返回的索引
            try:
                indices = [int(x.strip()) for x in ai_result.split(',') if x.strip().isdigit()]
                verified_smiles = [valid_smiles[i-1] for i in indices if 1 <= i <= len(valid_smiles)]
            except:
                logger.warning("AI返回格式解析失败，保留所有RDKit验证通过的分子")
                verified_smiles = valid_smiles
            
            result = {
                "success": True,
                "total_count": len(smiles_list),
                "rdkit_valid": len(valid_smiles),
                "ai_verified": len(verified_smiles),
                "invalid_count": len(smiles_list) - len(verified_smiles),
                "verified_smiles": verified_smiles
            }
            
            self._log_operation("ai_verify", {
                "formula": formula,
                "input_count": len(smiles_list),
                "verified_count": len(verified_smiles)
            })
            
            logger.info(f"AI验证完成: 最终保留 {len(verified_smiles)} 个有效分子")
            return result
            
        except Exception as e:
            logger.error(f"AI验证失败: {e}")
            return {
                "success": False,
                "error": str(e),
                "verified_smiles": valid_smiles
            }
    
    def ai_supplement_isomers(self, formula: str, existing_smiles: List[str], 
                              target_count: int = 10, max_attempts: int = 3) -> Dict[str, Any]:
        """使用DeepSeek AI补充缺失的异构体
        
        Args:
            formula: 目标分子式
            existing_smiles: 已存在的SMILES列表
            target_count: 目标补充数量
            max_attempts: 最大尝试次数
            
        Returns:
            补充结果字典
        """
        logger.info(f"开始AI补充异构体: {formula}, 目标数量: {target_count}")
        
        existing_normalized = set()
        for smi in existing_smiles:
            normalized = self._normalize_smiles(smi)
            if normalized:
                existing_normalized.add(normalized)
        
        new_smiles = []
        attempt = 0
        
        while len(new_smiles) < target_count and attempt < max_attempts:
            attempt += 1
            needed = target_count - len(new_smiles)
            
            # 显示部分已存在的分子
            existing_sample = list(existing_normalized)[:10]
            existing_str = "\n".join(existing_sample) if existing_sample else "（无）"
            
            prompt = f"""你是专业的有机化学家。为分子式 {formula} 生成 {needed} 个新的同分异构体。

已存在的异构体（请避免重复）：
{existing_str}

要求：
1. 生成化学上合理且唯一的异构体
2. 覆盖不同的结构类型（链式、支链、环状、芳香等，如果适用）
3. 只输出SMILES字符串，每行一个
4. 不要输出任何解释

输出格式：
每行一个SMILES字符串，例如：
C1CCCCC1
CCCCCC
...
"""
            
            try:
                response = self.client.chat.completions.create(
                    model=self.model,
                    messages=[
                        {"role": "system", "content": "你是专业的有机化学专家，擅长生成同分异构体。"},
                        {"role": "user", "content": prompt}
                    ],
                    temperature=0.7,
                    max_tokens=1000
                )
                
                candidate_smiles = [line.strip() for line in response.choices[0].message.content.strip().split('\n') if line.strip()]
                
                # 验证和去重
                for smi in candidate_smiles:
                    if len(new_smiles) >= target_count:
                        break
                    
                    normalized = self._normalize_smiles(smi)
                    if not normalized:
                        continue
                    
                    if normalized in existing_normalized or normalized in [self._normalize_smiles(s) for s in new_smiles]:
                        continue
                    
                    # RDKit验证
                    validation = self._validate_smiles_rdkit(smi)
                    if validation["valid"] and validation["formula"] == formula:
                        new_smiles.append(smi)
                        existing_normalized.add(normalized)
                
                logger.info(f"尝试 {attempt}: 获得 {len(new_smiles)} 个新异构体")
                
            except Exception as e:
                logger.error(f"第 {attempt} 次生成失败: {e}")
        
        result = {
            "success": True,
            "formula": formula,
            "target_count": target_count,
            "generated_count": len(new_smiles),
            "new_smiles": new_smiles,
            "attempts": attempt
        }
        
        self._log_operation("ai_supplement", result)
        
        logger.info(f"AI补充完成: 生成了 {len(new_smiles)} 个新异构体")
        return result
    
    def safe_add_molecules(self, smiles_list: List[str], formula: str, 
                          category: str = "ai_generated",
                          add_name: bool = False,
                          source: str = "deepseek") -> Dict[str, Any]:
        """安全地添加分子到数据库
        
        Args:
            smiles_list: 待添加的SMILES列表
            formula: 目标分子式
            category: 分子类别
            add_name: 是否自动生成名称
            source: 来源标识
            
        Returns:
            添加结果字典
        """
        logger.info(f"开始安全添加分子到数据库: {len(smiles_list)} 个分子")
        
        # 创建备份
        backup_path = self._create_backup("add")
        
        # AI验证
        verification = self.ai_verify_molecules(smiles_list, formula)
        
        if not verification["success"]:
            logger.error("AI验证失败，取消添加操作")
            return {
                "success": False,
                "error": "AI验证失败",
                "backup": backup_path
            }
        
        verified_smiles = verification["verified_smiles"]
        added_count = 0
        failed_count = 0
        added_molecules = []
        
        # 确保类别存在
        if category not in self.library_data["categories"]:
            self.library_data["categories"][category] = []
        
        # 添加分子
        for i, smiles in enumerate(verified_smiles):
            try:
                # 生成唯一ID
                molecule_id = f"mol_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{len(self.library_data['categories'][category])}_{i}"
                
                # 获取分子信息
                validation = self._validate_smiles_rdkit(smiles)
                
                # 生成名称
                name = f"{formula}_异构体_{len(self.library_data['categories'][category]) + 1}" if add_name else ""
                
                molecule_record = {
                    "id": molecule_id,
                    "name": name,
                    "smiles": smiles,
                    "formula": validation["formula"],
                    "molecular_weight": validation["molecular_weight"],
                    "carbon_count": validation["carbon_count"],
                    "hydrogen_count": validation["hydrogen_count"],
                    "oxygen_count": validation["oxygen_count"],
                    "description": f"由 {source} AI 生成的{formula}异构体",
                    "category": category,
                    "tags": ["ai_generated", source, formula],
                    "created_date": datetime.now().isoformat(),
                    "custom_properties": {
                        "source": source,
                        "verified_by_ai": True
                    }
                }
                
                self.library_data["categories"][category].append(molecule_record)
                added_count += 1
                added_molecules.append(molecule_record)
                
            except Exception as e:
                logger.error(f"添加分子失败: {e}")
                failed_count += 1
        
        # 保存数据库
        try:
            self.library_data["last_updated"] = datetime.now().isoformat()
            self.library_data["total_count"] = sum(len(v) for v in self.library_data["categories"].values())
            
            with open(self.library_file, 'w', encoding='utf-8') as f:
                json.dump(self.library_data, f, ensure_ascii=False, indent=2)
            
            self._log_operation("safe_add", {
                "input_count": len(smiles_list),
                "verified_count": len(verified_smiles),
                "added_count": added_count,
                "failed_count": failed_count,
                "category": category,
                "source": source
            })
            
            logger.info(f"安全添加完成: 成功 {added_count}, 失败 {failed_count}")
            
            return {
                "success": True,
                "added_count": added_count,
                "failed_count": failed_count,
                "added_molecules": added_molecules,
                "backup": backup_path
            }
            
        except Exception as e:
            logger.error(f"保存数据库失败: {e}")
            logger.info("尝试从备份恢复...")
            
            # 尝试恢复备份
            try:
                with open(backup_path, 'r', encoding='utf-8') as f:
                    self.library_data = json.load(f)
                logger.info("备份恢复成功")
            except:
                logger.error("备份恢复失败")
            
            return {
                "success": False,
                "error": f"保存失败: {str(e)}",
                "backup": backup_path
            }
    
    def ai_enhance_library(self, target_formula: str, target_total: int = None, 
                          current_total: int = None) -> Dict[str, Any]:
        """使用AI增强分子库
        
        Args:
            target_formula: 目标分子式
            target_total: 目标总数量（如果提供，将自动补充到这个数量）
            current_total: 当前已有数量（如果不提供，将从数据库查询）
            
        Returns:
            增强结果字典
        """
        logger.info(f"开始AI增强分子库: {target_formula}")
        
        # 获取当前已有的分子
        if current_total is None:
            existing_mols = []
            for category, molecules in self.library_data["categories"].items():
                for mol in molecules:
                    if mol.get("formula") == target_formula:
                        existing_mols.append(mol["smiles"])
            current_total = len(existing_mols)
        else:
            existing_mols = []
            for category, molecules in self.library_data["categories"].items():
                for mol in molecules:
                    if mol.get("formula") == target_formula:
                        existing_mols.append(mol["smiles"])
        
        logger.info(f"当前分子库中 {target_formula} 有 {current_total} 个分子")
        
        # 计算需要补充的数量
        if target_total is None:
            target_count = max(5, current_total // 2)  # 默认补充50%
        else:
            target_count = max(0, target_total - current_total)
        
        if target_count == 0:
            logger.info("已达到目标数量，无需补充")
            return {
                "success": True,
                "message": "已达到目标数量",
                "current_total": current_total,
                "added_count": 0
            }
        
        # 使用AI补充
        supplement_result = self.ai_supplement_isomers(target_formula, existing_mols, target_count)
        
        if not supplement_result["success"]:
            return supplement_result
        
        new_smiles = supplement_result["new_smiles"]
        
        # 安全添加到数据库
        add_result = self.safe_add_molecules(
            new_smiles, 
            target_formula, 
            category="ai_generated",
            add_name=True
        )
        
        result = {
            "success": add_result["success"],
            "formula": target_formula,
            "original_count": current_total,
            "added_count": add_result["added_count"],
            "new_total": current_total + add_result["added_count"],
            "details": {
                "supplement": supplement_result,
                "add": add_result
            }
        }
        
        self._log_operation("ai_enhance", result)
        
        logger.info(f"AI增强完成: {target_formula} 从 {current_total} 增加到 {result['new_total']}")
        
        return result
    
    def rollback_operation(self, backup_file: str) -> Dict[str, Any]:
        """回滚到指定的备份
        
        Args:
            backup_file: 备份文件名（完整路径或仅文件名）
            
        Returns:
            回滚结果
        """
        if not os.path.isabs(backup_file):
            backup_path = os.path.join(self.backup_dir, backup_file)
        else:
            backup_path = backup_file
        
        if not os.path.exists(backup_path):
            return {
                "success": False,
                "error": f"备份文件不存在: {backup_path}"
            }
        
        logger.info(f"开始回滚到备份: {backup_path}")
        
        try:
            # 创建当前状态的备份
            current_backup = self._create_backup("before_rollback")
            
            # 读取备份
            with open(backup_path, 'r', encoding='utf-8') as f:
                backup_data = json.load(f)
            
            # 写入数据库
            with open(self.library_file, 'w', encoding='utf-8') as f:
                json.dump(backup_data, f, ensure_ascii=False, indent=2)
            
            self.library_data = backup_data
            
            self._log_operation("rollback", {
                "backup_file": backup_file,
                "pre_rollback_backup": current_backup
            })
            
            logger.info("回滚成功")
            
            return {
                "success": True,
                "restored_from": backup_path,
                "new_backup": current_backup
            }
            
        except Exception as e:
            logger.error(f"回滚失败: {e}")
            return {
                "success": False,
                "error": str(e)
            }
    
    def get_backup_list(self) -> List[Dict[str, Any]]:
        """获取所有备份列表
        
        Returns:
            备份信息列表
        """
        backups = []
        
        if not os.path.exists(self.backup_dir):
            return backups
        
        for filename in sorted(os.listdir(self.backup_dir), reverse=True):
            if filename.endswith('.json'):
                filepath = os.path.join(self.backup_dir, filename)
                try:
                    stat = os.stat(filepath)
                    backups.append({
                        "filename": filename,
                        "path": filepath,
                        "size": stat.st_size,
                        "created": datetime.fromtimestamp(stat.st_ctime).isoformat()
                    })
                except:
                    pass
        
        return backups
    
    def get_library_statistics(self) -> Dict[str, Any]:
        """获取分子库统计信息
        
        Returns:
            统计信息字典
        """
        stats = {
            "total_molecules": 0,
            "by_category": {},
            "by_formula": {},
            "by_source": {}
        }
        
        for category, molecules in self.library_data["categories"].items():
            stats["by_category"][category] = len(molecules)
            stats["total_molecules"] += len(molecules)
            
            for mol in molecules:
                formula = mol.get("formula", "unknown")
                stats["by_formula"][formula] = stats["by_formula"].get(formula, 0) + 1
                
                if "custom_properties" in mol and "source" in mol["custom_properties"]:
                    source = mol["custom_properties"]["source"]
                    stats["by_source"][source] = stats["by_source"].get(source, 0) + 1
        
        return stats

# 便捷函数
def quick_ai_update(library_file: str = "molecule_library.json", 
                   formulas: List[str] = None,
                   target_per_formula: int = 20,
                   api_key: str = None) -> Dict[str, Any]:
    """快速AI更新分子库
    
    Args:
        library_file: 分子库文件路径
        formulas: 要更新的分子式列表
        target_per_formula: 每个分子的目标数量
        api_key: API密钥
        
    Returns:
        更新结果字典
    """
    if formulas is None:
        formulas = ["C8H18", "C9H20", "C10H22"]
    
    results = {}
    
    try:
        integration = DeepSeekDatabaseIntegration(api_key=api_key, library_file=library_file)
        
        for formula in formulas:
            result = integration.ai_enhance_library(
                target_formula=formula,
                target_total=target_per_formula
            )
            results[formula] = result
        
        return {
            "success": True,
            "results": results,
            "summary": {
                "total_formulas": len(formulas),
                "total_added": sum(r.get("added_count", 0) for r in results.values() if r.get("success"))
            }
        }
        
    except Exception as e:
        return {
            "success": False,
            "error": str(e)
        }

# 快捷函数 - 兼容 MolGenPlus_Main.py 的调用接口
def ai_deduplicate_smiles(smiles_list: List[str], use_smart_analysis: bool = False) -> Dict[str, Any]:
    """AI智能去重SMILES列表

    Args:
        smiles_list: SMILES字符串列表
        use_smart_analysis: 是否使用智能分析（已废弃，保留兼容性）

    Returns:
        去重结果字典
    """
    try:
        api_key = os.environ.get('DEEPSEEK_API_KEY')
        if not api_key:
            raise ValueError("DEEPSEEK_API_KEY环境变量未设置")

        integration = DeepSeekDatabaseIntegration(api_key=api_key)
        verification = integration.ai_verify_molecules(smiles_list, formula=None)

        if verification["success"]:
            return {
                "success": True,
                "unique_smiles": verification["verified_smiles"],
                "removed_count": len(smiles_list) - len(verification["verified_smiles"]),
                "duplicates": verification.get("duplicates", {})
            }
        else:
            return {
                "success": False,
                "error": verification.get("error", "未知错误"),
                "unique_smiles": smiles_list,  # 返回原始列表作为降级处理
                "removed_count": 0,
                "duplicates": {}
            }
    except Exception as e:
        logger.error(f"AI去重失败: {e}")
        # 降级处理：使用RDKit进行简单去重
        unique_set = set()
        duplicates = {}

        for i, smiles in enumerate(smiles_list):
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
                    if canonical in unique_set:
                        duplicates.setdefault(canonical, []).append(i)
                    else:
                        unique_set.add(canonical)
            except:
                pass

        unique_smiles = list(unique_set)
        return {
            "success": True,
            "unique_smiles": unique_smiles,
            "removed_count": len(smiles_list) - len(unique_smiles),
            "duplicates": duplicates
        }


def ai_enhance_isomer_generation(formula: str, current_smiles: List[str]) -> Dict[str, Any]:
    """AI增强异构体生成

    Args:
        formula: 目标分子式
        current_smiles: 当前的SMILES列表

    Returns:
        增强结果字典
    """
    try:
        api_key = os.environ.get('DEEPSEEK_API_KEY')
        if not api_key:
            raise ValueError("DEEPSEEK_API_KEY环境变量未设置")

        integration = DeepSeekDatabaseIntegration(api_key=api_key, library_file="molecule_library.json")

        # 请求AI生成新的异构体候选
        result = integration.ai_supplement_isomers(
            target_formula=formula,
            target_count=max(20, len(current_smiles) * 2)  # 至少生成20个或翻倍
        )

        if result["success"]:
            # 验证新生成的分子是否与已有分子重复
            verification = integration.ai_verify_molecules(
                result["new_smiles"] + current_smiles,
                formula=formula
            )

            if verification["success"]:
                final_smiles = verification["verified_smiles"]
                new_candidates = [s for s in final_smiles if s not in current_smiles]

                return {
                    "success": True,
                    "original_count": len(current_smiles),
                    "new_candidates": new_candidates,
                    "final_smiles": final_smiles
                }
            else:
                return {
                    "success": False,
                    "error": verification.get("error", "验证失败"),
                    "new_candidates": [],
                    "final_smiles": current_smiles
                }
        else:
            return {
                "success": False,
                "error": result.get("error", "生成失败"),
                "new_candidates": [],
                "final_smiles": current_smiles
            }
    except Exception as e:
        logger.error(f"AI增强失败: {e}")
        return {
            "success": False,
            "error": str(e),
            "new_candidates": [],
            "final_smiles": current_smiles
        }


def quick_dedup(smiles_list: List[str]) -> List[str]:
    """快速去重（简化版）"""
    result = ai_deduplicate_smiles(smiles_list, use_smart_analysis=False)
    return result.get("unique_smiles", smiles_list)


def quick_enhance(formula: str, smiles_list: List[str]) -> List[str]:
    """快速增强（简化版）"""
    result = ai_enhance_isomer_generation(formula, smiles_list)
    return result.get("final_smiles", smiles_list)


if __name__ == "__main__":
    # 示例使用
    print("DeepSeek数据库集成模块")
    print("=" * 60)

    # 示例1: 快速更新分子库
    print("\n示例1: 快速更新C8H18到20个异构体")
    result = quick_ai_update(
        library_file="MolGenPlus_Project/molecule_library.json",
        formulas=["C8H18"],
        target_per_formula=20,
        api_key=os.environ.get('DEEPSEEK_API_KEY')
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))

    # 示例2: 手动使用集成类
    print("\n示例2: 手动补充C10H22异构体")
    integration = DeepSeekDatabaseIntegration(
        api_key=os.environ.get('DEEPSEEK_API_KEY'),
        library_file="MolGenPlus_Project/molecule_library.json"
    )

    result = integration.ai_enhance_library(
        target_formula="C10H22",
        target_total=30
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))
