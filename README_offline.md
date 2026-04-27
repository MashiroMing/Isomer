# MolGenPlus 有机分子同分异构体分析工具

> 简体中文 | [English](./docs/README_EN.md)

---

## 目录

1. [项目简介](#1-项目简介)
2. [功能特点](#2-功能特点)
3. [环境要求](#3-环境要求)
4. [快速开始](#4-快速开始)
5. [使用教程](#5-使用教程)
6. [核心代码示例](#6-核心代码示例)
7. [常见问题](#7-常见问题)
8. [技术架构](#8-技术架构)

---

## 1. 项目简介

### 1.1 什么是 MolGenPlus？

**MolGenPlus** 是一个**有机分子同分异构体分析工具**，可以帮助你：

- 🔬 **生成同分异构体**：输入分子式，自动生成所有可能的结构
- 🧪 **结构筛选**：按环结构、键类型、芳香性、官能团筛选
- 💾 **分子库管理**：保存和管理常用分子
- 🤖 **AI 增强**：使用 DeepSeek AI 智能去重和补充异构体
- 📊 **机器学习分析**：聚类分析、PCA降维、异常检测

### 1.2 支持的分子式

| 碳数 | 支持的分子式示例 |
|:----:|------------------|
| C1-C5 | CH₄, C₂H₆, C₃H₈, C₄H₁₀, C₅H₁₂ |
| C6 | C₆H₆ (苯), C₆H₁₂, C₆H₆O |
| C7-C10 | C₇H₁₆, C₈H₁₈, C₉H₂₀, C₁₀H₂₂ |
| C11+ | 通用算法生成 |

### 1.3 典型应用场景

- 有机化学教学和实验设计
- 化学研究人员探索分子结构
- 药物化学前期分子筛选
- 化学竞赛辅助工具

---

## 2. 功能特点

### 2.1 核心功能

```

 MolGenPlus 功能矩阵                       

│   功能模块     │

异构体生成     
结构筛选          
分子库管理     
AI 智能去重    
AI 补充异构体  
机器学习分析   

```

### 2.2 结构筛选选项

| 筛选类型 | 可选条件 |
|---------|---------|
| 环结构 | 无环、单环、双环、三环、多环 |
| 键类型 | 单键、双键、三键、混合键 |
| 芳香性 | 无芳香、含苯环、其他芳香环 |
| 官能团 | 羟基(-OH)、羰基(C=O)、羧基(-COOH)、酯基(-COO-) |

### 2.3 导出格式

- **TXT** - 纯文本格式
- **XYZ** - 通用化学文件格式（兼容 GaussView）
- **GJF** - Gaussian 软件输入文件
- **CSV** - Excel 兼容格式

---

## 3. 环境要求

### 3.1 操作系统

- ✅ Windows 10/11 (推荐)
- ✅ macOS 10.14+
- ✅ Linux (Ubuntu 20.04+)

### 3.2 Python 版本

- Python 3.8 或更高版本

### 3.3 必需依赖

| 包名 | 用途 | 安装命令 |
|------|------|---------|
| rdkit | 化学结构处理 | `conda install -c conda-forge rdkit` |
| pillow | 图像处理 | `pip install pillow` |
| numpy | 数值计算 | `pip install numpy` |

### 3.4 可选依赖

| 包名 | 用途 | 安装命令 |
|------|------|---------|
| openai | DeepSeek AI 功能 | `pip install openai` |
| scikit-learn | 机器学习分析 | `pip install scikit-learn` |
| pandas | 数据处理 | `pip install pandas` |
| matplotlib | 数据可视化 | `pip install matplotlib` |

---

## 4. 快速开始

### 4.1 安装步骤

#### 第一步：下载项目

```bash
# 克隆或下载项目到本地
git clone <项目地址>
cd MolGenPlus_Project
```

#### 第二步：安装依赖

**方法 A：使用安装脚本（推荐）**

```bash
# Windows
cd docs
install_dependencies.bat

# 或直接双击运行 install_dependencies.bat
```

**方法 B：使用 Conda（推荐安装 rdkit）**

```bash
# 1. 安装 rdkit（推荐）
conda install -c conda-forge rdkit

# 2. 安装其他依赖
pip install openai pillow scikit-learn numpy pandas matplotlib
```

**方法 C：使用 pip**

```bash
pip install rdkit openai pillow scikit-learn numpy pandas matplotlib
```

#### 第三步：验证安装

```bash
cd utils
python check_dependencies.py
```

如果看到以下输出，说明安装成功：
```
======================================
MolGenPlus 项目依赖检查
======================================

✓ rdkit            - 已安装 (用于化学结构处理和SMILES验证)
✓ openai           - 已安装 (用于DeepSeek AI API集成)
✓ PIL              - 已安装 (Pillow - 用于图像处理)
✓ sklearn          - 已安装 (scikit-learn - 用于机器学习分析)
✓ numpy            - 已安装 (用于数值计算)
✓ pandas           - 已安装 (用于数据处理)

======================================
✓ 所有依赖都已安装，项目可以正常运行！
======================================
```

#### 第四步：运行程序

```bash
# 方法 1：使用启动脚本
run_molgen.bat

# 方法 2：直接运行 Python
python MolGenPlus_Main.py
```

### 4.2 首次运行设置

1. **启动程序**：双击 `run_molgen.bat`
2. **自动初始化**：程序会自动检查分子库，如果为空则导入50个预定义分子
3. **开始使用**：在界面中输入分子式，点击"生成"按钮

---

## 5. 使用教程

### 5.1 主界面介绍

```

```

### 5.2 基础使用流程

#### 步骤 1：输入分子式

在"分子式输入"框中输入化学式，例如：
- `C8H18` - 辛烷
- `C6H6` - 苯
- `C3H6O` - 丙酮

#### 步骤 2：选择结构类型

勾选需要的结构类型：
- ☑ 链烷烃
- ☑ 环烷烃
- ☑ 烯烃
- ☑ 炔烃

#### 步骤 3：生成异构体

点击"生成"按钮，程序会自动生成所有符合条件 的同分异构体。

#### 步骤 4：查看和筛选结果

生成的结果显示在结果区域，你可以通过：
- 点击单条结果查看详细信息
- 使用"结构筛选"进一步过滤
- 使用"AI去重"去除重复结构

#### 步骤 5：导出数据

点击"导出"按钮，选择格式：
- TXT - 纯文本
- XYZ - 可视化软件格式
- GJF - Gaussian 计算输入

### 5.3 高级功能

#### AI 增强（需要 DeepSeek API）

1. 点击"AI增强"按钮
2. 在弹窗中输入 DeepSeek API 密钥
3. 选择增强选项：
   - 补充缺失异构体
   - 生成新候选结构
   - 智能去重

#### 机器学习分析（需要 scikit-learn）

1. 确保已安装 scikit-learn
2. 点击"ML分析"按钮
3. 选择分析类型：
   - 聚类分析
   - PCA 降维可视化
   - 异常检测
   - 物性预测

---

## 6. 核心代码示例

> 这部分面向需要二次开发或学习代码的用户

### 6.1 基本使用 - 生成异构体

```python
# 导入必要的模块
import sys
import os

# 添加项目路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 导入异构体生成器
from generators.isomer_c1_c5_combined import generate_isomers_c1_c5
from generators.isomer_c6_combined import generate_isomers_c6
from generators.isomer_c8_c10_combined import generate_isomers_c8_c10
from generators.isomer_general import generate_isomers_general

# 生成 C8H18 的所有异构体
formula = "C8H18"
isomers = generate_isomers_c8_c10(formula)

print(f"分子式 {formula} 的同分异构体:")
print(f"共生成 {len(isomers)} 个结构")

# 打印前5个
for i, isomer in enumerate(isomers[:5], 1):
    print(f"{i}. {isomer}")
```

**输出示例：**
```
分子式 C8H18 的同分异构体:
共生成 18 个结构
1. CCCCCCCC
2. CCC(C)CCC
3. CC(C)(C)CCC
4. CCCCC(C)C
5. CCC(C)(C)CC
```

### 6.2 使用结构筛选

```python
from analysis.structure_filter import StructureFilter

# 假设这是我们生成的异构体列表
smiles_list = [
    "CCCCCCCC",      # 正辛烷
    "CCC(C)CCC",     # 2-甲基庚烷
    "CC(C)(C)CCC",   # 2,2-二甲基己烷
    "C1CCCCCC1",     # 环辛烷
    "C=CCCCCCC",     # 辛烯
]

# 创建筛选器
sf = StructureFilter()

# 筛选：无环、仅单键的烷烃
filtered = sf.filter_by_structure_type(smiles_list, allow_rings=False)
print(f"无环结构: {filtered}")

# 筛选：仅包含单键的烷烃
alkanes_only = sf.filter_by_bond_type(smiles_list, single=True, double=False, triple=False)
print(f"仅烷烃: {alkanes_only}")
```

### 6.3 使用 DeepSeek AI 进行智能去重

```python
import os
from database.deepseek_database_integration import ai_deduplicate_smiles, ai_enhance_isomer_generation

# 设置 API 密钥（方式1：环境变量）
os.environ['DEEPSEEK_API_KEY'] = 'your-api-key-here'

# 或者方式2：直接传入
api_key = 'your-api-key-here'

# 示例 SMILES 列表（有重复）
smiles_with_duplicates = [
    "CCCCCCCC",      # 正辛烷
    "CCCCCCCC",      # 重复！
    "CCC(C)CCC",     # 2-甲基庚烷
    "c1ccccc1",      # 苯
    "C1=CC=CC=C1",   # 苯（另一种表示）
]

# AI 智能去重
result = ai_deduplicate_smiles(smiles_with_duplicates, use_smart_analysis=False)

print(f"原始数量: {len(smiles_with_duplicates)}")
print(f"去重后: {len(result['unique_smiles'])}")
print(f"移除重复: {result['removed_count']}")
print(f"唯一结构: {result['unique_smiles']}")
```

### 6.4 使用 AI 增强异构体生成

```python
import os
from database.deepseek_database_integration import ai_enhance_isomer_generation

# 设置 API 密钥
os.environ['DEEPSEEK_API_KEY'] = 'your-api-key-here'

# 当前已有的异构体
current_isomers = ["CCCCCCCC", "CCC(C)CCC", "CC(C)(C)CCC"]

# 使用 AI 增强 - 生成更多候选
formula = "C8H18"
result = ai_enhance_isomer_generation(formula, current_isomers)

if result['success']:
    print(f"原始数量: {result['original_count']}")
    print(f"AI新增候选: {len(result['new_candidates'])}")
    print(f"最终结构数: {len(result['final_smiles'])}")

    if result['new_candidates']:
        print("\n新生成的候选结构:")
        for smiles in result['new_candidates']:
            print(f"  + {smiles}")
else:
    print(f"增强失败: {result.get('error')}")
```

### 6.5 分子库管理

```python
from database.molecule_library import molecule_library

# 添加分子到库中
result = molecule_library.add_molecule(
    smiles="CCCC",
    name="正丁烷",
    description="C4H10 直链烷烃",
    category="alkanes",
    tags=["烷烃", "直链"]
)

if result['success']:
    print("分子添加成功！")
else:
    print(f"添加失败: {result['error']}")

# 查询分子
butane = molecule_library.get_molecule_by_name("正丁烷")
if butane:
    print(f"找到分子: {butane['smiles']}")

# 获取统计信息
stats = molecule_library.get_statistics()
print(f"分子库总数: {stats['total_molecules']}")
print(f"按类别: {stats['by_category']}")
```

### 6.6 机器学习分析

```python
from analysis.ml_isomer_analyzer import analyze_isomers

# 待分析的 SMILES 列表
smiles_list = [
    "CCCCCCCC", "CCC(C)CCC", "CC(C)(C)CCC",
    "C1CCCCCC1", "C=CCCCCCC", "C#CCCCCC"
]

# 执行完整分析
report = analyze_isomers(smiles_list, analysis_type='all')

print("=" * 50)
print("分子分析报告")
print("=" * 50)

# 分子描述符
print("\n【分子描述符】")
for desc, value in report['descriptors'].items():
    print(f"  {desc}: {value:.4f}")

# 聚类结果
print("\n【聚类分析】")
for cluster_id, members in report['clustering'].items():
    print(f"  簇 {cluster_id}: {len(members)} 个分子")

# 物性预测
print("\n【物性预测】")
for prop, value in report['properties'].items():
    print(f"  {prop}: {value:.2f}")
```

### 6.7 完整示例：批量生成并导出

```python
import os
import json
from generators.isomer_c8_c10_combined import generate_isomers_c8_c10
from analysis.structure_filter import StructureFilter

def batch_generate_and_export():
    """批量生成并导出示例"""

    # 要生成的分子式列表
    formulas = ["C8H18", "C9H20", "C10H22"]

    all_results = {}

    for formula in formulas:
        print(f"\n正在生成 {formula}...")

        # 生成异构体
        isomers = generate_isomers_c8_c10(formula)

        # 结构筛选 - 只保留烷烃（无环）
        sf = StructureFilter()
        alkanes = sf.filter_by_structure_type(isomers, allow_rings=False)

        all_results[formula] = {
            'total': len(isomers),
            'alkanes': len(alkanes),
            'smiles': alkanes[:10]  # 保存前10个
        }

        print(f"  总数: {len(isomers)}, 烷烃: {len(alkanes)}")

    # 保存到 JSON 文件
    output_file = "batch_results.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(all_results, f, ensure_ascii=False, indent=2)

    print(f"\n结果已保存到: {output_file}")
    return all_results

# 运行
if __name__ == "__main__":
    results = batch_generate_and_export()
```

---

## 7. 常见问题

### 7.1 安装问题

#### Q: rdkit 安装失败

**原因**：rdkit 在 Windows 上使用 pip 安装可能出现问题

**解决方案**：使用 conda 安装
```bash
conda install -c conda-forge rdkit
```

#### Q: 提示 "ModuleNotFoundError"

**原因**：依赖未安装

**解决方案**：
```bash
pip install -r requirements.txt
# 或
python docs/install_dependencies.bat
```

### 7.2 使用问题

#### Q: 生成结果为空

**原因**：分子式格式不正确或不支持

**解决方案**：
- 检查分子式格式，如 `C8H18`（不是 `C8 H18`）
- 确认分子式是有效的化学式

#### Q: AI 功能不可用

**原因**：未安装 openai 或未设置 API 密钥

**解决方案**：
```bash
pip install openai
# 然后设置 API 密钥
set DEEPSEEK_API_KEY=your-key
# 或在程序中输入
```

#### Q: 导出文件无法用 GaussView 打开

**原因**：XYZ/GJF 格式问题

**解决方案**：尝试其他导出格式（TXT/CSV），或检查分子坐标是否有效

### 7.3 性能问题

#### Q: 生成大分子时很慢

**原因**：C10+ 异构体数量巨大，计算复杂

**解决方案**：
- 使用结构筛选减少输出
- 分批生成
- 升级硬件（更多 CPU 核心）

---

## 8. 技术架构

### 8.1 项目结构

```
MolGenPlus_Project/
├── MolGenPlus_Main.py          # 🎯 主程序入口（GUI）
├── run_molgen.bat             # Windows 启动脚本
│
├── generators/                 # 🧬 异构体生成模块
│   ├── isomer_c1_c5_combined.py   # C1-C5 生成
│   ├── isomer_c6_combined.py       # C6 生成
│   ├── isomer_c8_c10_combined.py  # C8-C10 生成
│   └── isomer_general.py          # C7+ 通用生成
│
├── analysis/                  # 📊 分析模块
│   ├── structure_filter.py        # 结构筛选
│   ├── ml_isomer_analyzer.py      # 机器学习分析
│   └── alkane_isomer_counter.py   # 烷烃统计
│
├── database/                  # 💾 数据库模块
│   ├── molecule_library.py        # 分子库管理
│   ├── deepseek_database_integration.py  # DeepSeek AI
│   ├── deepseek_database_gui.py   # 数据库 GUI
│   └── predefined_molecules.py    # 预定义分子
│
├── utils/                     # 🛠️ 工具模块
│   ├── mol_utils.py               # 分子工具
│   ├── output_optimizer.py        # 输出优化
│   └── check_dependencies.py      # 依赖检查
│
├── data/                      # 📁 数据目录
│   ├── molecule_library.json      # 分子库数据
│   ├── database_backups/         # 备份
│   └── ml_results/               # ML 结果
│
└── docs/                      # 📚 文档
    ├── README.md                  # 本文档
    └── install_dependencies.bat   # 安装脚本
```

### 8.2 核心技术栈

| 技术 | 用途 | 说明 |
|------|------|------|
| Python 3.8+ | 编程语言 | 核心开发语言 |
| RDKit | 化学计算 | 分子结构处理、SMILES 操作 |
| tkinter | GUI 界面 | 图形用户界面 |
| DeepSeek API | AI 增强 | 智能去重、结构补充 |
| scikit-learn | 机器学习 | 聚类、PCA、预测 |
| NumPy/Pandas | 数据处理 | 数值计算、数据管理 |

### 8.3 工作流程图

```
用户输入分子式
       │
       ▼
┌──────────────────┐
│  解析分子式      │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│  选择生成器      │──► C1-C5 / C6 / C8-C10 / C7+
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│  生成异构体      │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│  RDKit 验证      │──► SMILES 规范化 & 去重
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│  结构筛选（可选）│
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│  AI 增强（可选） │──► DeepSeek API 智能处理
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│  结果展示 & 导出 │
└──────────────────┘
```

---

## 附录

### A. 常用 SMILES 示例

| 分子 | SMILES | 说明 |
|------|--------|------|
| 甲烷 | C | 最简单的烷烃 |
| 乙烷 | CC | 两个碳 |
| 丙烷 | CCC | 三个碳 |
| 正丁烷 | CCCC | 直链丁烷 |
| 异丁烷 | CC(C)C | 2-甲基丙烷 |
| 苯 | c1ccccc1 | 芳香烃 |
| 环己烷 | C1CCCCC1 | 环烷烃 |

### B. 联系方式

- 项目作者：[作者名称]
- 邮箱：[邮箱地址]
- 问题反馈：[GitHub Issues]

---

**祝您使用愉快！** 🎉

> 最后更新：2025年
