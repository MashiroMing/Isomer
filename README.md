# Isomer - 有机分子同分异构体分析工具

<div align="center">

![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![RDKit](https://img.shields.io/badge/RDKit-Cheminformatics-orange.svg)
![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20macOS%20%7C%20Linux-lightgrey.svg)

**一个基于 RDKit 的高性能有机分子同分异构体生成与智能分析平台**

</div>

---

## 目录

- [项目简介](#项目简介)
- [功能特点](#功能特点)
- [技术架构](#技术架构)
- [环境配置](#环境配置)
- [快速开始](#快速开始)
- [使用指南](#使用指南)
- [项目结构](#项目结构)
- [核心技术](#核心技术)
- [许可证](#许可证)

---

## 项目简介

**Isomer** 是一款专为有机化学领域设计的同分异构体分析工具，基于 RDKit 化学信息学库开发。该工具能够根据输入的分子式自动生成所有符合化学价键规则的同分异构体结构，并提供多种智能分析和筛选功能。

### 1.1 设计目标

- 提供快速、准确的同分异构体枚举能力
- 支持多种结构筛选和过滤条件
- 集成人工智能技术实现智能去重与补充
- 提供机器学习驱动的分子分析功能
- 支持主流化学文件格式的导入导出

### 1.2 主要应用场景

| 应用领域 | 具体用途 |
|----------|----------|
| 化学教学 | 有机化学实验设计与辅助教学 |
| 药物研发 | 药物先导化合物的结构探索与筛选 |
| 化学研究 | 分子结构表征与构效关系分析 |
| 竞赛辅助 | 化学奥赛题目验证与答案生成 |

---

## 功能特点

### 2.1 核心功能模块

| 功能模块 | 描述 |
|----------|------|
| **异构体生成** | 支持 C1-C10 碳链的烷烃、烯烃、环烷烃及含氧衍生物的自动生成 |
| **结构筛选** | 支持按环结构、键类型、芳香性、官能团等条件精细筛选 |
| **分子库管理** | 提供分子收藏、分类存储、快速检索功能 |
| **AI 智能增强** | 集成 DeepSeek AI API，实现 SMILES 智能去重与异构体补充 |
| **机器学习分析** | 支持分子聚类、PCA 降维可视化、异常检测等分析 |
| **多格式导出** | 支持 TXT、XYZ、GJF（Gaussian）、CSV 等格式导出 |

### 2.2 结构筛选选项

#### 环结构筛选
- 无环（链状结构）
- 单环
- 双环
- 三环
- 多环

#### 键类型筛选
- 单键（饱和键）
- 双键（烯烃）
- 三键（炔烃）
- 混合键（含多种键型）

#### 芳香性筛选
- 无芳香（非芳香族）
- 含苯环（苯及其衍生物）
- 其他芳香环

#### 官能团筛选
- 羟基 (-OH)
- 羰基 (C=O)
- 羧基 (-COOH)
- 酯基 (-COO-)

### 2.3 支持的分子式范围

| 碳数范围 | 支持的分子式示例 | 算法类型 |
|:--------:|------------------|----------|
| C1-C5 | CH₄, C₂H₆, C₃H₈, C₄H₁₀, C₅H₁₂ | 专用优化算法 |
| C6 | C₆H₆ (苯), C₆H₁₂, C₆H₆O | 专用优化算法 |
| C7-C10 | C₇H₁₆, C₈H₁₈, C₉H₂₀, C₁₀H₂₂ | 通用生成算法 |
| C11+ | 任意烷烃分子式 | 通用生成算法 |

---

## 技术架构

### 3.1 技术栈

```
┌─────────────────────────────────────────────────────────────┐
│                        用户界面层                             │
│                    (Tkinter GUI)                            │
├─────────────────────────────────────────────────────────────┤
│                       业务逻辑层                             │
│  ┌─────────────┬──────────────┬──────────────┐              │
│  │  异构体生成  │   结构筛选   │  分子库管理  │              │
│  │  generators │  analysis    │  database    │              │
│  └─────────────┴──────────────┴──────────────┘              │
├─────────────────────────────────────────────────────────────┤
│                        数据层                                │
│     RDKit (化学信息学)  │  JSON (数据持久化)                 │
├─────────────────────────────────────────────────────────────┤
│                       外部服务                                │
│         OpenAI/DeepSeek API (AI 增强)                        │
└─────────────────────────────────────────────────────────────┘
```

### 3.2 核心依赖

| 依赖包 | 版本要求 | 用途说明 |
|--------|----------|----------|
| Python | ≥ 3.8 | 运行环境 |
| RDKit | 最新版 | 化学结构处理与分子指纹计算 |
| NumPy | ≥ 1.20 | 数值计算 |
| Pandas | ≥ 1.3 | 数据处理与分析 |
| scikit-learn | ≥ 0.24 | 机器学习分析 |
| Pillow | ≥ 8.0 | 分子结构图像渲染 |
| OpenAI | ≥ 1.0 | AI 功能集成（可选） |

---

## 环境配置

### 4.1 系统要求

- **操作系统**: Windows 10/11, macOS 10.14+, Linux (Ubuntu 20.04+)
- **内存**: 建议 4GB 及以上
- **磁盘**: 建议 500MB 可用空间

### 4.2 安装步骤

#### 步骤一：安装 Python 环境

建议使用 Anaconda 或 Miniconda 管理 Python 环境：

```bash
# 使用 conda 创建虚拟环境
conda create -n isomer python=3.10
conda activate isomer
```

#### 步骤二：安装 RDKit（关键依赖）

RDKit 推荐通过 conda 安装：

```bash
conda install -c conda-forge rdkit
```

#### 步骤三：安装其他依赖

```bash
# 核心依赖
pip install numpy pandas pillow

# 机器学习依赖（可选）
pip install scikit-learn matplotlib

# AI 功能依赖（可选）
pip install openai
```

#### 步骤四：验证安装

```bash
python -c "from rdkit import Chem; print('RDKit 安装成功')"
```

---

## 快速开始

### 5.1 下载与安装

```bash
# 克隆项目仓库
git clone https://github.com/MashiroMing/Isomer.git
cd Isomer

# 安装依赖
conda install -c conda-forge rdkit
pip install openai pillow scikit-learn numpy pandas matplotlib
```

### 5.2 运行程序

**方式一：命令行运行**

```bash
python MolGenPlus_Main.py
```

**方式二：使用启动脚本（Windows）**

```bash
# 双击运行
run_molgen.bat

# 或在命令行中运行
cmd /k run_molgen.bat
```

### 5.3 首次使用

1. 运行程序后，主界面将显示分子式输入框
2. 输入有效的分子式（如 `C4H10`）
3. 选择筛选条件（可选）
4. 点击"生成"按钮获取异构体列表
5. 选择感兴趣的异构体进行查看或导出

---

## 使用指南

### 6.1 基本操作流程

```
输入分子式 → 选择筛选条件 → 生成异构体 → 浏览结果 → 导出保存
```

### 6.2 分子式输入规范

- 使用标准分子式格式（如 `C4H10`、`C6H6`）
- 支持含氧有机物（如 `C6H6O` 表示苯酚）
- 不区分大小写

### 6.3 导出格式说明

| 格式 | 文件扩展名 | 适用场景 |
|------|------------|----------|
| TXT | `.txt` | 纯文本记录，通用性强 |
| XYZ | `.xyz` | 通用化学文件格式，兼容多数可视化工具 |
| GJF | `.gjf` | Gaussian 软件输入文件，用于量化计算 |
| CSV | `.csv` | 电子表格格式，便于数据分析 |

---

## 项目结构

```
Isomer/
├── MolGenPlus_Main.py          # 主程序入口（GUI 应用）
├── run_molgen.bat             # Windows 启动脚本
│
├── analysis/                  # 分析模块
│   ├── alkane_isomer_counter.py    # 烷烃异构体计数
│   ├── ml_isomer_analyzer.py        # 机器学习分析器
│   └── structure_filter.py          # 结构过滤器
│
├── core/                      # 核心算法
│
├── database/                  # 数据库模块
│   ├── molecule_library.py         # 分子库管理
│   ├── predefined_molecules.py      # 预定义分子
│   └── deepseek_database_integration.py  # AI 集成
│
├── generators/                # 异构体生成器
│   ├── isomer_c1_c5_combined.py    # C1-C5 异构体生成
│   ├── isomer_c6_combined.py       # C6 异构体生成
│   ├── isomer_c8_c10_combined.py   # C8-C10 异构体生成
│   └── isomer_general.py           # 通用异构体生成
│
├── utils/                    # 工具函数
│   ├── check_dependencies.py      # 依赖检查
│   ├── mol_utils.py               # 分子工具
│   └── output_optimizer.py         # 输出优化
│
├── data/                     # 数据目录
│   ├── molecule_library.json      # 分子库数据
│   └── database_backups/          # 数据备份
│
├── docs/                     # 文档目录
│
└── tests/                    # 测试目录
```

---

## 核心技术

### 8.1 Morgan 分子指纹

本项目使用 **Morgan 指纹**（ECFP 的 RDKit 实现）进行分子结构编码：

```python
# 半径=2 的 Morgan 指纹，1024位
fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
```

**特点**：
- 结构相似的分子具有相似的指纹向量
- 支持分子聚类和相似性搜索
- 计算效率高，适合大规模分子分析

### 8.2 机器学习分析

基于提取的分子特征（分子量、Wiener 指数、指纹等），支持：

- **K-Means 聚类**：将异构体按结构特征分组
- **PCA 降维**：高维特征可视化
- **异常检测**：识别结构异常的分子

### 8.3 AI 增强功能

集成 DeepSeek API 实现智能功能：

- **SMILES 去重**：智能识别化学等价结构
- **异构体补充**：AI 分析可能遗漏的结构

---

## 许可证

本项目基于 [MIT License](LICENSE) 开源。

---

## 致谢

- [RDKit](https://www.rdkit.org/) - 开源化学信息学工具包
- [scikit-learn](https://scikit-learn.org/) - 机器学习库
- [DeepSeek](https://www.deepseek.com/) - 大语言模型 API
