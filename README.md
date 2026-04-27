# Isomer - 有机分子同分异构体分析工具

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![RDKit](https://img.shields.io/badge/RDKit-Cheminformatics-orange.svg)

一个强大的有机分子同分异构体分析工具，支持 AI 智能增强和机器学习分析。

## 功能特点

- 🔬 **同分异构体生成** - 输入分子式，自动生成所有可能的结构
- 🧪 **结构筛选** - 按环结构、键类型、芳香性、官能团筛选
- 💾 **分子库管理** - 保存和管理常用分子
- 🤖 **AI 智能增强** - 使用 DeepSeek AI 智能去重和补充异构体
- 📊 **机器学习分析** - 聚类分析、PCA 降维、异常检测
- 📁 **多格式导出** - 支持 TXT、XYZ、GJF、CSV 格式

## 支持的分子式

| 碳数 | 示例分子式 |
|:----:|------------|
| C1-C5 | CH₄, C₂H₆, C₃H₈, C₄H₁₀, C₅H₁₂ |
| C6 | C₆H₆ (苯), C₆H₁₂, C₆H₆O |
| C7-C10 | C₇H₁₆, C₈H₁₈, C₉H₂₀, C₁₀H₂₂ |
| C11+ | 通用算法生成 |

## 快速开始

### 环境要求

- Python 3.8+
- [RDKit](https://www.rdkit.org/docs/Install.html) (通过 conda 安装)

### 安装

```bash
# 1. 克隆项目
git clone https://github.com/MashiroMing/Isomer.git
cd Isomer

# 2. 安装 RDKit (推荐使用 conda)
conda install -c conda-forge rdkit

# 3. 安装其他依赖
pip install openai pillow scikit-learn numpy pandas matplotlib
```

### 运行

```bash
# Windows
python MolGenPlus_Main.py

# 或双击运行 run_molgen.bat
```

## 项目结构

```
Isomer/
├── analysis/           # 分析模块 (ML分析、结构过滤)
├── core/              # 核心算法
├── database/          # 数据库和 AI 集成
├── docs/              # 文档
├── generators/        # 同分异构体生成器
├── utils/             # 工具函数
├── tests/             # 测试文件
├── MolGenPlus_Main.py # 主程序入口
└── run_molgen.bat     # Windows 启动脚本
```

## 技术栈

- **化学信息学**: RDKit
- **AI**: OpenAI / DeepSeek API
- **机器学习**: scikit-learn
- **GUI**: Tkinter
- **数据处理**: NumPy, Pandas

## 许可证

本项目基于 MIT 许可证开源。

## 贡献

欢迎提交 Issue 和 Pull Request！
