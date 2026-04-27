#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_dependencies.py - 依赖检查工具
检查项目运行所需的Python包是否已安装
"""

import sys
import importlib

# 需要检查的依赖包
DEPENDENCIES = {
    'rdkit': '用于化学结构处理和SMILES验证',
    'openai': '用于DeepSeek AI API集成',
    'PIL': 'Pillow - 用于图像处理',
    'sklearn': 'scikit-learn - 用于机器学习分析',
    'numpy': '用于数值计算',
    'pandas': '用于数据处理',
}

def check_dependency(package_name, description):
    """检查单个依赖包"""
    try:
        if package_name == 'PIL':
            importlib.import_module('PIL')
        else:
            importlib.import_module(package_name)
        print(f"✓ {package_name:15s} - 已安装 ({description})")
        return True
    except ImportError:
        print(f"✗ {package_name:15s} - 未安装 ({description})")
        return False

def check_dependencies():
    """检查所有依赖包"""
    print("=" * 70)
    print("MolGenPlus 项目依赖检查")
    print("=" * 70)
    print()

    missing_packages = []
    all_installed = True

    for package_name, description in DEPENDENCIES.items():
        if not check_dependency(package_name, description):
            missing_packages.append(package_name)
            all_installed = False

    print()
    print("=" * 70)

    if all_installed:
        print("✓ 所有依赖都已安装，项目可以正常运行！")
        print("=" * 70)
        return True
    else:
        print("✗ 部分依赖缺失，请安装以下包：")
        print()
        print("安装方法1（推荐conda安装rdkit）：")
        print("  conda install -c conda-forge rdkit")
        print("  pip install openai pillow scikit-learn numpy pandas")
        print()
        print("安装方法2（纯pip安装）：")
        print(f"  pip install {' '.join(missing_packages)}")
        print()
        print("或运行安装脚本：")
        print("  install_dependencies.bat")
        print("=" * 70)
        return False

def main():
    """主检查函数（用于命令行）"""
    result = check_dependencies()
    return 0 if result else 1

if __name__ == "__main__":
    sys.exit(main())
