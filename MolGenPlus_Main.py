# MolGenPlus_Main.py

import tkinter as tk
from tkinter import messagebox, scrolledtext, filedialog
import sys
import re
import os
from typing import Any
from datetime import datetime
import warnings

# 添加项目根目录到 Python 路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 导入输出优化器
try:
    from utils.output_optimizer import setup_silent_mode, get_clean_summary
    setup_silent_mode()
except ImportError:
    # 如果优化器不可用，使用基本静默设置
    warnings.filterwarnings('ignore')
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
    except:
        pass
    os.environ['OMP_NUM_THREADS'] = '1'

# 动态导入项目文件 (使用合并后的模块)
try:
    from generators.isomer_c1_c5_combined import generate_isomers_c1_c5, generate_isomers_c1_c5_hetero, generate_isomers_c1_to_c5
    from generators.isomer_c6_combined import generate_isomers_c6, generate_c6h6_general, generate_isomers_c6_hetero
    from generators.isomer_c8_c10_combined import generate_isomers_c8_c10
    from generators.isomer_general import generate_isomers_general
    from analysis.structure_filter import StructureFilter
    from database.molecule_library import molecule_library
    try:
        from database.molecule_library_gui_legacy import MoleculeLibraryGUI
        HAS_LEGACY_GUI = True
    except ImportError:
        HAS_LEGACY_GUI = False
        print("警告: 传统分子库GUI不可用", file=sys.stderr)

    # 安全导入机器学习模块
    try:
        from analysis.ml_isomer_analyzer import MLIsomerAnalyzer, analyze_isomers
        HAS_ML_ANALYZER = True
    except ImportError:
        HAS_ML_ANALYZER = False  # pyright: ignore[reportConstantRedefinition]
        print("提示: 机器学习分析模块不可用，请安装 scikit-learn: pip install scikit-learn", file=sys.stderr)

    # 导入 DeepSeek AI 模块
    try:
        from database.deepseek_database_integration import (
            DeepSeekDatabaseIntegration,
            quick_ai_update,
            ai_deduplicate_smiles,
            ai_enhance_isomer_generation,
            quick_dedup,
            quick_enhance
        )
        HAS_DEEPSEEK = True
    except ImportError as e:
        HAS_DEEPSEEK = False
        print(f"提示: DeepSeek AI 模块不可用，请安装 openai: pip install openai", file=sys.stderr)
        print(f"详细错误: {e}", file=sys.stderr)

    # 导入烷烃统计模块
    try:
        from analysis.alkane_isomer_counter import get_alkane_isomer_statistics
        HAS_ALKANE_COUNTER = True
    except ImportError:
        HAS_ALKANE_COUNTER = False

except ImportError as e:
    messagebox.showerror("导入错误", f"无法找到核心模块。请确保项目文件结构正确。\n详细错误: {e}")
    sys.exit(1)

# 初始化预定义分子
try:
    from database.predefined_molecules import check_and_initialize_predefined
    init_result = check_and_initialize_predefined()
    if init_result["failed_count"] > 0:
        print(f"警告: {init_result['failed_count']} 个预定义分子初始化失败")
except ImportError as e:
    print(f"警告: 无法加载预定义分子模块: {e}")
except Exception as e:
    print(f"警告: 初始化预定义分子时出错: {e}")


# 导入绘图相关的库
try:
    from rdkit.Chem.Draw import MolToImage
    from PIL import ImageTk
    HAS_DRAWING = True
except ImportError:
    HAS_DRAWING = False
    print("警告: 缺少 Pillow (PIL) 库，无法绘制分子图形。请运行 'pip install Pillow'。", file=sys.stderr)


# ==============================================================================
# 分子式格式化工具函数
# ==============================================================================

def format_formula(formula: str) -> str:
    """
    格式化分子式，符合化学规范（碳数为1时不显示数字）

    Args:
        formula: 原始分子式，如 C1H4, C2H6, C3H6O2

    Returns:
        格式化后的分子式，如 CH4, C2H6, C3H6O2
    """
    formula = formula.strip().upper()

    # 替换 C1 为 C（碳数为1时省略数字）
    # 使用正则确保只替换独立的 C1，而不是 C10, C11 等
    formula = re.sub(r'\bC1\b', 'C', formula)

    return formula


# ==============================================================================
# 核心调度逻辑 (分发到纯烃或含氧模块)
# ==============================================================================

def parse_formula_and_dispatch(formula: str, structure_prefs: dict[str, str] = None) -> dict[str, str] | list[dict[str, Any]]:
    """
    解析分子式 (CnHmOp)，根据碳数和是否含氧调用对应的生成模块。
    现在支持C7及以上的通用算法和结构筛选。
    支持 CH4, C, C2H6, C3H6O2 等多种格式。

    Args:
        formula: 分子式，如 CH4, C, C2H6O, C3H6O2
        structure_prefs: 结构偏好字典，可选
    """
    formula = formula.strip().upper().replace(' ', '')

    # 规范化分子式：处理 C, CH4 等格式为 CnHm 格式
    # 情况1: 只有 C (如甲烷) -> 转换为 C1H4
    if formula == 'C':
        formula = 'C1H4'
    # 情况2: C 后直接跟 H 或 O (如 CH4, CO) -> 插入数字1
    elif re.match(r'^C[H,O]', formula):
        formula = 'C1' + formula[1:]
    # 情况3: 其他格式 (如 C2H6, C3H6O2) -> 不需要修改

    # 检查格式：C开头，后接数字，可选H接数字，可选O接数字
    if not re.match(r'^C\d+(H\d*)?(O\d*)?$', formula):
        return {"error": "分子式格式错误。请使用 CH4, C, CnHmOp 格式 (例如: CH4, C2H6O, C3H6O2)。"}

    match_c = re.search(r'C(\d+)', formula)
    if not match_c:
        return {"error": "分子式解析失败。未找到碳原子数。"}

    try:
        n_c = int(match_c.group(1))
    except ValueError:
        return {"error": "碳原子数无效。"}
        
    contains_o = 'O' in formula

    # 生成原始异构体
    if 1 <= n_c <= 5:
        if contains_o:
            # 含氧，调度到 C1-C5 杂原子文件
            isomers = generate_isomers_c1_c5_hetero(formula)
        else:
            # 纯烃，调度到 C1-C5 纯烃文件
            isomers = generate_isomers_c1_c5(formula)
            
    elif n_c == 6:
        if contains_o:
            # 含氧，调度到 C6 杂原子文件
            isomers = generate_isomers_c6_hetero(formula)
        else:
            # 获取氢数
            match_h = re.search(r'H(\d*)', formula)
            n_h = int(match_h.group(1)) if match_h and match_h.group(1) else 0
            
            # C6纯烃，特殊处理C6H6（使用通用枚举算法，包含所有127种异构体）
            if n_h == 6:  # C6H6 - 苯及所有可能的价键异构体，使用通用枚举生成器
                isomers = generate_c6h6_general(formula)
            else:
                # 其他C6纯烃使用增强算法
                isomers = generate_c6_enhanced(formula)
            
    elif n_c == 7:
        # C7纯烃，特殊处理C7H16（使用已知的9种异构体）
        if not contains_o:
            isomers = generate_isomers_general(formula)
        else:
            # 含氧C7使用通用算法
            isomers = generate_isomers_general(formula)
            
    elif n_c == 8:
        # C8纯烃，特殊处理C8H18（使用算法生成18种异构体）
        if not contains_o:
            # 获取氢数以区分烷烃和烯烃/炔烃
            match_h = re.search(r'H(\d*)', formula)
            n_h = int(match_h.group(1)) if match_h and match_h.group(1) else 0
            
            # 只有C8H18（辛烷）使用专用算法，其他C8化合物使用通用算法
            if n_h == 18:  # C8H18 - 辛烷
                try:
                    # 使用C8-C10生成器
                    from generators.isomer_c8_c10_combined import generate_isomers_c8_c10
                    isomers = generate_isomers_c8_c10(formula)
                    print("使用C8-C10生成器")
                except Exception as e:
                    print(f"警告: C8H18生成器出错: {e}，使用通用算法")
                    isomers = generate_isomers_general(formula)
            else:
                # C8H16等非烷烃使用通用算法
                print(f"C8化合物{formula}使用通用算法")
                isomers = generate_isomers_general(formula)
        else:
            # 含氧C8使用通用算法
            isomers = generate_isomers_general(formula)
            
    elif n_c == 9:
        # C9纯烃，特殊处理C9H20（使用算法生成35种异构体）
        if not contains_o:
            # 获取氢数以区分烷烃和其他化合物
            match_h = re.search(r'H(\d*)', formula)
            n_h = int(match_h.group(1)) if match_h and match_h.group(1) else 0
            
            # 只有C9H20（壬烷）使用专用算法，其他C9化合物使用通用算法
            if n_h == 20:  # C9H20 - 壬烷
                try:
                    from generators.isomer_c8_c10_combined import generate_isomers_c8_c10
                    isomers = generate_isomers_c8_c10(formula)
                    print("使用C8-C10生成器")
                except Exception as e:
                    print(f"警告: C9H20生成器出错: {e}，使用通用算法")
                    isomers = generate_isomers_general(formula)
            else:
                # C9H18等非烷烃使用通用算法
                print(f"C9化合物{formula}使用通用算法")
                isomers = generate_isomers_general(formula)
        else:
            # 含氧C9使用通用算法
            isomers = generate_isomers_general(formula)
            
    elif n_c == 10:
        # C10纯烃，特殊处理C10H22（使用算法生成75种异构体）
        if not contains_o:
            # 获取氢数以区分烷烃和其他化合物
            match_h = re.search(r'H(\d*)', formula)
            n_h = int(match_h.group(1)) if match_h and match_h.group(1) else 0
            
            # 只有C10H22（癸烷）使用专用算法，其他C10化合物使用通用算法
            if n_h == 22:  # C10H22 - 癸烷
                try:
                    from generators.isomer_c8_c10_combined import generate_isomers_c8_c10
                    isomers = generate_isomers_c8_c10(formula)
                    print("使用C8-C10生成器")
                except Exception as e:
                    print(f"警告: C10H22生成器出错: {e}，使用通用算法")
                    isomers = generate_isomers_general(formula)
            else:
                # C10H20等非烷烃使用通用算法
                print(f"C10化合物{formula}使用通用算法")
                isomers = generate_isomers_general(formula)
        else:
            # 含氧C10使用通用算法
            isomers = generate_isomers_general(formula)
            
    elif n_c >= 11:
        # C11及以上使用通用算法
        isomers = generate_isomers_general(formula)
            
    else:
        return {"error": f"不支持的碳原子数。输入碳数: {n_c}"}
    
    # 检查是否需要结构筛选
    if structure_prefs and isinstance(isomers, list):
        filter = StructureFilter()
        filtered_isomers = filter.filter_isomers(isomers, structure_prefs)
        
        # 如果筛选后没有结果，给出提示
        if not filtered_isomers:
            return {
                "error": f"根据指定的结构筛选条件，未找到符合条件的异构体。\n"
                        f"原始生成 {len(isomers)} 个异构体，请尝试放宽筛选条件。"
            }
        
        return filtered_isomers
    
    return isomers


# ==============================================================================
# GUI 界面 (Tkinter) - (其余代码保持不变)
# ==============================================================================

class IsomerAnalyzerGUI:
    def __init__(self, master: tk.Tk) -> None:
        self.master = master
        master.title("MolGen+ 有机分子同分异构体分析工具")
        master.geometry("1200x800")  # 设置默认窗口大小
        
        self.formula_var = tk.StringVar(master, value="C3H6O2") 
        self.results_data: list[dict[str, Any]] = []
        self.image_references: list[Any] = []
        self.selected_isomers: list[int] = []  # 存储选中的异构体索引
        self.isomer_checkboxes: list[tk.Checkbutton] = []  # 存储复选框组件 

        # 创建主窗口布局 - 使用PanedWindow实现左右分栏
        main_paned = tk.PanedWindow(master, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # 左侧面板
        left_panel = tk.Frame(main_paned)
        right_panel = tk.Frame(main_paned)
        
        main_paned.add(left_panel, minsize=400)
        main_paned.add(right_panel, minsize=350)

        # === 左侧面板内容 ===
        # --- 输入控制区域 ---
        input_frame = tk.Frame(left_panel, padx=10, pady=10)
        input_frame.pack(fill=tk.X)
        
        # 第一行：分子式输入和分析按钮
        input_row1 = tk.Frame(input_frame)
        input_row1.pack(fill=tk.X, pady=5)
        
        tk.Label(input_row1, text="分子式 (C+, 兼容 O):").pack(side=tk.LEFT, padx=5)
        self.formula_entry = tk.Entry(input_row1, textvariable=self.formula_var, width=15)
        self.formula_entry.pack(side=tk.LEFT, padx=5)
        
        self.analyze_button = tk.Button(input_row1, text="分析同分异构体", command=self.analyze_isomers, bg="#4CAF50", fg="white")
        self.analyze_button.pack(side=tk.LEFT, padx=5)

        # 合并分子库和快速搜索按钮
        self.library_button = tk.Button(input_row1, text="🔍 分子库搜索", command=self.open_molecule_library, bg="#2196F3", fg="white")
        self.library_button.pack(side=tk.LEFT, padx=5)

        # 机器学习分析按钮
        if HAS_ML_ANALYZER:
            self.ml_analysis_button = tk.Button(input_row1, text="🤖 机器学习分析", command=self.open_ml_analysis_window, bg="#9C27B0", fg="white")
            self.ml_analysis_button.pack(side=tk.LEFT, padx=5)
        
        # 第二行：结构选择面板
        structure_frame = tk.LabelFrame(input_frame, text="结构筛选 (可选)", padx=10, pady=5)
        structure_frame.pack(fill=tk.X, pady=10)
        
        # 结构变量
        self.structure_vars = {
            'ring_type': tk.StringVar(value='all'),  # all, acyclic, monocyclic, polycyclic
            'bond_type': tk.StringVar(value='all'),   # all, single, double, triple, mixed
            'aromatic': tk.StringVar(value='all'),     # all, none, benzene, aromatic
            'functional_groups': tk.StringVar(value='all')  # all, alcohol, carbonyl, carboxyl, ester
        }
        
        # 第一列：环结构
        col1 = tk.Frame(structure_frame)
        col1.pack(side=tk.LEFT, padx=10, fill=tk.Y)
        tk.Label(col1, text="环结构:", font=('Arial', 9, 'bold')).pack(anchor=tk.W)
        ring_options = [('全部', 'all'), ('无环', 'acyclic'), ('单环', 'monocyclic'), ('多环', 'polycyclic')]
        for text, value in ring_options:
            tk.Radiobutton(col1, text=text, variable=self.structure_vars['ring_type'], 
                          value=value, font=('Arial', 8)).pack(anchor=tk.W)
        
        # 第二列：键类型
        col2 = tk.Frame(structure_frame)
        col2.pack(side=tk.LEFT, padx=10, fill=tk.Y)
        tk.Label(col2, text="键类型:", font=('Arial', 9, 'bold')).pack(anchor=tk.W)
        bond_options = [('全部', 'all'), ('单键', 'single'), ('双键', 'double'), ('三键', 'triple'), ('混合', 'mixed')]
        for text, value in bond_options:
            tk.Radiobutton(col2, text=text, variable=self.structure_vars['bond_type'], 
                          value=value, font=('Arial', 8)).pack(anchor=tk.W)
        
        # 第三列：芳香性
        col3 = tk.Frame(structure_frame)
        col3.pack(side=tk.LEFT, padx=10, fill=tk.Y)
        tk.Label(col3, text="芳香性:", font=('Arial', 9, 'bold')).pack(anchor=tk.W)
        aromatic_options = [('全部', 'all'), ('无芳香', 'none'), ('含苯环', 'benzene'), ('其他芳香', 'aromatic')]
        for text, value in aromatic_options:
            tk.Radiobutton(col3, text=text, variable=self.structure_vars['aromatic'], 
                          value=value, font=('Arial', 8)).pack(anchor=tk.W)
        
        # 第四列：官能团
        col4 = tk.Frame(structure_frame)
        col4.pack(side=tk.LEFT, padx=10, fill=tk.Y)
        tk.Label(col4, text="官能团:", font=('Arial', 9, 'bold')).pack(anchor=tk.W)
        fg_options = [('全部', 'all'), ('羟基(-OH)', 'alcohol'), ('羰基(C=O)', 'carbonyl'), 
                    ('羧基(-COOH)', 'carboxyl'), ('酯基(-COO-)', 'ester')]
        for text, value in fg_options:
            tk.Radiobutton(col4, text=text, variable=self.structure_vars['functional_groups'], 
                          value=value, font=('Arial', 8)).pack(anchor=tk.W)
        
        # 第二行：导出按钮组
        input_row2 = tk.Frame(input_frame)
        input_row2.pack(fill=tk.X)
        
        self.export_button = tk.Button(input_row2, text="导出全部 (.txt)", command=self.export_results, state=tk.DISABLED)
        self.export_button.pack(side=tk.LEFT, padx=2)
        
        self.export_selected_button = tk.Button(input_row2, text="导出选中 (.txt)", command=self.export_selected_results, state=tk.DISABLED)
        self.export_selected_button.pack(side=tk.LEFT, padx=2)
        
        self.export_gjf_button = tk.Button(input_row2, text="导出GaussView (.gjf)", command=self.export_gjf_format, state=tk.DISABLED)
        self.export_gjf_button.pack(side=tk.LEFT, padx=2)
        
        self.export_xyz_button = tk.Button(input_row2, text="导出XYZ (.xyz)", command=self.export_xyz_format, state=tk.DISABLED)
        self.export_xyz_button.pack(side=tk.LEFT, padx=2)

        # 机器学习快速分析按钮
        if HAS_ML_ANALYZER:
            self.ml_quick_button = tk.Button(input_row2, text="📊 ML快速分析", command=self.run_quick_ml_analysis, state=tk.DISABLED, bg="#9C27B0", fg="white")
            self.ml_quick_button.pack(side=tk.LEFT, padx=2)

        # DeepSeek AI 按钮
        if HAS_DEEPSEEK:
            self.ai_dedup_button = tk.Button(input_row2, text="🔄 AI去重", command=self.ai_deduplicate, state=tk.DISABLED, bg="#FF9800", fg="white")
            self.ai_dedup_button.pack(side=tk.LEFT, padx=2)

            self.ai_enhance_button = tk.Button(input_row2, text="✨ AI增强", command=self.ai_enhance, state=tk.DISABLED, bg="#E91E63", fg="white")
            self.ai_enhance_button.pack(side=tk.LEFT, padx=2)

        # --- 文本结果显示区域 ---
        text_frame = tk.Frame(left_panel, padx=10, pady=5)
        text_frame.pack(fill=tk.BOTH, expand=True)
        
        tk.Label(text_frame, text="分析结果:", font=("Arial", 10, "bold")).pack(anchor=tk.W)
        
        self.results_text = scrolledtext.ScrolledText(text_frame, wrap=tk.WORD, height=20, width=60, font=("Courier", 9))
        self.results_text.pack(fill=tk.BOTH, expand=True, pady=5)
        self.results_text.insert(tk.END, "请在上方输入分子式 (C1+) 并点击 '分析同分异构体'。\n\n注意：C7及以上使用通用算法生成，可能需要较长时间。")
        self.results_text.config(state=tk.DISABLED)

        # --- 选择控制区域 ---
        selection_frame = tk.Frame(left_panel, padx=10, pady=5)
        selection_frame.pack(fill=tk.X)
        
        self.select_all_button = tk.Button(selection_frame, text="全选", command=self.select_all_isomers)
        self.select_all_button.pack(side=tk.LEFT, padx=(0, 5))
        
        self.select_none_button = tk.Button(selection_frame, text="清空", command=self.select_none_isomers)
        self.select_none_button.pack(side=tk.LEFT, padx=(0, 5))
        
        self.selection_label = tk.Label(selection_frame, text="已选择: 0 个异构体", font=("Arial", 10))
        self.selection_label.pack(side=tk.LEFT, padx=10)

        # === 右侧面板内容 ===
        right_header = tk.Frame(right_panel, padx=10, pady=10)
        right_header.pack(fill=tk.X)
        tk.Label(right_header, text="分子图形展示", font=("Arial", 12, "bold")).pack()
        
        # 使用Listbox + Scrollbar的替代方案来显示分子图形
        self.create_molecule_display(right_panel)
        
    def create_molecule_display(self, parent):
        """创建分子显示区域，使用更可靠的滚动实现"""
        # 创建带滚动条的Frame
        container = tk.Frame(parent, padx=10, pady=0)
        container.pack(fill=tk.BOTH, expand=True)
        
        # 滚动条
        scrollbar = tk.Scrollbar(container)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # 创建Text widget作为容器，支持滚动
        self.molecule_text = tk.Text(container, wrap=tk.WORD, yscrollcommand=scrollbar.set, 
                                   width=40, bg='white', padx=5, pady=5)
        self.molecule_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.config(command=self.molecule_text.yview)
        
        # 禁用编辑，但允许选择
        self.molecule_text.config(state=tk.DISABLED)
        
        # 用于存储实际的分子组件
        self.molecule_widgets = []

    def analyze_isomers(self) -> None:
        formula = self.formula_var.get().strip().upper().replace(' ', '')

        # 获取结构偏好设置
        structure_prefs = {}
        structure_prefs = {}
        for key, var in self.structure_vars.items():
            structure_prefs[key] = var.get()
        
        self.results_data = []
        self.export_button.config(state=tk.DISABLED)
        self.export_selected_button.config(state=tk.DISABLED)
        self.export_gjf_button.config(state=tk.DISABLED)
        self.export_xyz_button.config(state=tk.DISABLED)

        # 禁用 DeepSeek AI 按钮
        if HAS_DEEPSEEK:
            self.ai_dedup_button.config(state=tk.DISABLED)
            self.ai_enhance_button.config(state=tk.DISABLED)

        self.image_references.clear()
        self.selected_isomers.clear()
        self.isomer_checkboxes.clear()
        self.molecule_widgets.clear()

        # 清空分子显示区域
        self.molecule_text.config(state=tk.NORMAL)
        self.molecule_text.delete(1.0, tk.END)
        
        self.results_text.config(state=tk.NORMAL)
        self.results_text.delete(1.0, tk.END)
        
        # 显示分析信息
        self.results_text.insert(tk.END, f"--- 正在分析分子式: {format_formula(formula)} ---\n")
        
        # 如果有结构筛选，显示筛选条件
        active_filters = [f"{k}: {v}" for k, v in structure_prefs.items() if v != 'all']
        if active_filters:
            self.results_text.insert(tk.END, f"结构筛选条件: {', '.join(active_filters)}\n")
        else:
            self.results_text.insert(tk.END, "无结构筛选条件\n")
        
        self.results_text.insert(tk.END, "\n")

        results = parse_formula_and_dispatch(formula, structure_prefs)
        
        if isinstance(results, dict) and "error" in results:
            messagebox.showerror("错误", results["error"])
            self.results_text.insert(tk.END, results["error"])
        else:
            self.results_data = results  # type: ignore
            self.display_results()
            if self.results_data:
                self.export_button.config(state=tk.NORMAL)
                self.export_selected_button.config(state=tk.NORMAL)
                self.export_gjf_button.config(state=tk.NORMAL)
                self.export_xyz_button.config(state=tk.NORMAL)

                # 启用ML分析按钮
                if HAS_ML_ANALYZER and len(self.results_data) >= 2:
                    self.ml_quick_button.config(state=tk.NORMAL)

                # 启用 DeepSeek AI 按钮
                if HAS_DEEPSEEK:
                    self.ai_dedup_button.config(state=tk.NORMAL)
                    self.ai_enhance_button.config(state=tk.NORMAL)


                # 显示筛选统计信息
                self._display_filter_statistics()

                # 显示烷烃理论统计信息（如果适用）
                self._display_alkane_statistics(formula)
        
        self.results_text.config(state=tk.DISABLED)
        self.molecule_text.config(state=tk.DISABLED)
        self.update_selection_label()
    
    def open_molecule_library(self) -> None:
        """打开分子库管理窗口"""
        try:
            # 优先使用优化的分子库GUI
            from database.molecule_library_gui_optimized import open_optimized_molecule_library
            open_optimized_molecule_library(self.master)
        except ImportError:
            # 回退到传统版本
            if HAS_LEGACY_GUI:
                try:
                    library_gui = MoleculeLibraryGUI(self.master)
                except Exception as e:
                    messagebox.showerror("错误", f"打开分子库失败: {e}")
            else:
                messagebox.showerror("错误", "分子库GUI模块不可用")
        except Exception as e:
            messagebox.showerror("错误", f"打开分子库失败: {e}")

    def open_ml_analysis_window(self) -> None:
        """打开机器学习分析窗口"""
        if not HAS_ML_ANALYZER:
            messagebox.showerror("错误", "机器学习模块不可用\n请安装: pip install scikit-learn matplotlib")
            return

        # 获取当前生成的SMILES列表
        current_smiles = self._get_current_smiles_list()

        if not current_smiles or len(current_smiles) < 2:
            messagebox.showwarning("提示", "请先生成至少2个异构体才能进行机器学习分析")
            return

        # 创建ML分析窗口
        ml_window = MLAnalysisWindow(self.master, current_smiles)

    def _get_current_smiles_list(self) -> list:
        """获取当前生成的SMILES列表"""
        if not self.results_data:
            return []

        return [item.get('SMILES', '') for item in self.results_data if item.get('SMILES')]

    def run_quick_ml_analysis(self) -> None:
        """快速运行ML分析并在结果区域显示"""
        if not HAS_ML_ANALYZER:
            messagebox.showerror("错误", "机器学习模块不可用")
            return

        current_smiles = self._get_current_smiles_list()

        if not current_smiles or len(current_smiles) < 2:
            messagebox.showwarning("提示", "请先生成至少2个异构体")
            return

        try:
            # 创建分析器并生成报告
            analyzer = MLIsomerAnalyzer()
            report = analyzer.generate_analysis_report(current_smiles)

            # 在结果文本框中显示ML分析结果
            self._display_ml_results(report)

            # 启用ML分析按钮
            if hasattr(self, 'ml_analysis_button'):
                self.ml_analysis_button.config(state=tk.NORMAL)

        except Exception as e:
            messagebox.showerror("错误", f"ML分析失败: {str(e)}")

    def _display_ml_results(self, report: dict) -> None:
        """在结果区域显示ML分析报告"""
        self.results_text.config(state=tk.NORMAL)
        self.results_text.insert(tk.END, "\n")
        self.results_text.insert(tk.END, "=" * 70 + "\n")
        self.results_text.insert(tk.END, "🤖 机器学习分析报告\n")
        self.results_text.insert(tk.END, "=" * 70 + "\n\n")

        # 显示摘要信息
        summary = report.get('summary', {})
        self.results_text.insert(tk.END, f"📊 分析摘要:\n")
        self.results_text.insert(tk.END, f"  总输入: {summary.get('total_input', 0)} 个分子\n")
        self.results_text.insert(tk.END, f"  有效分子: {summary.get('valid_molecules', 0)} 个\n")
        self.results_text.insert(tk.END, f"  无效分子: {summary.get('invalid_molecules', 0)} 个\n\n")

        # 显示聚类分析
        clustering = report.get('clustering', {})
        if clustering:
            self.results_text.insert(tk.END, f"🎯 聚类分析 ({clustering.get('method', 'N/A')}):\n")
            self.results_text.insert(tk.END, f"  聚类数量: {clustering.get('n_clusters', 0)}\n")
            if clustering.get('silhouette_score') is not None:
                self.results_text.insert(tk.END, f"  轮廓系数: {clustering['silhouette_score']:.4f}\n")
            self.results_text.insert(tk.END, "  聚类分布:\n")
            for cluster_name, count in clustering.get('cluster_sizes', {}).items():
                self.results_text.insert(tk.END, f"    {cluster_name}: {count} 个分子\n")
            self.results_text.insert(tk.END, "\n")

        # 显示多样性分析
        diversity = report.get('diversity', {})
        if diversity:
            self.results_text.insert(tk.END, f"🌈 多样性分析:\n")
            self.results_text.insert(tk.END, f"  多样性得分: {diversity.get('diversity_score', 0):.4f}\n")
            self.results_text.insert(tk.END, f"  平均成对距离: {diversity.get('mean_pairwise_distance', 0):.4f}\n\n")

        # 显示异常检测
        outliers = report.get('outliers', {})
        if outliers:
            self.results_text.insert(tk.END, f"⚠️  异常检测:\n")
            self.results_text.insert(tk.END, f"  异常分子数: {outliers.get('n_outliers', 0)}\n")
            self.results_text.insert(tk.END, f"  异常比例: {outliers.get('outlier_ratio', 0):.2%}\n\n")

        self.results_text.see(tk.END)
        self.results_text.config(state=tk.DISABLED)

    # ========================================
    # DeepSeek AI 功能方法
    # ========================================

    def ai_deduplicate(self) -> None:
        """AI 智能去重"""
        if not HAS_DEEPSEEK:
            messagebox.showerror("错误", "DeepSeek AI 模块不可用")
            return

        current_smiles = self._get_current_smiles_list()

        if not current_smiles or len(current_smiles) < 1:
            messagebox.showwarning("提示", "请先生成至少1个异构体")
            return

        try:
            result = ai_deduplicate_smiles(current_smiles, use_smart_analysis=False)

            # 显示结果
            self.results_text.config(state=tk.NORMAL)
            self.results_text.insert(tk.END, "\n")
            self.results_text.insert(tk.END, "=" * 70 + "\n")
            self.results_text.insert(tk.END, "🤖 DeepSeek AI 去重结果\n")
            self.results_text.insert(tk.END, "=" * 70 + "\n\n")
            self.results_text.insert(tk.END, f"原始数量: {len(current_smiles)} 个\n")
            self.results_text.insert(tk.END, f"去重后: {len(result['unique_smiles'])} 个\n")
            self.results_text.insert(tk.END, f"移除重复: {result['removed_count']} 个\n\n")

            if result['duplicates']:
                self.results_text.insert(tk.END, "重复项:\n")
                for norm_smiles, indices in result['duplicates'].items():
                    duplicates = [current_smiles[i] for i in indices]
                    self.results_text.insert(tk.END, f"  {norm_smiles}: {len(duplicates)} 个重复\n")
                self.results_text.insert(tk.END, "\n")

            self.results_text.insert(tk.END, "去重后的 SMILES:\n")
            for i, smi in enumerate(result['unique_smiles'], 1):
                self.results_text.insert(tk.END, f"  {i}. {smi}\n")

            self.results_text.see(tk.END)
            self.results_text.config(state=tk.DISABLED)

            messagebox.showinfo(
                "AI 去重完成",
                f"原始: {len(current_smiles)} 个\n去重后: {len(result['unique_smiles'])} 个\n移除: {result['removed_count']} 个"
            )

        except Exception as e:
            messagebox.showerror("错误", f"AI 去重失败: {str(e)}")

    def ai_enhance(self) -> None:
        """AI 增强异构体"""
        if not HAS_DEEPSEEK:
            messagebox.showerror("错误", "DeepSeek AI 模块不可用")
            return

        current_smiles = self._get_current_smiles_list()
        formula = self.formula_var.get().strip()

        if not current_smiles or len(current_smiles) < 1:
            messagebox.showwarning("提示", "请先生成至少1个异构体")
            return

        try:
            self.results_text.config(state=tk.NORMAL)
            self.results_text.insert(tk.END, "\n")
            self.results_text.insert(tk.END, "=" * 70 + "\n")
            self.results_text.insert(tk.END, "✨ DeepSeek AI 增强分析\n")
            self.results_text.insert(tk.END, "=" * 70 + "\n")
            self.results_text.insert(tk.END, "\n正在调用 DeepSeek API，请稍候...\n")
            self.results_text.update()

            result = ai_enhance_isomer_generation(formula, current_smiles)

            # 清除等待消息
            lines = self.results_text.get("1.0", tk.END).split('\n')
            if '正在调用 DeepSeek API，请稍候...' in lines[-2]:
                self.results_text.delete("end-2l linestart", "end-1l")

            # 显示结果
            self.results_text.insert(tk.END, f"\n原始数量: {len(current_smiles)} 个\n")
            self.results_text.insert(tk.END, f"AI 生成新候选: {len(result['new_candidates'])} 个\n")
            self.results_text.insert(tk.END, f"最终数量: {len(result['final_smiles'])} 个\n\n")

            if result['new_candidates']:
                self.results_text.insert(tk.END, "🎯 AI 生成的新候选:\n")
                for i, smi in enumerate(result['new_candidates'], 1):
                    self.results_text.insert(tk.END, f"  {i}. {smi}\n")
                self.results_text.insert(tk.END, "\n")

            self.results_text.insert(tk.END, f"所有异构体（共 {len(result['final_smiles'])} 个）:\n")
            for i, smi in enumerate(result['final_smiles'], 1):
                self.results_text.insert(tk.END, f"  {i}. {smi}\n")

            self.results_text.see(tk.END)
            self.results_text.config(state=tk.DISABLED)

            messagebox.showinfo(
                "AI 增强完成",
                f"原始: {len(current_smiles)} 个\n新增: {len(result['new_candidates'])} 个\n最终: {len(result['final_smiles'])} 个"
            )

        except Exception as e:
            messagebox.showerror("错误", f"AI 增强失败: {str(e)}")

    def _display_filter_statistics(self) -> None:
        """显示筛选统计信息"""
        if not self.results_data:
            return
            
        try:
            from analysis.structure_filter import StructureFilter
            filter = StructureFilter()
            stats = filter.get_structure_statistics(self.results_data)
            
            if stats:
                self.results_text.config(state=tk.NORMAL)
                self.results_text.insert(tk.END, f"\n--- 筛选结果统计 ---\n")
                self.results_text.insert(tk.END, f"总计找到: {stats['total_count']} 个异构体\n")
                
                # 显示主要统计
                if any(count > 0 for count in stats['ring_types'].values()):
                    self.results_text.insert(tk.END, f"\n环结构分布:\n")
                    for ring_type, count in stats['ring_types'].items():
                        if count > 0:
                            self.results_text.insert(tk.END, f"  {ring_type}: {count} 个\n")
                
                if any(count > 0 for count in stats['aromatic_types'].values()):
                    self.results_text.insert(tk.END, f"\n芳香性分布:\n")
                    for aromatic_type, count in stats['aromatic_types'].items():
                        if count > 0:
                            self.results_text.insert(tk.END, f"  {aromatic_type}: {count} 个\n")
                
                self.results_text.config(state=tk.DISABLED)
                
        except Exception as e:
            print(f"显示统计信息时出错: {e}", file=sys.stderr)
    
    def _display_alkane_statistics(self, formula: str) -> None:
        """显示烷烃同分异构体理论统计信息"""
        if not HAS_ALKANE_COUNTER:
            return
        
        # 检查是否为纯烷烃分子式
        import re
        if not re.match(r'^C\d+H\d+$', formula):
            return
        
        try:
            stats = get_alkane_isomer_statistics(formula)
            if "error" in stats:
                return
            
            self.results_text.config(state=tk.NORMAL)
            self.results_text.insert(tk.END, f"\n--- 烷烃理论统计 ---\n")
            self.results_text.insert(tk.END, f"理论异构体总数: {stats['isomer_count']:,}\n")
            self.results_text.insert(tk.END, f"当前生成数量: {len(self.results_data)}\n")
            
            if stats['isomer_count'] > 0:
                coverage = (len(self.results_data) / stats['isomer_count'] * 100)
                self.results_text.insert(tk.END, f"覆盖率: {coverage:.1f}%\n")
                
                missing = max(0, stats['isomer_count'] - len(self.results_data))
                if missing > 0:
                    self.results_text.insert(tk.END, f"未生成: {missing} 个异构体\n")
                else:
                    self.results_text.insert(tk.END, f"✓ 已生成所有理论异构体\n")
            
            self.results_text.insert(tk.END, f"异构体密度: {stats['density']:.3f} 个/碳原子\n")
            self.results_text.config(state=tk.DISABLED)
            
        except Exception as e:
            print(f"显示烷烃统计时出错: {e}", file=sys.stderr)

    def display_results(self) -> None:
        total_isomers = len(self.results_data)
        self.results_text.insert(tk.END, f"共找到 {total_isomers} 种同分异构体 (基于优化的预设结构集)\n\n")
        
        for i, item in enumerate(self.results_data):
            # 1. 文本输出
            header = f"============= 异构体 {i+1} =============\n"
            self.results_text.insert(tk.END, header)
            self.results_text.insert(tk.END, f"SMILES: {item.get('SMILES', '')}\n")
            self.results_text.insert(tk.END, f"结构类型: {item.get('StructureType', '')}\n")
            self.results_text.insert(tk.END, f"\n--- 原子坐标 ({item.get('Dimension', '')} Angstrom) ---\n")
            self.results_text.insert(tk.END, item.get('Coords', ''))
            self.results_text.insert(tk.END, "\n\n") 
        
        self.results_text.see(1.0)

        # 2. 启用按钮
        if total_isomers > 0:
            self.export_button.config(state=tk.NORMAL)
            self.export_selected_button.config(state=tk.NORMAL)
            self.export_gjf_button.config(state=tk.NORMAL)
            self.export_xyz_button.config(state=tk.NORMAL)

            # 启用 DeepSeek AI 按钮
            if HAS_DEEPSEEK:
                self.ai_dedup_button.config(state=tk.NORMAL)
                self.ai_enhance_button.config(state=tk.NORMAL)

        # 3. 绘制所有分子图形
        if total_isomers > 0 and HAS_DRAWING:
            self.draw_all_molecules()

    def draw_all_molecules(self) -> None:
        # 清空之前的分子显示
        self.molecule_widgets.clear()
        self.molecule_text.config(state=tk.NORMAL)
        self.molecule_text.delete(1.0, tk.END)
        
        if not HAS_DRAWING:
            self.molecule_text.insert(tk.END, "警告: 请安装 Pillow 库以显示分子图形\n")
            self.molecule_text.config(state=tk.DISABLED)
            return
        
        # 创建主框架来包含所有分子
        main_frame = tk.Frame(self.molecule_text, bg='white')
        
        for i, item in enumerate(self.results_data):
            mol = item.get('Mol')
            smiles = item.get('SMILES', '')
            
            # 为每个异构体创建一个框架
            isomer_frame = tk.Frame(main_frame, bg='white', relief=tk.RIDGE, bd=2)
            isomer_frame.pack(pady=10, padx=5, fill=tk.X)
            
            # 添加复选框和标题
            header_frame = tk.Frame(isomer_frame, bg='white')
            header_frame.pack(fill=tk.X, padx=5, pady=5)
            
            var = tk.BooleanVar(value=False)
            checkbox = tk.Checkbutton(header_frame, text=f"异构体 {i+1}", 
                                     variable=var, 
                                     command=lambda idx=i, v=var: self.toggle_isomer_selection(idx, v))
            checkbox.pack(side=tk.LEFT)
            checkbox.variable = var  # 保存变量引用
            self.isomer_checkboxes.append(checkbox)
            
            # 分子信息
            info_text = f"{item.get('StructureType', '')}"
            info_label = tk.Label(header_frame, text=info_text, font=('Arial', 9, 'italic'), bg='white')
            info_label.pack(side=tk.LEFT, padx=(10, 0))
            
            # SMILES 标签
            smiles_label = tk.Label(isomer_frame, text=f"SMILES: {smiles}", 
                                   font=('Arial', 8), bg='white', fg='blue')
            smiles_label.pack(anchor=tk.W, padx=5, pady=5)
            
            # 分子图像
            if mol:
                try:
                    pil_img = MolToImage(mol, size=(280, 280))
                    tk_img = ImageTk.PhotoImage(pil_img)  # type: ignore
                    self.image_references.append(tk_img) 
                    
                    image_label = tk.Label(isomer_frame, image=tk_img, bg='white')
                    image_label.pack(pady=5)
                    
                except Exception as e:
                    msg = str(e)
                    if 'Cairo' in msg or 'MolToImage requires' in msg:
                        hint = "请安装: conda install -c conda-forge rdkit pillow"
                        error_text = f"绘图失败: RDKit 未启用 Cairo 支持\n{hint}"
                    else:
                        error_text = f"绘图错误: {e}"
                    
                    error_label = tk.Label(isomer_frame, text=error_text, 
                                         fg="red", bg='white', font=('Arial', 8), 
                                         justify=tk.LEFT, wraplength=280)
                    error_label.pack(pady=10)
                    print(f"异构体 #{i+1} 绘图失败: {e}", file=sys.stderr)
            else:
                no_img_label = tk.Label(isomer_frame, text="无法生成分子图像", 
                                       fg="gray", bg='white', font=('Arial', 9))
                no_img_label.pack(pady=10)
        
        # 将主框架添加到Text widget
        self.molecule_text.window_create(tk.END, window=main_frame)
        self.molecule_text.config(state=tk.DISABLED)
        
        # 滚动到顶部
        self.molecule_text.see(tk.END)

    def toggle_isomer_selection(self, index: int, var: tk.BooleanVar) -> None:
        """切换异构体选择状态"""
        if var.get():
            if index not in self.selected_isomers:
                self.selected_isomers.append(index)
        else:
            if index in self.selected_isomers:
                self.selected_isomers.remove(index)
        self.update_selection_label()

    def update_selection_label(self) -> None:
        """更新选择状态标签"""
        count = len(self.selected_isomers)
        self.selection_label.config(text=f"已选择: {count} 个异构体")

    def select_all_isomers(self) -> None:
        """全选所有异构体"""
        self.selected_isomers = list(range(len(self.results_data)))
        for i, checkbox in enumerate(self.isomer_checkboxes):
            if hasattr(checkbox, 'variable'):
                checkbox.variable.set(True)
        self.update_selection_label()

    def select_none_isomers(self) -> None:
        """清空所有选择"""
        self.selected_isomers.clear()
        for checkbox in self.isomer_checkboxes:
            if hasattr(checkbox, 'variable'):
                checkbox.variable.set(False)
        self.update_selection_label()

    def export_selected_results(self) -> None:
        """导出选中的异构体"""
        if not self.results_data:
            messagebox.showwarning("警告", "没有可导出的数据。请先执行分析。")
            return
        
        if not self.selected_isomers:
            messagebox.showwarning("警告", "请先选择要导出的异构体。")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="保存选中的同分异构体分析结果"
        )

        if not file_path:
            return
        
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(f"*** MolGen+ 选中的同分异构体分析结果 ***\n")
                f.write(f"分子式: {format_formula(self.formula_var.get())}\n")
                f.write(f"文件名: {os.path.basename(file_path)}\n")
                f.write(f"生成时间: {tk.Tcl().eval('clock format [clock seconds]')}\n")
                f.write(f"总异构体数: {len(self.results_data)}\n")
                f.write(f"选中导出: {len(self.selected_isomers)} 个异构体\n")
                f.write(f"选中序号: {[i+1 for i in sorted(self.selected_isomers)]}\n\n")

                for idx in sorted(self.selected_isomers):
                    if 0 <= idx < len(self.results_data):
                        item = self.results_data[idx]
                        f.write(f"=========================================\n")
                        f.write(f"异构体 {idx+1} (已选中):\n")
                        f.write(f"SMILES: {item.get('SMILES', '')}\n")
                        f.write(f"InChIKey: {item.get('InChIKey', '')}\n")
                        f.write(f"结构类型: {item.get('StructureType', '')}\n")
                        f.write(f"-----------------------------------------\n")
                        f.write(f"原子坐标 ({item.get('Dimension', '')} Angstrom):\n")
                        f.write(item.get('Coords', '') + "\n")
                        f.write("=========================================\n\n")
            
            messagebox.showinfo("成功", f"选中的 {len(self.selected_isomers)} 个异构体已成功导出至:\n{file_path}")

        except Exception as e:
            messagebox.showerror("导出错误", f"导出文件时发生错误: {e}")

    def parse_coordinates(self, coords_str: str) -> list[dict[str, Any]]:
        """解析坐标字符串，返回原子列表"""
        atoms = []
        lines = coords_str.strip().split('\n')
        
        for line in lines:
            line = line.strip()
            if not line or not line.startswith('Atom'):
                continue
            
            # 解析格式: "Atom 1 (C): x=1.2345, y=2.3456, z=3.4567"
            try:
                # 提取原子符号
                if '(' in line and ')' in line:
                    symbol_part = line.split('(')[1].split(')')[0]
                    atom_symbol = symbol_part
                else:
                    continue
                
                # 提取坐标
                x = y = z = 0.0
                if 'x=' in line:
                    x_str = line.split('x=')[1].split(',')[0].strip()
                    x = float(x_str)
                if 'y=' in line:
                    y_str = line.split('y=')[1].split(',')[0].strip()
                    y = float(y_str)
                if 'z=' in line:
                    z_str = line.split('z=')[1].strip()
                    z = float(z_str)
                else:
                    # 2D结构，设置z=0
                    z = 0.0
                
                atoms.append({
                    'symbol': atom_symbol,
                    'x': x,
                    'y': y,
                    'z': z
                })
            except (ValueError, IndexError) as e:
                print(f"解析坐标行失败: {line}, 错误: {e}", file=sys.stderr)
                continue
        
        return atoms

    def export_gjf_format(self) -> None:
        """导出Gaussian输入格式(.gjf)，兼容GaussView"""
        if not self.results_data:
            messagebox.showwarning("警告", "没有可导出的数据。请先执行分析。")
            return
        
        if not self.selected_isomers:
            messagebox.showwarning("警告", "请先选择要导出的异构体。")
            return

        # 如果只选中一个异构体，导出单个文件；如果选中多个，分别导出多个文件
        if len(self.selected_isomers) == 1:
            self.export_single_gjf()
        else:
            self.export_multiple_gjf()

    def detect_bond_pattern(self, atoms: list[dict[str, Any]], smiles: str = "") -> list[dict[str, Any]]:
        """检测键连接模式，使用距离检测方法"""
        bonds = []
        
        # 检测苯环：6个碳原子且所有C-C键长相近（约1.4Å）
        if self.is_benzene_ring(atoms):
            bonds = self.create_benzene_bonds(atoms)
        else:
            bonds = self.create_general_bonds(atoms)
            
        return bonds
    

    
    def is_benzene_ring(self, atoms: list[dict[str, Any]]) -> bool:
        """判断是否为苯环结构：6个碳原子，键长基本一致"""
        # 检查是否有6个碳原子（可能包含氢原子）
        carbon_atoms = [a for a in atoms if a['symbol'] == 'C']
        if len(carbon_atoms) != 6:
            return False
        
        # 检查是否有6个氢原子（苯的分子式C6H6）
        hydrogen_atoms = [a for a in atoms if a['symbol'] == 'H']
        if len(hydrogen_atoms) != 6:
            return False
        
        # 计算所有碳-碳键长（只计算碳原子之间的距离）
        cc_distances = []
        for i in range(len(carbon_atoms)):
            for j in range(i+1, len(carbon_atoms)):
                dist = ((carbon_atoms[i]['x'] - carbon_atoms[j]['x'])**2 + 
                       (carbon_atoms[i]['y'] - carbon_atoms[j]['y'])**2 + 
                       (carbon_atoms[i]['z'] - carbon_atoms[j]['z'])**2)**0.5
                if dist < 2.0:  # 只考虑小于2Å的C-C距离
                    cc_distances.append(dist)
        
        # 苯环应该有6个C-C键，键长应该在1.4Å附近
        if len(cc_distances) == 6:
            avg_distance = sum(cc_distances) / len(cc_distances)
            # 检查键长是否在苯环范围内，且变化很小
            if 1.35 < avg_distance < 1.45:  # 苯环C-C键长约1.395Å
                max_deviation = max(abs(d - avg_distance) for d in cc_distances)
                if max_deviation < 0.1:  # 所有键长变化小于0.1Å
                    return True
        
        return False
    
    def create_benzene_bonds(self, atoms: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """为苯环创建芳香键连接"""
        bonds = []
        
        # 找到最接近的6个碳-碳连接（成环）
        distances = []
        for i in range(len(atoms)):
            for j in range(i+1, len(atoms)):
                dist = ((atoms[i]['x'] - atoms[j]['x'])**2 + 
                       (atoms[i]['y'] - atoms[j]['y'])**2 + 
                       (atoms[i]['z'] - atoms[j]['z'])**2)**0.5
                distances.append((dist, i, j))
        
        # 按距离排序，取前6个最小的距离（成环的键）
        distances.sort()
        for dist, i, j in distances[:6]:
            bonds.append({
                'atom1': i + 1,
                'atom2': j + 1,
                'bond_order': 1.5  # 芳香键
            })
        
        return bonds
    
    def create_general_bonds(self, atoms: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """为一般分子创建键连接，简化和修复的键检测算法"""
        bonds = []
        
        # 标准键长阈值表（保守设置，确保不错过键）
        bond_thresholds = {
            ('H', 'H'): 1.3,
            ('H', 'C'): 1.4,
            ('H', 'N'): 1.3,
            ('H', 'O'): 1.5,  # 扩大O-H键检测范围，确保RDKit生成的O-H键不被遗漏
            ('H', 'S'): 1.6,
            ('C', 'C'): 1.9,
            ('C', 'N'): 1.8,
            ('C', 'O'): 1.8,
            ('C', 'S'): 2.1,
            ('N', 'N'): 1.7,
            ('N', 'O'): 1.7,
            ('O', 'O'): 1.7,
            ('N', 'S'): 2.1,
            ('O', 'S'): 2.1,
        }
        
        # 遍历所有原子对，检测键
        for i in range(len(atoms)):
            for j in range(i+1, len(atoms)):
                symbol1 = atoms[i]['symbol']
                symbol2 = atoms[j]['symbol']
                
                # 计算原子间距离
                dist = ((atoms[i]['x'] - atoms[j]['x'])**2 + 
                       (atoms[i]['y'] - atoms[j]['y'])**2 + 
                       (atoms[i]['z'] - atoms[j]['z'])**2)**0.5
                
                # 获取对应原子对的键长阈值
                key = tuple(sorted([symbol1, symbol2]))
                bond_threshold = bond_thresholds.get(key, 2.0)  # 默认阈值2.0Å
                
                # 判断是否成键
                if dist <= bond_threshold:
                    # 确定键序
                    bond_order = self.determine_bond_order(symbol1, symbol2, dist)
                    
                    bonds.append({
                        'atom1': i + 1,
                        'atom2': j + 1,
                        'bond_order': bond_order,
                        'distance': dist,
                        'symbol1': symbol1,
                        'symbol2': symbol2
                    })
                    
                    # 输出调试信息
                    print(f"检测到键: {symbol1}-{symbol2}, 距离={dist:.3f}Å, 键序={bond_order}", file=sys.stderr)
        
        # 应用羧基键修正，确保O-H键正确检测
        bonds = self.fix_carboxyl_bonds(bonds, atoms)
        
        return bonds
    

    
    def fix_carboxyl_bonds(self, potential_bonds: list[dict[str, Any]], atoms: list[dict[str, Any]]) -> list[dict[str, Any]]:
        """修正羧基中的键序，确保一个C=O双键和一个C-O单键，并检测O-H键"""
        bonds = potential_bonds.copy()
        
        # 统计O-H键
        oh_bonds_before = [bond for bond in bonds if (bond['symbol1'] == 'O' and bond['symbol2'] == 'H') or (bond['symbol1'] == 'H' and bond['symbol2'] == 'O')]
        print(f"fix_carboxyl_bonds 开始: 检测到 {len(oh_bonds_before)} 个O-H键", file=sys.stderr)
        
        # 首先确保所有O-H键都被正确检测，扩大搜索范围确保羧基O-H键不被遗漏
        for i in range(len(atoms)):
            for j in range(i+1, len(atoms)):
                symbol1 = atoms[i]['symbol']
                symbol2 = atoms[j]['symbol']
                
                # 检查O-H键
                if (symbol1 == 'O' and symbol2 == 'H') or (symbol1 == 'H' and symbol2 == 'O'):
                    dist = ((atoms[i]['x'] - atoms[j]['x'])**2 + 
                           (atoms[i]['y'] - atoms[j]['y'])**2 + 
                           (atoms[i]['z'] - atoms[j]['z'])**2)**0.5
                    
                    # 扩大O-H键检测范围，包含可能的羧基O-H键
                    if 0.80 <= dist <= 1.40:  # 扩大范围以兼容RDKit生成的坐标
                        # 检查是否已存在这个键
                        atom1_idx = i + 1
                        atom2_idx = j + 1
                        exists = False
                        for bond in bonds:
                            if ((bond['atom1'] == atom1_idx and bond['atom2'] == atom2_idx) or
                                (bond['atom1'] == atom2_idx and bond['atom2'] == atom1_idx)):
                                exists = True
                                bond['bond_order'] = 1.0  # 确保是单键
                                break
                        
                        if not exists:
                            bonds.append({
                                'atom1': atom1_idx,
                                'atom2': atom2_idx,
                                'bond_order': 1.0,
                                'distance': dist,
                                'symbol1': symbol1,
                                'symbol2': symbol2
                            })
                            print(f"检测到O-H键: 原子{atom1_idx}({symbol1})-原子{atom2_idx}({symbol2}), 距离={dist:.3f}Å", file=sys.stderr)
        
        oh_bonds_after_first_pass = [bond for bond in bonds if (bond['symbol1'] == 'O' and bond['symbol2'] == 'H') or (bond['symbol1'] == 'H' and bond['symbol2'] == 'O')]
        print(f"fix_carboxyl_bonds 第一次扫描后: {len(oh_bonds_after_first_pass)} 个O-H键", file=sys.stderr)
        
        # 查找所有C-O键
        co_bonds = [bond for bond in bonds if 
                   ((bond['symbol1'] == 'C' and bond['symbol2'] == 'O') or 
                    (bond['symbol1'] == 'O' and bond['symbol2'] == 'C'))]
        
        # 如果有多个C-O键连接到同一个碳原子，可能是羧基
        if len(co_bonds) >= 2:
            # 找出连接到同一个碳原子的C-O键对
            for i, bond1 in enumerate(co_bonds):
                for bond2 in co_bonds[i+1:]:
                    # 找出碳原子
                    if bond1['symbol1'] == 'C':
                        c_atom1 = bond1['atom1']
                        o_atom1 = bond1['atom2']
                    else:
                        c_atom1 = bond1['atom2']
                        o_atom1 = bond1['atom1']
                    
                    if bond2['symbol1'] == 'C':
                        c_atom2 = bond2['atom1']
                        o_atom2 = bond2['atom2']
                    else:
                        c_atom2 = bond2['atom2']
                        o_atom2 = bond2['atom1']
                    
                    # 如果两个C-O键连接到同一个碳原子
                    if c_atom1 == c_atom2:
                        # 检查哪个氧原子连接了氢（形成O-H键）
                        o1_has_h = False
                        o2_has_h = False
                        
                        for bond in bonds:
                            if ((bond['atom1'] == o_atom1 or bond['atom2'] == o_atom1) and
                                ((atoms[bond['atom1']-1]['symbol'] == 'H' and atoms[bond['atom2']-1]['symbol'] == 'O') or
                                 (atoms[bond['atom1']-1]['symbol'] == 'O' and atoms[bond['atom2']-1]['symbol'] == 'H'))):
                                o1_has_h = True
                            if ((bond['atom1'] == o_atom2 or bond['atom2'] == o_atom2) and
                                ((atoms[bond['atom1']-1]['symbol'] == 'H' and atoms[bond['atom2']-1]['symbol'] == 'O') or
                                 (atoms[bond['atom1']-1]['symbol'] == 'O' and atoms[bond['atom2']-1]['symbol'] == 'H'))):
                                o2_has_h = True
                        
                        # 修正逻辑：连接氢的氧原子总是形成单键，不连接氢的形成双键
                        if o1_has_h and not o2_has_h:
                            # o1连接氢，所以o1是单键，o2是双键
                            single_bond = bond1 if bond1['atom1'] == o_atom1 or bond1['atom2'] == o_atom1 else bond2
                            double_bond = bond2 if bond1['atom1'] == o_atom1 or bond1['atom2'] == o_atom1 else bond1
                        elif o2_has_h and not o1_has_h:
                            # o2连接氢，所以o2是单键，o1是双键
                            single_bond = bond2 if bond2['atom1'] == o_atom2 or bond2['atom2'] == o_atom2 else bond1
                            double_bond = bond1 if bond2['atom1'] == o_atom2 or bond2['atom2'] == o_atom2 else bond2
                        else:
                            # 如果都没有连接氢或者都连接了氢，按距离判断（距离短的为双键）
                            if bond1['distance'] < bond2['distance']:
                                double_bond = bond1
                                single_bond = bond2
                            else:
                                double_bond = bond2
                                single_bond = bond1
                        
                        # 更新键序
                        double_bond['bond_order'] = 2.0
                        single_bond['bond_order'] = 1.0
                        
                        # 更新原始bonds列表中的对应项
                        for bond in bonds:
                            if ((bond['atom1'] == double_bond['atom1'] and bond['atom2'] == double_bond['atom2']) or
                                (bond['atom1'] == double_bond['atom2'] and bond['atom2'] == double_bond['atom1'])):
                                bond['bond_order'] = 2.0
                            elif ((bond['atom1'] == single_bond['atom1'] and bond['atom2'] == single_bond['atom2']) or
                                  (bond['atom1'] == single_bond['atom2'] and bond['atom2'] == single_bond['atom1'])):
                                bond['bond_order'] = 1.0
                        break
        
        oh_bonds_final = [bond for bond in bonds if (bond['symbol1'] == 'O' and bond['symbol2'] == 'H') or (bond['symbol1'] == 'H' and bond['symbol2'] == 'O')]
        print(f"fix_carboxyl_bonds 结束: 最终有 {len(oh_bonds_final)} 个O-H键", file=sys.stderr)
        
        return bonds
    
    def determine_bond_order(self, symbol1: str, symbol2: str, distance: float) -> float:
        """根据原子类型和距离判断键序，优化版本"""
        # 按原子对排序，便于比较
        pair = tuple(sorted([symbol1, symbol2]))
        
        # 更精确的键长范围表（最小值, 最大值, 对应键序）
        bond_ranges = {
            ('C', 'C'): [
                (1.15, 1.30, 2.0),   # 双键 (C=C)
                (1.35, 1.40, 1.5),   # 芳香键 (苯环)
                (1.45, 1.60, 1.0)    # 单键 (C-C)
            ],
            ('C', 'O'): [
                (1.15, 1.25, 2.0),   # 双键 (C=O)
                (1.30, 1.50, 1.0)    # 单键 (C-O)
            ],
            ('C', 'N'): [
                (1.20, 1.35, 2.0),   # 双键 (C=N)
                (1.35, 1.55, 1.0)    # 单键 (C-N)
            ],
            ('O', 'H'): [(0.85, 1.30, 1.0)],     # O-H单键（扩大范围以兼容RDKit生成的坐标）
            ('N', 'H'): [(0.85, 1.10, 1.0)],     # N-H单键
            ('C', 'H'): [(1.00, 1.20, 1.0)],     # C-H单键
            ('O', 'O'): [(1.40, 1.50, 1.0)],     # O-O单键 (罕见)
            ('N', 'N'): [(1.20, 1.35, 2.0), (1.40, 1.50, 1.0)],  # N=N双键, N-N单键
        }
        
        # 特殊处理：如果是有机分子中的C-O键，需要根据上下文判断
        if pair == ('C', 'O'):
            # C=O双键通常较短，C-O单键较长
            if distance < 1.25:
                return 2.0  # C=O双键
            else:
                return 1.0  # C-O单键
        
        # 其他键类型的判断
        if pair in bond_ranges:
            # 按优先级检查，先检查更严格的范围
            for min_dist, max_dist, order in bond_ranges[pair]:
                if min_dist <= distance <= max_dist:
                    return order
        
        # 默认返回单键
        return 1.0

    def export_single_gjf(self) -> None:
        """导出单个GJF文件，包含连接性信息以正确显示芳香性"""
        idx = self.selected_isomers[0]
        item = self.results_data[idx]

        file_path = filedialog.asksaveasfilename(
            defaultextension=".gjf",
            filetypes=[("Gaussian files", "*.gjf *.com"), ("All files", "*.*")],
            title=f"保存异构体{idx+1}的Gaussian输入文件",
            initialfile=f"{format_formula(self.formula_var.get())}_Isomer{idx+1}.gjf"
        )

        if not file_path:
            return
        
        try:
            base_name = os.path.splitext(os.path.basename(file_path))[0]
            atoms = self.parse_coordinates(item.get('Coords', ''))
            bonds = self.detect_bond_pattern(atoms, item.get('SMILES', ''))
            
            with open(file_path, 'w', encoding='utf-8') as f:
                # 标准GJF文件头部
                f.write(f"%chk={base_name}.chk\n")
                f.write("%mem=1GB\n")
                f.write("%nprocshared=4\n")
                f.write("#p B3LYP/6-31G(d) opt freq\n\n")
                f.write(f"{item.get('SMILES', '')} - {format_formula(self.formula_var.get())} - Isomer {idx+1}\n\n")
                f.write("0 1\n")  # 电荷和自旋多重度
                
                # 写入原子坐标
                for atom in atoms:
                    f.write(f"{atom['symbol']:>2}   {atom['x']:>10.6f}   {atom['y']:>10.6f}   {atom['z']:>10.6f}\n")
                
                f.write("\n")
                
                # 写入连接性信息（重要：确保GaussView正确显示）
                if bonds:
                    # 按GaussView标准格式写入连接性
                    # 格式：每行以原子序号开头，后面跟着相连的原子索引和键序
                    # 注意：每个键连接只写一次，不重复写
                    connectivity_dict = {}
                    used_connections = set()  # 记录已使用的连接，避免重复
                    
                    # 统计O-H键
                    oh_bonds_in_connectivity = [bond for bond in bonds if ((bond['symbol1'] == 'O' and bond['symbol2'] == 'H') or (bond['symbol1'] == 'H' and bond['symbol2'] == 'O'))]
                    print(f"GJF写入: 连接性中有 {len(oh_bonds_in_connectivity)} 个O-H键", file=sys.stderr)
                    for oh_bond in oh_bonds_in_connectivity:
                        print(f"  O-H键详情: 原子{oh_bond['atom1']}-{oh_bond['atom2']}, 距离={oh_bond['distance']:.3f}Å, 键序={oh_bond['bond_order']}", file=sys.stderr)
                    
                    for bond in bonds:
                        atom1 = bond['atom1']
                        atom2 = bond['atom2']
                        bond_order = bond['bond_order']
                        
                        # 创建连接的唯一标识（小的序号在前）
                        conn_key = tuple(sorted([atom1, atom2]))
                        
                        # 避免重复写入同一个连接
                        if conn_key in used_connections:
                            continue
                        used_connections.add(conn_key)
                        
                        # 双向添加连接信息（每个原子都要知道自己的连接）
                        if atom1 not in connectivity_dict:
                            connectivity_dict[atom1] = []
                        connectivity_dict[atom1].append((atom2, bond_order))
                        
                        if atom2 not in connectivity_dict:
                            connectivity_dict[atom2] = []
                        connectivity_dict[atom2].append((atom1, bond_order))
                    
                    # 按GaussView 5标准格式写入连接性信息
                    # 规则：碳原子写入所有连接，氧原子只写O-H键，氢原子空白
                    written_connections = set()  # 记录已写入的连接
                    
                    for atom_idx in range(1, len(atoms) + 1):
                        atom_symbol = atoms[atom_idx-1]['symbol']
                        
                        if atom_idx in connectivity_dict:
                            connections = connectivity_dict[atom_idx]
                            # 按照GaussView标准格式：1个空格开头，原子序号占3位（2个空格+序号），连接信息用1个空格分隔
                            line_parts = [f"{atom_idx:3d}"]  # 原子序号，右对齐占3位

                            # 根据原子类型应用不同规则
                            if atom_symbol == 'C':
                                # 碳原子：写入所有连接到序号更大的原子
                                for conn_atom, bond_order in connections:
                                    if conn_atom > atom_idx:
                                        conn_key = tuple(sorted([atom_idx, conn_atom]))
                                        if conn_key not in written_connections:
                                            # 格式：原子序号 键序，每个占4位（3位原子+1位空格，3位键序）
                                            line_parts.append(f"{conn_atom:3d}{bond_order:4.1f}")
                                            written_connections.add(conn_key)
                            
                            elif atom_symbol == 'O':
                                # 氧原子：根据连接的原子类型决定写入策略
                                is_carbonyl_oxygen = False  # 是否为羰基氧(C=O)
                                has_oh_bond = False          # 是否有O-H键

                                for conn_atom, bond_order in connections:
                                    conn_atom_symbol = atoms[conn_atom-1]['symbol']

                                    # 检查是否连接到碳且为双键（羰基氧）
                                    if conn_atom_symbol == 'C' and bond_order >= 2.0:
                                        is_carbonyl_oxygen = True

                                    # 检查是否有O-H键
                                    if conn_atom_symbol == 'H':
                                        has_oh_bond = True

                                # 根据氧原子类型决定写入策略
                                if is_carbonyl_oxygen:
                                    # 羰基氧：不写入任何连接（在碳原子行已写入C=O双键）
                                    # 这是GaussView的标准格式
                                    pass
                                elif has_oh_bond:
                                    # 羟基氧：强制写入O-H键，不受written_connections限制
                                    for conn_atom, bond_order in connections:
                                        conn_atom_symbol = atoms[conn_atom-1]['symbol']
                                        if conn_atom_symbol == 'H':
                                            # 对于O-H键，不检查written_connections，直接写入
                                            line_parts.append(f"{conn_atom:3d}{bond_order:4.1f}")
                                            # 仍然记录到written_connections以保持一致性
                                            conn_key = tuple(sorted([atom_idx, conn_atom]))
                                            written_connections.add(conn_key)
                                            print(f"  羟基氧{atom_idx}写入O-H键到氢原子{conn_atom}", file=sys.stderr)
                                else:
                                    # 其他氧原子（如醚键氧）：写入所有连接
                                    for conn_atom, bond_order in connections:
                                        conn_key = tuple(sorted([atom_idx, conn_atom]))
                                        if conn_key not in written_connections:
                                            line_parts.append(f"{conn_atom:3d}{bond_order:4.1f}")
                                            written_connections.add(conn_key)
                                            print(f"  氧原子{atom_idx}写入连接到{conn_atom}({atoms[conn_atom-1]['symbol']}), 键序={bond_order}", file=sys.stderr)

                            elif atom_symbol == 'H':
                                # 氢原子：不写入连接（已在对面原子行写入）
                                print(f"  氢原子{atom_idx}的所有连接: {connections}", file=sys.stderr)
                                # 保持空白，不写入任何连接

                            # 如果没有写入任何连接，只写原子序号
                            if len(line_parts) == 1:  # 只有原子序号
                                f.write(f"{line_parts[0]}\n")
                            else:
                                line = ''.join(line_parts)
                                print(f"  写入连接性行: {line}", file=sys.stderr)
                                f.write(line + "\n")
                        else:
                            # 没有连接的原子只写原子序号
                            print(f"  写入连接性行:  {atom_idx:3d}", file=sys.stderr)
                            f.write(f"{atom_idx:3d}\n")
                
                f.write("\n")
            
            messagebox.showinfo("成功", f"GaussView兼容文件已成功导出至:\n{file_path}\n\n注：芳香环将以实线+虚线形式显示")

        except Exception as e:
            messagebox.showerror("导出错误", f"导出GaussView文件时发生错误: {e}")

    def export_multiple_gjf(self) -> None:
        """导出多个GJF文件到指定目录，包含连接性信息"""
        directory = filedialog.askdirectory(title="选择保存多个GJF文件的目录")
        
        if not directory:
            return
        
        try:
            exported_files = []
            
            for idx in sorted(self.selected_isomers):
                if 0 <= idx < len(self.results_data):
                    item = self.results_data[idx]
                    atoms = self.parse_coordinates(item.get('Coords', ''))
                    bonds = self.detect_bond_pattern(atoms, item.get('SMILES', ''))
                    
                    filename = f"{format_formula(self.formula_var.get())}_Isomer{idx+1}.gjf"
                    file_path = os.path.join(directory, filename)
                    
                    with open(file_path, 'w', encoding='utf-8') as f:
                        # 标准GJF文件头部
                        base_name = os.path.splitext(filename)[0]
                        f.write(f"%chk={base_name}.chk\n")
                        f.write("%mem=1GB\n")
                        f.write("%nprocshared=4\n")
                        f.write("#p B3LYP/6-31G(d) opt freq\n\n")
                        f.write(f"{item.get('SMILES', '')} - {format_formula(self.formula_var.get())} - Isomer {idx+1}\n\n")
                        f.write("0 1\n")  # 电荷和自旋多重度
                        
                        # 写入原子坐标
                        for atom in atoms:
                            f.write(f"{atom['symbol']:>2}   {atom['x']:>10.6f}   {atom['y']:>10.6f}   {atom['z']:>10.6f}\n")
                        
                        f.write("\n")
                        
                        # 写入连接性信息
                        if bonds:
                            # 按GaussView标准格式写入连接性
                            connectivity_dict = {}
                            used_connections = set()  # 记录已使用的连接，避免重复
                            
                            for bond in bonds:
                                atom1 = bond['atom1']
                                atom2 = bond['atom2']
                                bond_order = bond['bond_order']
                                
                                # 创建连接的唯一标识（小的序号在前）
                                conn_key = tuple(sorted([atom1, atom2]))
                                
                                # 避免重复写入同一个连接
                                if conn_key in used_connections:
                                    continue
                                used_connections.add(conn_key)
                                
                                # 双向添加连接信息
                                if atom1 not in connectivity_dict:
                                    connectivity_dict[atom1] = []
                                connectivity_dict[atom1].append((atom2, bond_order))
                                
                                if atom2 not in connectivity_dict:
                                    connectivity_dict[atom2] = []
                                connectivity_dict[atom2].append((atom1, bond_order))
                            
                            # 按GaussView 5标准格式写入连接性信息
                            written_connections = set()
                            
                            for atom_idx in range(1, len(atoms) + 1):
                                atom_symbol = atoms[atom_idx-1]['symbol']
                                
                                if atom_idx in connectivity_dict:
                                    connections = connectivity_dict[atom_idx]
                                    line = f" {atom_idx:2d}"
                                    
                                    # 根据原子类型应用不同规则
                                    if atom_symbol == 'C':
                                        # 碳原子：写入所有连接到序号更大的原子
                                        for conn_atom, bond_order in connections:
                                            if conn_atom > atom_idx:
                                                conn_key = tuple(sorted([atom_idx, conn_atom]))
                                                if conn_key not in written_connections:
                                                    line += f" {conn_atom:2d} {bond_order:.1f}"
                                                    written_connections.add(conn_key)
                                    
                                    elif atom_symbol == 'O':
                                        # 氧原子：需要区分羰基氧和其他含氧情况
                                        is_carbonyl_oxygen = False  # 是否为羰基氧(C=O)
                                        has_oh_bond = False          # 是否有O-H键
                                        
                                        for conn_atom, bond_order in connections:
                                            conn_atom_symbol = atoms[conn_atom-1]['symbol']
                                            
                                            # 检查是否连接到碳且为双键（羰基氧）
                                            if conn_atom_symbol == 'C' and bond_order >= 2.0:
                                                is_carbonyl_oxygen = True
                                            
                                            # 检查是否有O-H键
                                            if conn_atom_symbol == 'H':
                                                has_oh_bond = True
                                        
                                        # 根据氧原子类型决定写入策略
                                        if is_carbonyl_oxygen:
                                            # 羰基氧：不写入任何连接（在碳原子行已写入C=O双键）
                                            pass
                                        elif has_oh_bond:
                                            # 羟基氧：不写入O-H键（改为在氢原子行写入）
                                            pass
                                        else:
                                            # 其他氧原子（如醚键氧）：写入连接到序号更大的原子
                                            for conn_atom, bond_order in connections:
                                                if conn_atom > atom_idx:
                                                    conn_key = tuple(sorted([atom_idx, conn_atom]))
                                                    if conn_key not in written_connections:
                                                        line += f" {conn_atom:2d} {bond_order:.1f}"
                                                        written_connections.add(conn_key)
                                    
                                    elif atom_symbol == 'H':
                                        # 氢原子：如果连接到氧，写入连接（特别是羧基的O-H键）
                                        for conn_atom, bond_order in connections:
                                            conn_atom_symbol = atoms[conn_atom-1]['symbol']
                                            if conn_atom_symbol == 'O':
                                                conn_key = tuple(sorted([atom_idx, conn_atom]))
                                                if conn_key not in written_connections:
                                                    line += f" {conn_atom:2d} {bond_order:.1f}"
                                                    written_connections.add(conn_key)
                                                break  # 氢原子只写一个连接（通常是H-O或H-C）
                                    
                                    # 如果没有写入任何连接，只写原子序号
                                    if len(line.strip()) == len(f"{atom_idx:2d}"):
                                        f.write(f" {atom_idx:2d}\n")
                                    else:
                                        f.write(line + "\n")
                                else:
                                    # 没有连接的原子只写原子序号
                                    f.write(f" {atom_idx:2d}\n")
                        
                        f.write("\n")
                    
                    exported_files.append(filename)
            
            messagebox.showinfo("成功", f"已成功导出{len(exported_files)}个GaussView兼容文件到:\n{directory}\n\n文件列表:\n" + "\n".join(exported_files))

        except Exception as e:
            messagebox.showerror("导出错误", f"导出GaussView文件时发生错误: {e}")

    def export_results(self) -> None:
        """导出所有异构体（原功能）"""
        if not self.results_data:
            messagebox.showwarning("警告", "没有可导出的数据。请先执行分析。")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="保存同分异构体分析结果"
        )

        if not file_path:
            return
        
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(f"*** MolGen+ 简易同分异构体分析结果 ***\n")
                f.write(f"分子式: {format_formula(self.formula_var.get())}\n")
                f.write(f"文件名: {os.path.basename(file_path)}\n")
                f.write(f"生成时间: {tk.Tcl().eval('clock format [clock seconds]')}\n")
                f.write(f"共生成 {len(self.results_data)} 种同分异构体 (基于优化的模拟结构集)\n\n")

                for i, item in enumerate(self.results_data):
                    f.write(f"=========================================\n")
                    f.write(f"异构体 {i+1}:\n")
                    f.write(f"SMILES: {item.get('SMILES', '')}\n")
                    f.write(f"InChIKey: {item.get('InChIKey', '')}\n")
                    f.write(f"结构类型: {item.get('StructureType', '')}\n")
                    f.write(f"-----------------------------------------\n")
                    f.write(f"原子坐标 ({item.get('Dimension', '')} Angstrom):\n")
                    f.write(item.get('Coords', '') + "\n")
                    f.write("=========================================\n\n")
            
            messagebox.showinfo("成功", f"结果已成功导出至:\n{file_path}")

        except Exception as e:
            messagebox.showerror("导出错误", f"导出文件时发生错误: {e}")

    def export_xyz_format(self) -> None:
        """导出XYZ格式，兼容GaussView和其他分子可视化软件"""
        if not self.results_data:
            messagebox.showwarning("警告", "没有可导出的数据。请先执行分析。")
            return
        
        if not self.selected_isomers:
            messagebox.showwarning("警告", "请先选择要导出的异构体。")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".xyz",
            filetypes=[("XYZ files", "*.xyz"), ("All files", "*.*")],
            title="保存XYZ格式文件"
        )

        if not file_path:
            return
        
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                for idx in sorted(self.selected_isomers):
                    if 0 <= idx < len(self.results_data):
                        item = self.results_data[idx]
                        atoms = self.parse_coordinates(item.get('Coords', ''))
                        
                        if not atoms:
                            continue
                        
                        # 写入XYZ格式头部
                        f.write(f"{len(atoms)}\n")
                        f.write(f"Isomer {idx+1} - {item.get('SMILES', '')} - {format_formula(self.formula_var.get())}\n")
                        
                        # 写入原子坐标
                        for atom in atoms:
                            f.write(f"  {atom['symbol']:>2}   {atom['x']:>12.6f}   {atom['y']:>12.6f}   {atom['z']:>12.6f}\n")
                        
                        # 在多个分子之间添加空行分隔
                        if idx != sorted(self.selected_isomers)[-1]:
                            f.write("\n")
            
            messagebox.showinfo("成功", f"XYZ格式文件已成功导出至:\n{file_path}")

        except Exception as e:
            messagebox.showerror("导出错误", f"导出XYZ文件时发生错误: {e}")


# ==============================================================================
# 机器学习分析窗口
# ==============================================================================

class MLAnalysisWindow:
    """机器学习分析专用窗口"""

    def __init__(self, parent: tk.Tk, smiles_list: list):
        self.window = tk.Toplevel(parent)
        self.window.title("🤖 机器学习分子分析")
        self.window.geometry("1000x700")
        self.window.transient(parent)
        self.window.grab_set()

        self.smiles_list = smiles_list
        self.analyzer = MLIsomerAnalyzer()
        self.current_report = {}

        # 创建界面
        self._create_widgets()

        # 自动运行初始分析
        self._run_initial_analysis()

    def _create_widgets(self):
        """创建界面组件"""
        # 工具栏
        toolbar = tk.Frame(self.window, padx=10, pady=10, bg="#f0f0f0")
        toolbar.pack(fill=tk.X)

        tk.Label(toolbar, text="分析类型:", font=('Arial', 10, 'bold'), bg="#f0f0f0").pack(side=tk.LEFT, padx=5)

        tk.Button(toolbar, text="📊 综合报告", command=self._run_full_report,
                 bg="#4CAF50", fg="white", width=12).pack(side=tk.LEFT, padx=5)
        tk.Button(toolbar, text="🎯 聚类分析", command=self._run_clustering,
                 bg="#2196F3", fg="white", width=12).pack(side=tk.LEFT, padx=5)
        tk.Button(toolbar, text="📈 PCA分析", command=self._run_pca,
                 bg="#FF9800", fg="white", width=12).pack(side=tk.LEFT, padx=5)
        tk.Button(toolbar, text="⚠️ 异常检测", command=self._run_outlier_detection,
                 bg="#F44336", fg="white", width=12).pack(side=tk.LEFT, padx=5)
        tk.Button(toolbar, text="🎨 可视化", command=self._run_visualization,
                 bg="#9C27B0", fg="white", width=12).pack(side=tk.LEFT, padx=5)
        tk.Button(toolbar, text="💾 保存报告", command=self._save_report,
                 bg="#607D8B", fg="white", width=12).pack(side=tk.LEFT, padx=5)

        # 主内容区域
        content_paned = tk.PanedWindow(self.window, orient=tk.HORIZONTAL)
        content_paned.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # 左侧：结果文本
        left_frame = tk.Frame(content_paned)
        right_frame = tk.Frame(content_paned)

        content_paned.add(left_frame, minsize=400)
        content_paned.add(right_frame, minsize=300)

        # 结果文本区域
        tk.Label(left_frame, text="📋 分析结果", font=('Arial', 10, 'bold')).pack(anchor=tk.W, padx=5, pady=5)

        self.result_text = scrolledtext.ScrolledText(left_frame, wrap=tk.WORD,
                                                     font=('Consolas', 9), height=30)
        self.result_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # 右侧：图表和信息
        tk.Label(right_frame, text="📈 图表预览", font=('Arial', 10, 'bold')).pack(anchor=tk.W, padx=5, pady=5)

        self.chart_canvas = tk.Canvas(right_frame, bg='white', height=300)
        self.chart_canvas.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # 统计信息
        tk.Label(right_frame, text="📊 快速统计", font=('Arial', 10, 'bold')).pack(anchor=tk.W, padx=5, pady=5)

        self.stats_text = scrolledtext.ScrolledText(right_frame, wrap=tk.WORD,
                                                   font=('Consolas', 8), height=10)
        self.stats_text.pack(fill=tk.X, padx=5, pady=5)

    def _run_initial_analysis(self):
        """运行初始分析"""
        self.result_text.insert(tk.END, f"🤖 机器学习分子分析工具\n")
        self.result_text.insert(tk.END, f"{'='*60}\n\n")
        self.result_text.insert(tk.END, f"📝 输入: {len(self.smiles_list)} 个分子\n")
        self.result_text.insert(tk.END, f"⏳ 正在运行初始分析...\n\n")

        self.window.update()

        # 运行综合报告
        self._run_full_report()

    def _run_full_report(self):
        """生成综合报告"""
        try:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, "⏳ 正在生成综合分析报告...\n\n")
            self.window.update()

            report = self.analyzer.generate_analysis_report(self.smiles_list)
            self.current_report = report

            self._display_report(report)
            self._display_stats(report)

        except Exception as e:
            self.result_text.insert(tk.END, f"❌ 错误: {str(e)}\n")

    def _run_clustering(self):
        """运行聚类分析"""
        try:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, "⏳ 正在进行聚类分析...\n\n")
            self.window.update()

            result = self.analyzer.perform_clustering(self.smiles_list,
                                                    method='kmeans',
                                                    n_clusters=3)

            self._display_clustering_result(result)
            self._update_stats('clustering', result)

        except Exception as e:
            self.result_text.insert(tk.END, f"❌ 错误: {str(e)}\n")

    def _run_pca(self):
        """运行PCA分析"""
        try:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, "⏳ 正在进行PCA降维分析...\n\n")
            self.window.update()

            result = self.analyzer.perform_pca_analysis(self.smiles_list,
                                                       n_components=2)

            self._display_pca_result(result)
            self._update_stats('pca', result)

        except Exception as e:
            self.result_text.insert(tk.END, f"❌ 错误: {str(e)}\n")

    def _run_outlier_detection(self):
        """运行异常检测"""
        try:
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, "⏳ 正在检测异常分子...\n\n")
            self.window.update()

            result = self.analyzer.detect_outliers(self.smiles_list,
                                                  contamination=0.1)

            self._display_outlier_result(result)
            self._update_stats('outliers', result)

        except Exception as e:
            self.result_text.insert(tk.END, f"❌ 错误: {str(e)}\n")

    def _run_visualization(self):
        """运行可视化"""
        try:
            self.result_text.insert(tk.END, "⏳ 正在生成可视化图表...\n\n")
            self.window.update()

            viz_path = self.analyzer.visualize_clusters(self.smiles_list,
                                                       method='kmeans',
                                                       n_clusters=3)

            if viz_path:
                self.result_text.insert(tk.END, f"✅ 可视化已保存到:\n{viz_path}\n")
                # 尝试显示图片
                try:
                    from PIL import Image, ImageTk
                    img = Image.open(viz_path)
                    img.thumbnail((400, 300))
                    photo = ImageTk.PhotoImage(img)
                    self.chart_canvas.create_image(10, 10, anchor=tk.NW, image=photo)
                    self.chart_canvas.image = photo
                except ImportError:
                    self.result_text.insert(tk.END, "⚠️ 需要Pillow库来显示图片\n")
            else:
                self.result_text.insert(tk.END, "❌ 可视化生成失败\n")

        except Exception as e:
            self.result_text.insert(tk.END, f"❌ 错误: {str(e)}\n")

    def _save_report(self):
        """保存报告到文件"""
        try:
            file_path = filedialog.asksaveasfilename(
                defaultextension=".txt",
                filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
                title="保存ML分析报告"
            )

            if not file_path:
                return

            with open(file_path, 'w', encoding='utf-8') as f:
                f.write("🤖 机器学习分子分析报告\n")
                f.write(f"{'='*60}\n\n")
                f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"分析分子数: {len(self.smiles_list)}\n\n")
                f.write(self.result_text.get(1.0, tk.END))

            messagebox.showinfo("成功", f"报告已保存到:\n{file_path}")

        except Exception as e:
            messagebox.showerror("错误", f"保存失败: {str(e)}")

    def _display_report(self, report: dict):
        """显示完整报告"""
        # 摘要
        summary = report.get('summary', {})
        self.result_text.insert(tk.END, "📊 分析摘要\n")
        self.result_text.insert(tk.END, f"{'-'*40}\n")
        self.result_text.insert(tk.END, f"总输入: {summary.get('total_input', 0)}\n")
        self.result_text.insert(tk.END, f"有效分子: {summary.get('valid_molecules', 0)}\n")
        self.result_text.insert(tk.END, f"无效分子: {summary.get('invalid_molecules', 0)}\n\n")

        # 聚类
        clustering = report.get('clustering', {})
        if clustering:
            self.result_text.insert(tk.END, "🎯 聚类分析\n")
            self.result_text.insert(tk.END, f"{'-'*40}\n")
            self.result_text.insert(tk.END, f"方法: {clustering.get('method', 'N/A')}\n")
            self.result_text.insert(tk.END, f"聚类数: {clustering.get('n_clusters', 0)}\n")
            if clustering.get('silhouette_score'):
                self.result_text.insert(tk.END,
                                     f"轮廓系数: {clustering['silhouette_score']:.4f}\n")
            self.result_text.insert(tk.END, "\n聚类分布:\n")
            for name, count in clustering.get('cluster_sizes', {}).items():
                self.result_text.insert(tk.END, f"  {name}: {count} 个分子\n")
            self.result_text.insert(tk.END, "\n")

        # 多样性
        diversity = report.get('diversity', {})
        if diversity:
            self.result_text.insert(tk.END, "🌈 多样性分析\n")
            self.result_text.insert(tk.END, f"{'-'*40}\n")
            self.result_text.insert(tk.END,
                                 f"多样性得分: {diversity.get('diversity_score', 0):.4f}\n")
            self.result_text.insert(tk.END,
                                 f"平均成对距离: {diversity.get('mean_pairwise_distance', 0):.4f}\n\n")

        # 异常
        outliers = report.get('outliers', {})
        if outliers:
            self.result_text.insert(tk.END, "⚠️  异常检测\n")
            self.result_text.insert(tk.END, f"{'-'*40}\n")
            self.result_text.insert(tk.END, f"异常分子数: {outliers.get('n_outliers', 0)}\n")
            self.result_text.insert(tk.END,
                                 f"异常比例: {outliers.get('outlier_ratio', 0):.2%}\n\n")

    def _display_clustering_result(self, result: dict):
        """显示聚类结果"""
        self.result_text.insert(tk.END, "🎯 聚类分析结果\n")
        self.result_text.insert(tk.END, f"{'='*50}\n\n")

        if 'error' in result:
            self.result_text.insert(tk.END, f"❌ {result['error']}\n")
            return

        self.result_text.insert(tk.END, f"聚类方法: {result.get('method', 'N/A')}\n")
        self.result_text.insert(tk.END, f"聚类数量: {result.get('n_clusters', 0)}\n\n")

        if result.get('silhouette_score'):
            self.result_text.insert(tk.END,
                                 f"轮廓系数: {result['silhouette_score']:.4f}\n")

        self.result_text.insert(tk.END, "\n聚类详情:\n")
        for cluster_name, molecules in result.get('clusters', {}).items():
            self.result_text.insert(tk.END, f"\n{cluster_name} ({len(molecules)} 个分子):\n")
            for i, smiles in enumerate(molecules[:10]):  # 只显示前10个
                self.result_text.insert(tk.END, f"  {i+1}. {smiles}\n")
            if len(molecules) > 10:
                self.result_text.insert(tk.END, f"  ... 还有 {len(molecules)-10} 个\n")

    def _display_pca_result(self, result: dict):
        """显示PCA结果"""
        self.result_text.insert(tk.END, "📈 PCA降维分析\n")
        self.result_text.insert(tk.END, f"{'='*50}\n\n")

        if 'error' in result:
            self.result_text.insert(tk.END, f"❌ {result['error']}\n")
            return

        explained_var = result.get('explained_variance_ratio', [])
        cumulative_var = result.get('cumulative_variance', [])

        self.result_text.insert(tk.END, "主成分解释方差比:\n")
        for i, (exp, cum) in enumerate(zip(explained_var, cumulative_var)):
            self.result_text.insert(tk.END,
                                 f"  PC{i+1}: {exp:.2%} (累积: {cum:.2%})\n")

        importance = result.get('feature_importance', {})
        if importance:
            self.result_text.insert(tk.END, "\n特征重要性 (前5名):\n")
            for pc_name, top_features in importance.items():
                self.result_text.insert(tk.END, f"\n{pc_name}:\n")
                for feature, value in list(top_features.items())[:5]:
                    self.result_text.insert(tk.END,
                                         f"  {feature}: {value:.4f}\n")

    def _display_outlier_result(self, result: dict):
        """显示异常检测结果"""
        self.result_text.insert(tk.END, "⚠️  异常检测结果\n")
        self.result_text.insert(tk.END, f"{'='*50}\n\n")

        if 'error' in result:
            self.result_text.insert(tk.END, f"❌ {result['error']}\n")
            return

        outliers = result.get('outliers', [])
        normal = result.get('normal', [])

        self.result_text.insert(tk.END, f"异常分子数: {len(outliers)}\n")
        self.result_text.insert(tk.END, f"正常分子数: {len(normal)}\n")
        self.result_text.insert(tk.END, f"异常比例: {result.get('outlier_ratio', 0):.2%}\n\n")

        if outliers:
            self.result_text.insert(tk.END, "异常分子列表:\n")
            for i, smiles in enumerate(outliers[:10]):
                self.result_text.insert(tk.END, f"  {i+1}. {smiles}\n")
            if len(outliers) > 10:
                self.result_text.insert(tk.END, f"  ... 还有 {len(outliers)-10} 个\n")

    def _display_stats(self, report: dict):
        """在右侧显示统计信息"""
        self.stats_text.delete(1.0, tk.END)

        self.stats_text.insert(tk.END, "📊 快速统计\n")
        self.stats_text.insert(tk.END, f"{'='*30}\n\n")

        # 摘要统计
        summary = report.get('summary', {})
        self.stats_text.insert(tk.END, f"分子总数: {summary.get('valid_molecules', 0)}\n\n")

        # 聚类统计
        clustering = report.get('clustering', {})
        if clustering:
            self.stats_text.insert(tk.END, f"聚类数: {clustering.get('n_clusters', 0)}\n")
            if clustering.get('silhouette_score'):
                self.stats_text.insert(tk.END,
                                     f"轮廓系数: {clustering['silhouette_score']:.3f}\n")
            self.stats_text.insert(tk.END, "\n")

        # 多样性统计
        diversity = report.get('diversity', {})
        if diversity:
            self.stats_text.insert(tk.END,
                                 f"多样性得分: {diversity.get('diversity_score', 0):.3f}\n")
            self.stats_text.insert(tk.END,
                                 f"平均距离: {diversity.get('mean_pairwise_distance', 0):.3f}\n\n")

        # 异常统计
        outliers = report.get('outliers', {})
        if outliers:
            self.stats_text.insert(tk.END, f"异常分子: {outliers.get('n_outliers', 0)}\n")
            self.stats_text.insert(tk.END,
                                 f"异常比例: {outliers.get('outlier_ratio', 0):.1%}\n\n")

    def _update_stats(self, analysis_type: str, result: dict):
        """更新统计信息"""
        if analysis_type == 'clustering':
            self.stats_text.delete(1.0, tk.END)
            self.stats_text.insert(tk.END, f"聚类数: {result.get('n_clusters', 0)}\n")
            if result.get('silhouette_score'):
                self.stats_text.insert(tk.END,
                                     f"轮廓系数: {result['silhouette_score']:.3f}\n")
        elif analysis_type == 'pca':
            explained_var = result.get('explained_variance_ratio', [])
            if explained_var:
                self.stats_text.delete(1.0, tk.END)
                self.stats_text.insert(tk.END, "解释方差比:\n")
                for i, var in enumerate(explained_var):
                    self.stats_text.insert(tk.END, f"PC{i+1}: {var:.2%}\n")
        elif analysis_type == 'outliers':
            self.stats_text.delete(1.0, tk.END)
            self.stats_text.insert(tk.END, f"异常分子: {len(result.get('outliers', []))}\n")
            self.stats_text.insert(tk.END,
                                 f"异常比例: {result.get('outlier_ratio', 0):.1%}\n")


if __name__ == "__main__":
    root = tk.Tk()
    app = IsomerAnalyzerGUI(root)
    root.mainloop()