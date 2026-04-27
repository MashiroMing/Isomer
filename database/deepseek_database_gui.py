# deepseek_database_gui.py
# DeepSeek数据库集成GUI界面

import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext, filedialog
import threading
import json
import os
import sys
from datetime import datetime

# 添加项目根目录到 Python 路径（支持作为主程序运行）
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

# 导入 DeepSeekDatabaseIntegration（根据运行方式选择）
try:
    # 尝试相对导入（作为模块导入时）
    from .deepseek_database_integration import DeepSeekDatabaseIntegration
except ImportError:
    # 回退到绝对导入（作为主程序运行时）
    from database.deepseek_database_integration import DeepSeekDatabaseIntegration

# 检查依赖
try:
    import rdkit
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    from openai import OpenAI
    HAS_OPENAI = True
except ImportError:
    HAS_OPENAI = False

class DeepSeekDatabaseGUI:
    """DeepSeek数据库集成GUI类"""
    
    def __init__(self, master):
        self.master = master
        self.master.title("DeepSeek AI 数据库集成")
        self.master.geometry("1000x700")

        # 检查依赖
        if not HAS_OPENAI:
            messagebox.showwarning(
                "缺少依赖",
                "提示: DeepSeek AI 模块不可用，请安装 openai:\n"
                "  pip install openai\n"
                "或运行安装脚本:\n"
                "  install_dependencies.bat"
            )

        if not HAS_RDKIT:
            messagebox.showwarning(
                "缺少依赖",
                "提示: RDKit 模块不可用，请安装 rdkit:\n"
                "  conda install -c conda-forge rdkit\n"
                "或运行安装脚本:\n"
                "  install_dependencies.bat"
            )

        # 集成实例
        self.integration = None
        self.current_library_stats = None

        # 创建界面
        self.create_widgets()

        # 默认API密钥
        self.api_key_var.set(os.environ.get('DEEPSEEK_API_KEY', ''))
    
    def create_widgets(self):
        """创建界面组件"""
        # 主框架
        main_frame = ttk.Frame(self.master, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # 配置行列权重
        self.master.columnconfigure(0, weight=1)
        self.master.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(3, weight=1)
        
        # 标题
        title_label = ttk.Label(main_frame, text="DeepSeek AI 数据库集成", 
                               font=('Arial', 16, 'bold'))
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 20))
        
        # 连接区域
        connection_frame = ttk.LabelFrame(main_frame, text="连接配置", padding="10")
        connection_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))
        
        self.api_key_var = tk.StringVar()
        ttk.Label(connection_frame, text="API密钥:").grid(row=0, column=0, sticky=tk.W)
        api_key_entry = ttk.Entry(connection_frame, textvariable=self.api_key_var, 
                                  width=50, show="*")
        api_key_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=10)
        
        # 默认使用 data 目录中的文件
        default_library_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'data',
            'molecule_library.json'
        )
        self.library_path_var = tk.StringVar(value=default_library_path)
        ttk.Label(connection_frame, text="数据库路径:").grid(row=1, column=0, sticky=tk.W, pady=5)
        ttk.Entry(connection_frame, textvariable=self.library_path_var, 
                 width=50).grid(row=1, column=1, sticky=(tk.W, tk.E), padx=10, pady=5)
        
        browse_button = ttk.Button(connection_frame, text="浏览...", 
                                   command=self.browse_library)
        browse_button.grid(row=1, column=2, padx=5, pady=5)
        
        connect_button = ttk.Button(connection_frame, text="连接", 
                                    command=self.connect_to_database)
        connect_button.grid(row=0, column=2, padx=5)
        
        # 分子库统计
        stats_frame = ttk.LabelFrame(main_frame, text="分子库统计", padding="10")
        stats_frame.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))
        
        self.stats_text = scrolledtext.ScrolledText(stats_frame, height=5, width=80)
        self.stats_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        stats_frame.columnconfigure(0, weight=1)
        
        # 选项卡
        notebook = ttk.Notebook(main_frame)
        notebook.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=(0, 10))
        
        # AI增强选项卡
        enhance_frame = ttk.Frame(notebook, padding="10")
        notebook.add(enhance_frame, text="AI增强")
        self.create_enhance_tab(enhance_frame)
        
        # 备份管理选项卡
        backup_frame = ttk.Frame(notebook, padding="10")
        notebook.add(backup_frame, text="备份管理")
        self.create_backup_tab(backup_frame)
        
        # 状态栏
        self.status_var = tk.StringVar(value="就绪")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, 
                              relief=tk.SUNKEN, anchor=tk.W)
        status_bar.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E))
    
    def create_enhance_tab(self, parent):
        """创建AI增强选项卡"""
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(2, weight=1)
        
        # 分子式输入
        formula_frame = ttk.Frame(parent)
        formula_frame.grid(row=0, column=0, sticky=(tk.W, tk.E), pady=(0, 10))
        
        self.formula_var = tk.StringVar(value="C10H22")
        ttk.Label(formula_frame, text="分子式:").grid(row=0, column=0)
        ttk.Entry(formula_frame, textvariable=self.formula_var, width=20).grid(row=0, column=1, padx=5)
        
        self.target_count_var = tk.StringVar(value="30")
        ttk.Label(formula_frame, text="目标数量:").grid(row=0, column=2, padx=(20, 0))
        ttk.Entry(formula_frame, textvariable=self.target_count_var, width=10).grid(row=0, column=3, padx=5)
        
        # 快捷按钮
        quick_frame = ttk.Frame(parent)
        quick_frame.grid(row=1, column=0, sticky=tk.W, pady=(0, 10))
        
        ttk.Button(quick_frame, text="C8H18 → 20个",
                  command=lambda: self.quick_enhance("C8H18", 20)).grid(row=0, column=0, padx=5)
        ttk.Button(quick_frame, text="C9H20 → 25个",
                  command=lambda: self.quick_enhance("C9H20", 25)).grid(row=0, column=1, padx=5)
        ttk.Button(quick_frame, text="C10H22 → 30个",
                  command=lambda: self.quick_enhance("C10H22", 30)).grid(row=0, column=2, padx=5)
        ttk.Button(quick_frame, text="全部更新",
                  command=self.enhance_all).grid(row=0, column=3, padx=5)
        
        # 自定义增强
        ttk.Button(parent, text="增强选定的分子式",
                  command=self.enhance_selected).grid(row=0, column=1, padx=10)
        
        # 结果显示
        result_frame = ttk.LabelFrame(parent, text="增强结果", padding="10")
        result_frame.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S))
        result_frame.columnconfigure(0, weight=1)
        result_frame.rowconfigure(0, weight=1)
        
        self.enhance_result = scrolledtext.ScrolledText(result_frame, height=15, width=80)
        self.enhance_result.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
    
    def create_backup_tab(self, parent):
        """创建备份管理选项卡"""
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(1, weight=1)
        
        # 备份列表
        list_frame = ttk.Frame(parent)
        list_frame.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(0, 10))
        
        ttk.Button(list_frame, text="刷新备份列表",
                  command=self.refresh_backup_list).grid(row=0, column=0, padx=5)
        ttk.Button(list_frame, text="查看详情",
                  command=self.view_backup_detail).grid(row=0, column=1, padx=5)
        ttk.Button(list_frame, text="回滚到选定的备份",
                  command=self.rollback_backup).grid(row=0, column=2, padx=5)
        
        # 备份列表
        columns = ("filename", "size", "created")
        self.backup_tree = ttk.Treeview(parent, columns=columns, show="headings", height=10)
        self.backup_tree.heading("filename", text="文件名")
        self.backup_tree.heading("size", text="大小")
        self.backup_tree.heading("created", text="创建时间")
        self.backup_tree.column("filename", width=200)
        self.backup_tree.column("size", width=100)
        self.backup_tree.column("created", width=200)
        self.backup_tree.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E))
        
        # 滚动条
        scrollbar = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=self.backup_tree.yview)
        scrollbar.grid(row=1, column=2, sticky=(tk.N, tk.S))
        self.backup_tree.configure(yscrollcommand=scrollbar.set)
    
    def browse_library(self):
        """浏览数据库文件"""
        filename = filedialog.askopenfilename(
            title="选择分子库文件",
            filetypes=[("JSON文件", "*.json"), ("所有文件", "*.*")]
        )
        if filename:
            self.library_path_var.set(filename)
    
    def connect_to_database(self):
        """连接到数据库"""
        try:
            api_key = self.api_key_var.get().strip()
            library_path = self.library_path_var.get().strip()
            
            if not api_key:
                messagebox.showerror("错误", "请输入API密钥")
                return
            
            self.status_var.set("正在连接...")
            
            # 创建集成实例
            self.integration = DeepSeekDatabaseIntegration(
                api_key=api_key,
                library_file=library_path
            )
            
            # 获取统计信息
            self.current_library_stats = self.integration.get_library_statistics()
            
            # 更新显示
            self.update_stats_display()
            self.refresh_backup_list()
            
            self.status_var.set("连接成功")
            messagebox.showinfo("成功", "数据库连接成功！")
            
        except Exception as e:
            self.status_var.set("连接失败")
            messagebox.showerror("错误", f"连接失败: {str(e)}")
    
    def update_stats_display(self):
        """更新统计信息显示"""
        if not self.current_library_stats:
            return
        
        self.stats_text.delete(1.0, tk.END)
        
        stats = self.current_library_stats
        text = f"总分子数: {stats['total_molecules']}\n\n"
        
        text += "按类别统计:\n"
        for category, count in stats['by_category'].items():
            text += f"  {category}: {count}\n"
        
        text += "\n按分子式统计:\n"
        for formula, count in sorted(stats['by_formula'].items()):
            text += f"  {formula}: {count}\n"
        
        text += "\n按来源统计:\n"
        for source, count in sorted(stats['by_source'].items()):
            text += f"  {source}: {count}\n"
        
        self.stats_text.insert(1.0, text)
    
    def quick_enhance(self, formula, target_count):
        """快速增强指定分子式"""
        if not self.integration:
            messagebox.showerror("错误", "请先连接数据库")
            return
        
        self.formula_var.set(formula)
        self.target_count_var.set(str(target_count))
        self.enhance_selected()
    
    def enhance_selected(self):
        """增强选定的分子式"""
        if not self.integration:
            messagebox.showerror("错误", "请先连接数据库")
            return
        
        formula = self.formula_var.get().strip()
        try:
            target_count = int(self.target_count_var.get().strip())
        except ValueError:
            messagebox.showerror("错误", "目标数量必须是整数")
            return
        
        if not formula:
            messagebox.showerror("错误", "请输入分子式")
            return
        
        # 在新线程中执行
        def run_enhance():
            try:
                self.status_var.set(f"正在增强 {formula}...")
                result = self.integration.ai_enhance_library(
                    target_formula=formula,
                    target_total=target_count
                )
                
                # 更新结果显示
                self.master.after(0, lambda: self.update_enhance_result(result))
                
                # 更新统计
                self.current_library_stats = self.integration.get_library_statistics()
                self.master.after(0, self.update_stats_display)
                
                self.master.after(0, lambda: self.status_var.set("增强完成"))
                
            except Exception as e:
                error_msg = f"增强失败: {str(e)}"
                self.master.after(0, lambda: messagebox.showerror("错误", error_msg))
                self.master.after(0, lambda: self.status_var.set("增强失败"))
        
        threading.Thread(target=run_enhance, daemon=True).start()
    
    def enhance_all(self):
        """增强所有分子式"""
        if not self.integration:
            messagebox.showerror("错误", "请先连接数据库")
            return
        
        if not messagebox.askyesno("确认", "是否要增强所有主要分子式？\n\n这将更新: C8H18, C9H20, C10H22"):
            return
        
        def run_enhance_all():
            try:
                formulas = [("C8H18", 20), ("C9H20", 25), ("C10H22", 30)]
                all_results = {}
                
                for formula, target in formulas:
                    self.master.after(0, lambda f=formula: self.status_var.set(f"正在增强 {f}..."))
                    result = self.integration.ai_enhance_library(
                        target_formula=formula,
                        target_total=target
                    )
                    all_results[formula] = result
                
                # 显示汇总结果
                self.master.after(0, lambda: self.show_summary_result(all_results))
                
                # 更新统计
                self.current_library_stats = self.integration.get_library_statistics()
                self.master.after(0, self.update_stats_display)
                
                self.master.after(0, lambda: self.status_var.set("全部增强完成"))
                
            except Exception as e:
                error_msg = f"批量增强失败: {str(e)}"
                self.master.after(0, lambda: messagebox.showerror("错误", error_msg))
                self.master.after(0, lambda: self.status_var.set("批量增强失败"))
        
        threading.Thread(target=run_enhance_all, daemon=True).start()
    
    def update_enhance_result(self, result):
        """更新增强结果显示"""
        self.enhance_result.delete(1.0, tk.END)
        
        text = json.dumps(result, ensure_ascii=False, indent=2)
        self.enhance_result.insert(1.0, text)
    
    def show_summary_result(self, results):
        """显示汇总结果"""
        summary = {"总结果": "成功", "详情": results}
        self.update_enhance_result(summary)
        
        total_added = sum(r.get("added_count", 0) for r in results.values() if r.get("success"))
        messagebox.showinfo("完成", f"批量增强完成！\n\n总共添加了 {total_added} 个新分子。")
    
    def refresh_backup_list(self):
        """刷新备份列表"""
        if not self.integration:
            return
        
        # 清空列表
        for item in self.backup_tree.get_children():
            self.backup_tree.delete(item)
        
        # 获取备份列表
        backups = self.integration.get_backup_list()
        
        # 添加到列表
        for backup in backups:
            self.backup_tree.insert("", tk.END, values=(
                backup["filename"],
                f"{backup['size']} bytes",
                backup["created"]
            ))
    
    def view_backup_detail(self):
        """查看备份详情"""
        selection = self.backup_tree.selection()
        if not selection:
            messagebox.showwarning("提示", "请选择一个备份文件")
            return
        
        item = self.backup_tree.item(selection[0])
        filename = item["values"][0]
        
        try:
            backup_path = os.path.join(self.integration.backup_dir, filename)
            with open(backup_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # 显示详情
            detail_window = tk.Toplevel(self.master)
            detail_window.title(f"备份详情: {filename}")
            detail_window.geometry("600x400")
            
            text = scrolledtext.ScrolledText(detail_window, width=70, height=20)
            text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
            text.insert(1.0, json.dumps(data, ensure_ascii=False, indent=2))
            text.config(state=tk.DISABLED)
            
        except Exception as e:
            messagebox.showerror("错误", f"读取备份失败: {str(e)}")
    
    def rollback_backup(self):
        """回滚到选定的备份"""
        if not self.integration:
            messagebox.showerror("错误", "请先连接数据库")
            return
        
        selection = self.backup_tree.selection()
        if not selection:
            messagebox.showwarning("提示", "请选择一个备份文件")
            return
        
        item = self.backup_tree.item(selection[0])
        filename = item["values"][0]
        
        if not messagebox.askyesno("确认", f"确定要回滚到备份 {filename} 吗？\n\n当前数据库将被备份！"):
            return
        
        try:
            result = self.integration.rollback_operation(filename)
            
            if result["success"]:
                messagebox.showinfo("成功", "回滚成功！\n\n" + 
                                  f"新备份: {result.get('new_backup', 'N/A')}")
                
                # 刷新显示
                self.current_library_stats = self.integration.get_library_statistics()
                self.update_stats_display()
                self.refresh_backup_list()
            else:
                messagebox.showerror("失败", f"回滚失败: {result.get('error', 'Unknown')}")
                
        except Exception as e:
            messagebox.showerror("错误", f"回滚失败: {str(e)}")

def main():
    """主函数"""
    root = tk.Tk()
    app = DeepSeekDatabaseGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
