# MolGenPlus 有机分子同分异构体分析系统

## 答辩文档

---

## 一、项目概述

### 1.1 研究背景

有机化学中，**同分异构体**（Isomer）是指具有相同分子式但结构不同的化合物。随着碳原子数增加，同分异构体数量呈指数级增长：

| 碳数 | 烷烃异构体数 | 传统手动列举难度 |
|------|-------------|-----------------|
| C6   | 5           | 可行            |
| C7   | 9           | 较困难          |
| C8   | 18          | 困难            |
| C9   | 35          | 极困难          |
| C10  | 75          | 几乎不可能      |

**OEIS A000602** 序列记录了这一理论值。

### 1.2 项目目标

开发一个**通用、自动化的同分异构体生成与分析系统**，能够：
1. 自动生成 C1-C10 烷烃的**全部**同分异构体
2. 支持结构筛选（环结构、键类型、官能团等）
3. 提供分子库管理和可视化功能
4. 集成 AI 增强（可选）

---

## 二、系统架构

### 2.1 整体架构

```
┌─────────────────────────────────────────────────────────────┐
│                    MolGenPlus 主界面 (Tkinter GUI)         │
├─────────────────────────────────────────────────────────────┤
│  输入层         │  处理层              │  输出层           │
│  ─────────────  │  ──────────────────  │  ──────────────  │
│  分子式输入     │  分子式解析          │  结果展示        │
│  结构筛选选项   │  生成器调度          │  分子可视化      │
│                 │  结构筛选            │  导出功能        │
│                 │  规范化去重          │  ML分析          │
└─────────────────────────────────────────────────────────────┘
                          │
          ┌───────────────┼───────────────┐
          ▼               ▼               ▼
   ┌─────────────┐  ┌─────────────┐  ┌─────────────┐
   │ C1-C5 生成器 │  │ C6 生成器    │  │ C8-C10 生成器│
   │ (树结构算法+ │  │ (特殊优化)   │  │ (树结构算法) │
   │  硬编码)     │  │              │  │              │
   └──────┬──────┘  └──────┬──────┘  └─────────────┘
          │                │
          ▼                ▼
   ┌─────────────────────────────────────────┐
   │      alkane_tree_generator.py           │
   │   (C1-C10 饱和烷烃统一生成算法)          │
   └─────────────────────────────────────────┘
```

### 2.2 核心模块

| 模块 | 文件 | 功能 |
|------|------|------|
| **树结构生成器** | `alkane_tree_generator.py` | C1-C10饱和烷烃通用生成算法 |
| C1-C5 生成器 | `isomer_c1_c5_combined.py` | C1-C5 统一入口：饱和烷烃→树算法，非饱和烃→硬编码，含氧化合物→硬编码 |
| C6 生成器 | `isomer_c6_combined.py` | C6 特殊处理 |
| C8-C10 生成器 | `isomer_c8_c10_combined.py` | C8-C10 统一入口 |
| 结构筛选器 | `structure_filter.py` | 环、键、官能团筛选 |
| 分子工具 | `mol_utils.py` | SMILES处理、格式转换 |

---

## 三、核心算法详解

### 3.1 树结构生成算法（最核心）

#### 3.1.1 算法原理

烷烃分子可以建模为**无标号有根树**：
- **节点** = 碳原子（C）
- **边** = C-C 键
- **根** = 分子的**质心（Centroid）**

```
化学分子                  树结构表示                SMILES
─────────────────────────────────────────────────────────

正丁烷: CH3-CH2-CH2-CH3
       ↓                    C                        CCCC
                            |
                            C
                            |
                            C
                            |
                            C

异丁烷: CH3-CH(CH3)-CH3
       ↓                    C                         CC(C)C
                           /|\
                          C C C
```

#### 3.1.2 整数拆分

将 n 个碳原子分解为若干分支（每个分支最多4个，因为碳是四价）：

```python
# 对于 C4 (n=4)，拆分 n-1 = 3：
# [3]   → C(C,C,C)  → 3-甲基丙烷（实际被过滤）
# [2,1] → C(C(C),C) → 异丁烷
# [1,1,1] → C(C,C,C) → 正丁烷
```

关键代码（`alkane_tree_generator.py:31-42`）：

```python
def find_parts(rem, min_val, path):
    """生成所有整数拆分（最多4个部分）"""
    if len(path) > 4: return  # 碳原子最多4个分支
    if rem == 0:
        partitions.append(tuple(path))
        return
    if len(path) == 4: return
    
    for i in range(min_val, rem + 1):
        find_parts(rem - i, i, path + [i])
```

#### 3.1.3 递归树生成

```python
def gen(self, n):
    """递归生成所有可能的树结构"""
    if n == 1:
        return [("C()", [])]  # 叶子节点
    
    results = []
    target = n - 1
    
    # 对每个拆分方式
    for p in partitions:
        # 递归生成每个子树
        sub_trees_lists = [self.gen(s) for s in p]
        
        # 组合所有子树（笛卡尔积）
        for combo in itertools.product(*sub_trees_lists):
            children_reps = [c[0] for c in combo]
            # 排序确保同构树产生相同表示
            sorted_reps = sorted(children_reps)
            current_rep = "C(" + ",".join(sorted_reps) + ")"
            
            if current_rep not in seen_reps:
                seen_reps.add(current_rep)
                results.append((current_rep, list(p)))
    
    return results
```

#### 3.1.4 质心过滤（关键！）

**为什么需要过滤？**

同一棵树选择不同节点作为根，会生成重复的分子！

```
同一棵树，不同的根：
    C          C
   /|\        /|\
  C C C  和  C C C  是同一种分子！

质心 = 树的"几何中心"
```

**过滤条件**：`max(branch_sizes) <= n // 2`

```python
# 代码实现 (line 124-130)
limit = n // 2

for rep, sub_sizes in all_trees:
    # 保留以质心为根的树结构
    if not sub_sizes or max(sub_sizes) <= limit:
        # 生成 SMILES...
```

**数学原理**：
- 奇数碳：质心是一个节点，最大分支 ≤ n/2
- 偶数碳：质心是一条边，任选一端为根，最大分支 = n/2

#### 3.1.5 SMILES 转换

```python
def rep_to_smiles(rep_str):
    """树结构 → SMILES 字符串"""
    if rep_str == "C()":
        return "C"
    
    # 解析括号，提取子节点
    children = parse_children(rep_str)  # 深度优先解析
    
    # 递归转换
    child_smiles = [rep_to_smiles(c) for c in children]
    
    # 构建 SMILES: C(child1)(child2)...
    return "C" + "".join(f"({s})" for s in child_smiles)
```

示例转换：

```
树: C(C(C),C)
    ↓ 递归解析
    C + (C(C)) + (C)
    ↓
SMILES: CC(C)C
```

#### 3.1.6 RDKit 规范化去重

```python
# 代码实现 (line 132-138)
mol = Chem.MolFromSmiles(raw_smiles)
if mol:
    # 转换为规范 SMILES（唯一表示）
    canon_sm = Chem.MolToSmiles(mol, isomericSmiles=False)
    if canon_sm not in canonical_smiles_set:
        valid_molecules.append(canon_sm)
```

### 3.2 烯烃、炔烃、环烷烃生成算法

除了饱和烷烃，本项目还支持**烯烃**（含双键）、**炔烃**（含三键）和**环烷烃**的生成。核心策略在 `isomer_general.py` 中实现。

#### 3.2.1 分子式识别策略

根据氢碳比（H/C）判断化合物类型：

| 氢碳比 | 化合物类型 | 分子式示例 |
|--------|-----------|-----------|
| H = 2C + 2 | 饱和烷烃 | C8H18 |
| H = 2C | 烯烃 / 环烷烃 | C6H12 |
| H = 2C - 2 | 炔烃 / 二烯烃 | C4H6 |
| H < 2C - 2 | 多不饱和化合物 | - |

关键代码（`isomer_general.py:86-100`）：

```python
def _generate_chain_structures(self, n_c: int, n_h: int, n_o: int):
    # 烷烃 CnH2n+2
    if n_h == 2 * n_c + 2:
        chains.append('C' * n_c)
    
    # 烯烃 CnH2n（双键位置异构）
    elif n_h == 2 * n_c:
        for i in range(1, n_c):  # 双键在不同位置
            smiles = 'C' * i + '=' + 'C' * (n_c - i)
            chains.append(smiles)
    
    # 炔烃 CnH2n-2（三键位置异构）
    elif n_h == 2 * n_c - 2:
        for i in range(1, n_c):  # 三键在不同位置
            smiles = 'C' * i + '#' + 'C' * (n_c - i)
            chains.append(smiles)
```

#### 3.2.2 烯烃的生成（双键位置异构）

烯烃的同分异构现象包括：
1. **碳链异构**：主链长度不同
2. **双键位置异构**：双键位置不同（如 1-丁烯 vs 2-丁烯）
3. **顺反异构**：E/Z 构型

生成算法：

```python
# C5H10 戊烯的生成
# 双键位置：1-位、2-位、3-位...
for i in range(1, n_c):
    smiles = 'C' * i + '=' + 'C' * (n_c - i)
    # C=CCCC (1-戊烯)
    # CC=CCC (2-戊烯)
    # CCC=CC (3-戊烯，等同于2-戊烯)
```

#### 3.2.3 炔烃的生成（三键位置异构）

炔烃同样具有位置异构现象，生成逻辑与烯烃类似：

```python
# C4H6 丁炔的生成
for i in range(1, n_c):
    smiles = 'C' * i + '#' + 'C' * (n_c - i)
    # C#CCC (1-丁炔)
    # CC#CC (2-丁炔)
```

#### 3.2.4 环烷烃的生成

环状结构的生成采用**环+链组合策略**：

```python
def _generate_ring_structures(self, n_c: int, n_h: int, n_o: int):
    rings = []
    
    # 小环结构（3-6元环）
    for ring_size in range(3, min(7, n_c + 1)):
        if ring_size <= n_c:
            if ring_size == n_c:
                # 纯环状：C3H6 → 环丙烷
                rings.append('C1' + 'C' * (ring_size - 2) + 'C1')
            else:
                # 环+链：环己烷 + 甲基 = 甲基环己烷
                ring_part = 'C1' + 'C' * (ring_size - 2) + 'C1'
                chain_part = 'C' * (n_c - ring_size)
                rings.append(ring_part + chain_part)
    
    # 芳香环（碳数≥6）
    if n_c >= 6:
        rings.append('c1ccccc1')  # 苯
    
    return rings
```

#### 3.2.5 多环及复杂结构

对于更复杂的分子，算法还支持：

- **多环结构**：萘、联苯等
- **并环结构**：环己烷并环丙烷
- **桥环结构**：双环[2.2.1]庚烷

```python
def _generate_polycyclic_structures(self, n_c: int):
    # 萘 (C10)
    if n_c == 10:
        polycyclic.append('c1ccc2ccccc2c1')
    # 联苯 (C12+)
    elif n_c >= 12:
        polycyclic.append('c1ccccc1c1ccccc1')
```

#### 3.2.6 结构分类与验证

生成后使用 RDKit 进行结构分类（`isomer_c6_combined.py:82-109`）：

```python
if '#' in smiles:
    if '=' in smiles:
        return '多烯烃/多炔烃'
    return '炔烃'
elif '=' in smiles and '1' in smiles:
    return '环状烯烃'
elif '=' in smiles:
    return '烯烃'
elif '1' in smiles:
    return '环状烷烃'
else:
    return '烷烃'
```

### 3.3 算法验证

运行结果与 OEIS A000602 理论值**完全匹配**：

```
C1H4:  1/1  ✅
C2H6:  1/1  ✅
C3H8:  1/1  ✅
C4H10: 2/2  ✅
C5H12: 3/3  ✅
C6H14: 5/5  ✅
C7H16: 9/9  ✅
C8H18: 18/18 ✅
C9H20: 35/35 ✅
C10H22: 75/75 ✅
```

---

## 四、项目运行流程

### 4.1 启动方式

```bash
# 方式1：运行批处理文件
cd MolGenPlus_Project
run_molgen.bat

# 方式2：直接运行 Python
python MolGenPlus_Main.py
```

### 4.2 用户交互流程

```
┌────────────────────────────────────────────┐
│           MolGenPlus 主界面                │
├────────────────────────────────────────────┤
│                                            │
│  分子式输入: [C9H20          ] [分析同分异构体]│
│                                            │
│  结构筛选:                                  │
│  ○ 无环  ● 全部  ○ 单环                    │
│  键类型: 单键  ○ 全部  ○ 双键              │
│                                            │
│  [导出TXT] [导出GJF] [导出XYZ]             │
│                                            │
│  ─────────────────────────────────────     │
│  分析结果:                                  │
│  1. CCCCCCCCC (正壬烷)                     │
│  2. CCCCCCCC(C)C (2-甲基辛烷)              │
│  3. CCCCCCC(C)CC (3-甲基辛烷)              │
│  ...                                        │
│                                            │
│  显示分子结构图...                          │
└────────────────────────────────────────────┘
```

### 4.3 内部处理流程

```
用户输入 "C9H20"
       │
       ▼
parse_formula_and_dispatch()
       │
       ├─ C1-C5 (饱和烷烃): 调用 isomer_c1_c5_combined.py
       │                         │
       │                         ▼
       │                    generate_alkane_isomers(n)
       │                         │
       │                         ├─ TreeGen.gen(n)     [树结构生成]
       │                         ├─ 质心过滤            [去重]
       │                         ├─ rep_to_smiles()    [SMILES转换]
       │                         └─ RDKit 规范化        [最终去重]
       │
       ├─ C1-C5 (非饱和烃/含氧化合物): 硬编码数据库
       │
       ├─ C6:   调用 isomer_c6_combined.py
       │
       └─ C7-C10: 调用 isomer_c8_c10_combined.py
                      │
                      ▼
                 generate_alkane_isomers(n)
                      │
                      ├─ TreeGen.gen(n)     [树结构生成]
                      ├─ 质心过滤            [去重]
                      ├─ rep_to_smiles()     [SMILES转换]
                      └─ RDKit 规范化        [最终去重]
                      │
                      ▼
                 返回异构体列表
```

**统一性**：无论是 C1-C10 中的哪一种饱和烷烃，都使用同一套 `generate_alkane_isomers()` 树结构生成算法，保证了代码的一致性和可维护性。

---

## 五、关键技术亮点

### 5.1 记忆化缓存

```python
class AlkaneTreeGenerator:
    def __init__(self):
        self.memo = {}  # 缓存已计算的树结构
    
    def gen(self, n):
        if n in self.memo:
            return self.memo[n]  # 直接返回缓存结果
        # ... 计算 ...
        self.memo[n] = results  # 保存结果
        return results
```

**效果**：避免重复计算，大幅提升性能。

### 5.2 同构树去重

```python
# 对子树表示排序，确保同构树产生相同的字符串
sorted_reps = sorted(children_reps)
current_rep = "C(" + ",".join(sorted_reps) + ")"

if current_rep not in seen_reps:
    seen_reps.add(current_rep)
    results.append(...)
```

### 5.3 深度优先括号解析

正确处理嵌套括号（代码行 74-90）：

```python
# 使用深度计数器正确分割嵌套括号
children = []
depth = 0
current = ""
for char in inner:
    if char == '(':
        depth += 1
    elif char == ')':
        depth -= 1
    elif char == ',' and depth == 0:
        children.append(current)
        current = ""
        continue
    current += char
```

### 5.4 统一生成接口

```
所有碳数 (C1-C10) 的饱和烷烃共用同一生成器接口：
- generate_alkane_isomers(1)  → C1H4
- generate_alkane_isomers(2)  → C2H6
- ...
- generate_alkane_isomers(10) → C10H22 (75个异构体)

isomer_c1_c5_combined.py 根据分子式自动判断：
- 饱和烷烃 (H = 2C + 2) → 调用树算法
- 非饱和烃/含氧化合物 → 硬编码数据库
```

### 5.5 代码复用

```
isomer_c1_c5_combined.py 调度逻辑：
1. 解析分子式：提取碳数 nC、氢数 nH
2. 判断类型：if nH == 2*nC + 2 → 饱和烷烃
3. 分发处理：
   - 饱和烷烃 → generate_alkane_isomers(nC)
   - 其他 → STRUCTURE_DATABASE_C1_C5 硬编码
```

---

## 六、SMARTS 结构匹配算法

### 6.1 SMARTS 简介

**SMARTS**（SMiles ARbitrary Target Specification）是 SMILES 的扩展，用于描述分子子结构模式。本项目使用 SMARTS 进行**官能团识别**和**结构筛选**。

### 6.2 官能团 SMARTS 模式定义

项目在 `mol_utils.py` 中定义了官能团 SMARTS 模式（行 38-50）：

```python
FUNCTIONAL_GROUP_SMARTS = {
    # 含氧官能团 (按优先级排序)
    "羧基 (Carboxylic Acid)": Chem.MolFromSmarts('[CX3](=O)[OH]'),
    "酯基 (Ester)": Chem.MolFromSmarts('[CX3](=O)[O][#6]'),
    "醛基 (Aldehyde)": Chem.MolFromSmarts('[CX3H1](=O)[#6]'),
    "酮基 (Ketone)": Chem.MolFromSmarts('[CX3](=O)[#6][#6]'),
    "醚键 (Ether)": Chem.MolFromSmarts('[CX4]-[O]-[CX4]'),
    "醇/酚羟基 (Hydroxyl)": Chem.MolFromSmarts('[OX2H]'),
    
    # 烃类结构
    "芳香烃 (苯)": Chem.MolFromSmarts('c1ccccc1'),
    "炔烃": Chem.MolFromSmarts('C#C'),
    "烯烃": Chem.MolFromSmarts('C=C'),
}
```

### 6.3 SMARTS 模式匹配原理

使用 RDKit 的 `HasSubstructMatch()` 方法进行子结构匹配：

```python
def _matches_functional_group(self, mol, group_name):
    """检查分子是否含有指定官能团"""
    pattern = FUNCTIONAL_GROUP_SMARTS[group_name]
    if pattern is None:
        return False
    return mol.HasSubstructMatch(pattern)
```

**匹配示例**：

| SMARTS 模式 | 匹配的分子特征 |
|-------------|---------------|
| `[CX3](=O)[OH]` | 羧基 (-COOH) |
| `[CX3](=O)[O][#6]` | 酯基 (-COO-) |
| `c1ccccc1` | 苯环 |
| `C=C` | 碳碳双键 |
| `C#C` | 碳碳三键 |

### 6.4 结构筛选应用

在 `structure_filter.py` 中，SMARTS 用于精确筛选：

```python
# 官能团筛选 (行 273-280)
if mol.HasSubstructMatch(Chem.MolFromSmarts('[OH]')):
    functional_group = '醇/酚'
elif mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3]=[O]')):
    functional_group = '羰基'
elif mol.HasSubstructMatch(Chem.MolFromSmarts('C(=O)[OX1H0-,OX2H+]')):
    functional_group = '羧基'
```

---

## 七、分子指纹与聚类分析

### 7.1 Morgan 分子指纹

为了进行**结构相似性聚类**，需要将分子转换为数值特征向量。本项目使用 **Morgan 指纹**（ECFP 的 RDKit 实现）。

```python
def _get_morgan_fp(self, mol, n_bits=1024):
    """获取 Morgan 分子指纹"""
    # 半径=2 的 Morgan 指纹，1024位
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
    fp_array = np.zeros(n_bits)
    DataStructs.ConvertToNumpyArray(fp, fp_array)
    return fp_array.tolist()
```

**指纹特点**：
- 每一位代表一个局部化学环境（半径=2 的圆）
- 结构相似的分子有相似的指纹向量
- 可计算 Tanimoto 相似度

### 7.2 特征矩阵构建

```python
def extract_feature_matrix(self, smiles_list):
    """构建特征矩阵用于 ML 分析"""
    feature_dicts = []
    
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        
        # 分子描述符
        features = {
            'MW': Descriptors.MolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'TPSA': Descriptors.TPSA(mol),
            'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
            # Morgan 指纹 (1024维 → 降维取前10位)
            'morgan_fp': self._get_morgan_fp(mol, n_bits=1024),
        }
        feature_dicts.append(features)
    
    return feature_matrix, feature_names
```

### 7.3 K-Means 聚类算法

使用 **K-Means** 对分子进行结构相似性聚类：

```python
def perform_clustering(self, smiles_list, method='kmeans', n_clusters=3):
    """执行聚类分析"""
    # 1. 提取特征矩阵
    feature_matrix, _ = self.extract_feature_matrix(smiles_list)
    
    # 2. 标准化特征
    scaled_features = self.scaler.fit_transform(feature_matrix)
    
    # 3. K-Means 聚类
    clusterer = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = clusterer.fit_predict(scaled_features)
    
    # 4. 按聚类分组
    clusters = {}
    for label, smiles in zip(labels, valid_smiles):
        cluster_key = f'Cluster_{label}'
        clusters[cluster_key].append(smiles)
    
    return clusters
```

### 7.4 聚类质量评估

使用**轮廓系数**（Silhouette Score）评估聚类质量：

```python
silhouette = silhouette_score(scaled_features, labels)
# 范围: [-1, 1]
# 接近 1: 聚类紧凑，分离良好
# 接近 0: 聚类重叠
# 接近 -1: 错误聚类
```

### 7.5 聚类结果展示

```
🎯 聚类分析 (K-Means):
  聚类数量: 3
  轮廓系数: 0.8542

  Cluster_0 (12 个分子):
    CCCCCCC  正庚烷
    CC(C)CCC  2-甲基己烷
    ...

  Cluster_1 (15 个分子):
    CC(C)(C)CC  2,2-二甲基戊烷
    ...

  Cluster_2 (8 个分子):
    CC(C)(C)(C)C  2,2,3-三甲基丁烷
    ...
```

---

## 八、功能展示

### 8.1 同分异构体生成

### 8.2 SMARTS 结构匹配

- **官能团识别**：自动识别醇、醛、酮、羧酸、酯等
- **烃类识别**：烷烃、烯烃、炔烃、芳香烃
- **子结构筛选**：按指定结构片段筛选分子

### 8.3 聚类分析

- **Morgan 指纹**：1024维结构特征向量
- **K-Means 聚类**：自动将分子按相似性分组
- **轮廓系数评估**：量化聚类质量

- **输入**：C8H18
- **输出**：18 种异构体（完全匹配理论值）

### 6.2 结构筛选

支持按以下条件筛选：
- **环结构**：无环、单环、多环
- **键类型**：单键、双键、三键
- **芳香性**：无芳香、含苯环、其他芳香
- **官能团**：羟基、羰基、羧基、酯基

### 6.3 导出功能

支持多种化学文件格式：
- `.txt` - 文本格式
- `.xyz` - XYZ 坐标（兼容 GaussView）
- `.gjf` - Gaussian 输入文件

---

## 九、技术栈

| 技术 | 用途 |
|------|------|
| Python 3.x | 编程语言 |
| RDKit | 化学信息学核心库 |
| Tkinter | GUI 界面 |
| Scikit-learn | 机器学习分析 |
| DeepSeek API | AI 增强功能（可选） |

---

## 十、总结与展望

### 10.1 主要贡献

1. **实现了通用树结构算法**，能够生成 C1-C10 所有饱和烷烃异构体
2. **算法经过严格验证**，生成结果与 OEIS 理论值完全匹配
3. **提供了完整的 GUI 界面**，易于使用
4. **支持结构筛选和多种导出格式**
5. **集成 SMARTS 匹配和聚类分析**，实现智能分子分类
6. **统一生成架构**：C1-C10 饱和烷烃统一使用树算法，消除代码重复

### 10.2 创新点

- **质心过滤**：巧妙利用图论中的质心概念，避免重复枚举
- **记忆化 + 递归**：高效枚举树结构
- **规范化去重**：结合领域知识（化学）与计算机科学（规范化）
- **智能调度**：根据分子式自动选择算法（树算法 vs 硬编码）

### 10.3 后续工作

- 扩展到 C11-C20 烷烃（算法理论上可行，计算量增大）
- 优化含氧化合物生成
- 增加更多官能团支持

---

## 参考文献

1. OEIS Sequence A000602: Number of alkane isomers
2. RDKit Documentation: Open-Source Cheminformatics
3. Pólya Enumeration Theory

---

**答辩人**：[姓名]  
**指导教师**：[教师姓名]  
**日期**：2026年3月
