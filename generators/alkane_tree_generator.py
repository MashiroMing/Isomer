# alkane_tree_generator.py
# 通用烷烃同分异构体生成器 (C1-C10)
# 基于树结构生成 + RDKit 规范化

import itertools
from rdkit import Chem
from typing import List, Dict, Any

class AlkaneTreeGenerator:
    """基于树结构生成烷烃异构体的生成器"""
    
    def __init__(self):
        self.memo = {} 

    def gen(self, n):
        """生成所有节点数为n的树，每个节点度数<=4"""
        if n in self.memo:
            return self.memo[n]
        
        if n == 1:
            res = [("C()", [])]
            self.memo[n] = res
            return res
        
        results = []
        target = n - 1
        seen_reps = set()
        
        # 生成整数拆分 (Partitions of target into at most 4 parts)
        partitions = []
        def find_parts(rem, min_val, path):
            if len(path) > 4: return
            if rem == 0:
                if len(path) > 0:
                    partitions.append(tuple(path))
                return
            if len(path) == 4: return
            
            for i in range(min_val, rem + 1):
                find_parts(rem - i, i, path + [i])
        
        find_parts(target, 1, [])
        
        for p in partitions:
            sub_trees_lists = [self.gen(s) for s in p]
            
            for combo in itertools.product(*sub_trees_lists):
                children_reps = [c[0] for c in combo]
                # 排序以确保同构的树产生相同的表示字符串
                sorted_reps = sorted(children_reps)
                current_rep = "C(" + ",".join(sorted_reps) + ")"
                
                if current_rep not in seen_reps:
                    seen_reps.add(current_rep)
                    results.append((current_rep, list(p)))
        
        self.memo[n] = results
        return results


def rep_to_smiles(rep_str):
    """将内部树表示转换为 SMILES 字符串"""
    if rep_str == "C()":
        return "C"
    
    if not rep_str.startswith("C(") or not rep_str.endswith(")"):
        return "C"
    
    inner = rep_str[2:-1]
    if not inner:
        return "C"
    
    # 正确分割嵌套括号
    children = []
    depth = 0
    current = ""
    for char in inner:
        if char == '(':
            depth += 1
            current += char
        elif char == ')':
            depth -= 1
            current += char
        elif char == ',' and depth == 0:
            children.append(current)
            current = ""
        else:
            current += char
    if current:
        children.append(current)
    
    child_smiles_list = [rep_to_smiles(c) for c in children]
    
    if not child_smiles_list:
        return "C"
    
    # 将所有子树放在括号中，确保连接在根原子上
    # 按长度排序仅为了输出美观，不影响化学结构
    child_smiles_list.sort(key=len, reverse=True)
    
    result = "C"
    for child in child_smiles_list:
        result += f"({child})"
        
    return result


def generate_alkane_isomers(n: int) -> List[str]:
    """
    生成碳数为 n 的烷烃异构体
    
    Args:
        n: 碳原子数 (1-10)
    
    Returns:
        SMILES 字符串列表
    """
    tg = AlkaneTreeGenerator()
    all_trees = tg.gen(n)
    
    canonical_smiles_set = set()
    valid_molecules = []
    
    # 关键修正：过滤条件改为 max(sub_sizes) <= n // 2
    # 这同时适用于奇数 (质心为点) 和偶数 (质心为边，取其中一点为根时最大子树为 n/2)
    limit = n // 2
    
    for rep, sub_sizes in all_trees:
        # 处理空列表情况（只有一个节点时）
        if not sub_sizes or max(sub_sizes) <= limit:
            raw_smiles = rep_to_smiles(rep)
            try:
                mol = Chem.MolFromSmiles(raw_smiles)
                if mol:
                    canon_sm = Chem.MolToSmiles(mol, isomericSmiles=False)
                    if canon_sm not in canonical_smiles_set:
                        canonical_smiles_set.add(canon_sm)
                        valid_molecules.append(canon_sm)
            except Exception:
                continue
    
    return sorted(valid_molecules, key=lambda x: (len(x), x))


# 理论参考值 (OEIS A000602)
ALKANE_THEORETICAL_COUNTS = {
    1: 1, 2: 1, 3: 1, 4: 2, 5: 3, 
    6: 5, 7: 9, 8: 18, 9: 35, 10: 75
}


def get_theoretical_count(n: int) -> int:
    """获取碳数 n 的烷烃理论异构体数量"""
    return ALKANE_THEORETICAL_COUNTS.get(n, 0)


if __name__ == "__main__":
    # 测试代码
    print("=" * 60)
    print("烷烃同分异构体生成器 (C1 - C10)")
    print("=" * 60)
    
    for n in range(1, 11):
        isomers = generate_alkane_isomers(n)
        theory = get_theoretical_count(n)
        status = "✅" if len(isomers) == theory else "❌"
        print(f"C{n}H{2*n+2}: {len(isomers)}/{theory} {status}")
    
    print("\n详细列表 (C10):")
    c10 = generate_alkane_isomers(10)
    for i, sm in enumerate(c10, 1):
        print(f"{i:2d}: {sm}")
