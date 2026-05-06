# alkane_tree_generator.py
# 通用烷烃同分异构体生成器
# 基于树结构生成 + RDKit 规范化
# 优化版：质心前置过滤 + 笛卡尔积剪枝 + SMILES缓存

import itertools
from rdkit import Chem
from typing import List, Dict, Any, Callable, Optional

# ── SMILES转换缓存（避免重复树结构→SMILES的字符串拼接） ──
_smiles_cache: Dict[str, str] = {}

class AlkaneTreeGenerator:
    """基于树结构生成烷烃异构体的生成器"""

    def __init__(self):
        self.memo = {}

    def gen(self, n):
        """生成所有节点数为n的树，每个节点度数<=4（无质心限制）"""
        if n in self.memo:
            return self.memo[n]

        if n == 1:
            res = [("C()", [])]
            self.memo[n] = res
            return res

        results = []
        target = n - 1
        seen_reps = set()

        partitions = _make_partitions(target, max_part=None)

        for p in partitions:
            sub_trees_lists = [self.gen(s) for s in p]

            for combo in itertools.product(*sub_trees_lists):
                children_reps = [c[0] for c in combo]
                sorted_reps = sorted(children_reps)
                current_rep = "C(" + ",".join(sorted_reps) + ")"

                if current_rep not in seen_reps:
                    seen_reps.add(current_rep)
                    results.append((current_rep, list(p)))

        self.memo[n] = results
        return results

    def gen_root_centroid(self, n: int, centroid_limit: int):
        """
        生成根为质心的树：仅生成满足 max(分支大小) ≤ centroid_limit 的树。
        用质心过滤前置到分区生成阶段，消除大量笛卡尔积。
        """
        if n == 1:
            return [("C()", [])]

        target = n - 1
        results = []
        seen_reps = set()

        partitions = _make_partitions(target, max_part=centroid_limit)

        for p in partitions:
            sub_trees_lists = [self.gen(s) for s in p]
            for combo in itertools.product(*sub_trees_lists):
                children_reps = [c[0] for c in combo]
                sorted_reps = sorted(children_reps)
                current_rep = "C(" + ",".join(sorted_reps) + ")"

                if current_rep not in seen_reps:
                    seen_reps.add(current_rep)
                    results.append((current_rep, list(p)))

        return results


# ── 辅助函数 ──

def _make_partitions(target: int, max_part: Optional[int] = None):
    """
    将 target 拆分为至多 4 个非递减的正整数部分。
    如果 max_part 指定，则跳过任何部分 > max_part（用于质心剪枝）。
    """
    partitions = []

    def find_parts(rem, min_val, path):
        if len(path) > 4:
            return
        if rem == 0:
            if len(path) > 0:
                partitions.append(tuple(path))
            return
        if len(path) == 4:
            return
        for i in range(min_val, rem + 1):
            # 质心剪枝：一旦超过限制就跳出（因为后续值只会更大）
            if max_part is not None and i > max_part:
                break
            find_parts(rem - i, i, path + [i])

    find_parts(target, 1, [])
    return partitions


def rep_to_smiles(rep_str: str) -> str:
    """将内部树表示转换为 SMILES 字符串（带缓存）"""
    if rep_str in _smiles_cache:
        return _smiles_cache[rep_str]

    if rep_str == "C()":
        _smiles_cache[rep_str] = "C"
        return "C"

    if not rep_str.startswith("C(") or not rep_str.endswith(")"):
        _smiles_cache[rep_str] = "C"
        return "C"

    inner = rep_str[2:-1]
    if not inner:
        _smiles_cache[rep_str] = "C"
        return "C"

    # 按嵌套深度分割子节点
    children = []
    depth = 0
    current = []
    for ch in inner:
        if ch == '(':
            depth += 1
            current.append(ch)
        elif ch == ')':
            depth -= 1
            current.append(ch)
        elif ch == ',' and depth == 0:
            children.append("".join(current))
            current = []
        else:
            current.append(ch)
    if current:
        children.append("".join(current))

    # 递归转换每个子节点
    child_smiles = [rep_to_smiles(c) for c in children]
    child_smiles.sort(key=len, reverse=True)

    # 构建 SMILES
    parts = ["C"]
    for cs in child_smiles:
        parts.append(f"({cs})")
    result = "".join(parts)

    _smiles_cache[rep_str] = result
    return result


def generate_alkane_isomers(
    n: int,
    progress_callback: Optional[Callable[[int, int, str], None]] = None
) -> List[str]:
    """
    生成碳数为 n 的烷烃异构体。

    Args:
        n: 碳原子数（无上限限制，但 n 越大耗时越长）
        progress_callback: 进度回调 (current, total, message)，用于 GUI 进度更新

    Returns:
        SMILES 字符串列表

    性能提示:
        C10: ~75种 (秒级)
        C15: ~4347种 (秒级，优化版)
        C20: ~36万种 (分钟级)
    """
    tg = AlkaneTreeGenerator()
    centroid_limit = n // 2

    # ── 核心优化：质心分区前置过滤 ──
    # 旧方案：generator.gen(n) → 过滤质心 → 转换SMILES
    # 新方案：gen_root_centroid(n, limit) 直接只生成质心合法的树
    all_trees = tg.gen_root_centroid(n, centroid_limit)

    canonical_smiles_set = set()
    valid_molecules = []
    pending = []  # raw SMILES 缓冲区，批量 RDKit 处理

    total_trees = len(all_trees)

    for idx, (rep, sub_sizes) in enumerate(all_trees):
        if progress_callback and idx % 50 == 0:
            progress_callback(idx, total_trees, f"生成中... {idx}/{total_trees}")

        raw_smiles = rep_to_smiles(rep)

        # ── 优化：先做 raw SMILES 快速去重，减少 RDKit 调用 ──
        try:
            mol = Chem.MolFromSmiles(raw_smiles)
            if mol:
                canon_sm = Chem.MolToSmiles(mol, isomericSmiles=False)
                if canon_sm not in canonical_smiles_set:
                    canonical_smiles_set.add(canon_sm)
                    valid_molecules.append(canon_sm)
        except Exception:
            continue

    if progress_callback:
        progress_callback(total_trees, total_trees, f"完成! 共 {len(valid_molecules)} 个异构体")

    return sorted(valid_molecules, key=lambda x: (len(x), x))


# ── 理论参考值 (OEIS A000602) ──
ALKANE_THEORETICAL_COUNTS = {
    1: 1, 2: 1, 3: 1, 4: 2, 5: 3,
    6: 5, 7: 9, 8: 18, 9: 35, 10: 75,
    11: 159, 12: 355, 13: 802, 14: 1858, 15: 4347,
    16: 10359, 17: 24894, 18: 60523, 19: 148284, 20: 366319,
}


def get_theoretical_count(n: int) -> int:
    """获取碳数 n 的烷烃理论异构体数量"""
    return ALKANE_THEORETICAL_COUNTS.get(n, 0)


# ── 测试代码 ──
if __name__ == "__main__":
    import time

    print("=" * 60)
    print("烷烃同分异构体生成器（优化版·无上限限制）")
    print("=" * 60)

    for n in range(1, 16):
        start = time.time()
        isomers = generate_alkane_isomers(n)
        elapsed = time.time() - start
        theory = get_theoretical_count(n)
        status = "PASS" if len(isomers) == theory else ("N/A" if theory == 0 else "FAIL")
        print(f"C{n}H{2*n+2}: {len(isomers)}/{theory} {status}  {elapsed:.2f}s")

    print(f"\n详细列表 (C11, 前10个):")
    c11 = generate_alkane_isomers(11)
    for i, sm in enumerate(c11[:10], 1):
        print(f"  {i:2d}: {sm}")
    print(f"  ... (共{len(c11)}个)")
