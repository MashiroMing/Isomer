"""
烷烃同分异构体理论计数模块
提供已知烷烃的同分异构体数量统计功能
"""

def get_alkane_isomer_statistics(formula: str) -> dict:
    """
    获取烷烃同分异构体的理论统计信息

    Args:
        formula: 分子式，如 C2H6, C3H8, C4H10 等

    Returns:
        包含统计信息的字典:
        - isomer_count: 理论异构体数量
        - density: 每个碳原子的平均异构体数
        - error: 错误信息（如果有）
    """
    import re

    # 解析分子式
    match_c = re.search(r'C(\d+)', formula)
    match_h = re.search(r'H(\d+)', formula)

    if not match_c:
        return {"error": "无法解析碳原子数"}

    try:
        n_c = int(match_c.group(1))
        n_h = int(match_h.group(1)) if match_h else 0
    except ValueError:
        return {"error": "分子式格式错误"}

    # 验证是否为烷烃 (饱和烃: CnH(2n+2))
    expected_h = 2 * n_c + 2
    if n_h != expected_h:
        return {"error": f"不是标准烷烃 (C{n_c}H{expected_h})"}

    # 理论异构体数量表 (C1-C20)
    isomer_counts = {
        1: 1,    # 甲烷
        2: 1,    # 乙烷
        3: 1,    # 丙烷
        4: 2,    # 丁烷
        5: 3,    # 戊烷
        6: 5,    # 己烷
        7: 9,    # 庚烷
        8: 18,   # 辛烷
        9: 35,   # 壬烷
        10: 75,  # 癸烷
        11: 159,
        12: 355,
        13: 802,
        14: 1858,
        15: 4347,
        16: 10359,
        17: 24894,
        18: 60523,
        19: 148284,
        20: 366319
    }

    isomer_count = isomer_counts.get(n_c, 0)
    density = isomer_count / n_c if n_c > 0 else 0

    return {
        "isomer_count": isomer_count,
        "density": density
    }
