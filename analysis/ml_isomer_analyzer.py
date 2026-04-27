# -*- coding: utf-8 -*-
"""
机器学习分子分析模块
提供分子特征提取、聚类分析和预测功能
兼容MolGenPlus项目架构
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, rdMolDescriptors, DataStructs
from rdkit.Chem.Draw import MolToImage
import numpy as np
import warnings
import os

# 静默各种警告
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# 抑制OpenMP和sklearn的警告
os.environ['OMP_NUM_THREADS'] = '1'

# 静默RDKit警告
try:
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
except:
    pass

# 机器学习库导入
try:
    from sklearn.cluster import KMeans, DBSCAN
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestClassifier, IsolationForest
    from sklearn.svm import SVC
    from sklearn.model_selection import train_test_split, cross_val_score
    from sklearn.metrics import classification_report, silhouette_score
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
    warnings.warn("scikit-learn未安装，部分机器学习功能将不可用。请运行: pip install scikit-learn")

# 可视化库导入
try:
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    warnings.warn("matplotlib未安装，可视化功能将不可用。请运行: pip install matplotlib")

# PIL导入用于图像处理
try:
    from PIL import Image
    HAS_PIL = True
except ImportError:
    HAS_PIL = False
    warnings.warn("PIL/Pillow未安装，图像处理功能将不可用。请运行: pip install Pillow")


class MLIsomerAnalyzer:
    """机器学习分子分析器"""

    def __init__(self):
        """初始化分析器"""
        self.features_cache = {}
        self.scaler = StandardScaler() if HAS_SKLEARN else None
        self.model_cache = {}

    def extract_molecular_features(self, smiles: str) -> dict:
        """
        提取分子特征

        Args:
            smiles: SMILES字符串

        Returns:
            特征字典
        """
        if smiles in self.features_cache:
            return self.features_cache[smiles].copy()

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        features = {
            # 基础描述符
            'mol_weight': Descriptors.MolWt(mol),
            'num_atoms': mol.GetNumAtoms(),
            'num_heavy_atoms': mol.GetNumHeavyAtoms(),
            'num_bonds': mol.GetNumBonds(),
            'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'num_rings': Descriptors.RingCount(mol),

            # 电子性质
            'logp': Descriptors.MolLogP(mol),
            'tpsa': Descriptors.TPSA(mol),
            'num_h_donors': Descriptors.NumHDonors(mol),
            'num_h_acceptors': Descriptors.NumHAcceptors(mol),
            'formal_charge': Chem.GetFormalCharge(mol),

            # 拓扑描述符
            'bertz_ct': Descriptors.BertzCT(mol),
            'weiner_index': self._calculate_weiner_index(mol),

            # 碳骨架特征
            'num_carbon': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C'),
            'num_hydrogen': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'H'),
            'num_hetero': mol.GetNumHeavyAtoms() - sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C'),

            # 支链特征
            'num_branches': self._count_branches(mol),
            'max_branch_length': self._get_max_branch_length(mol),

            # 芳香性
            'num_aromatic_rings': Descriptors.NumAromaticRings(mol),
            'num_aromatic_bonds': sum(1 for bond in mol.GetBonds() if bond.GetIsAromatic()),

            # 指纹特征 (Morgan指纹的前10位作为示例)
            'morgan_fp': self._get_morgan_fp(mol, n_bits=10),

            # SMILES特征
            'smiles_length': len(smiles),
            'num_rings_in_smiles': smiles.count('1'),
        }

        self.features_cache[smiles] = features.copy()
        return features.copy()

    def _calculate_weiner_index(self, mol) -> float:
        """计算Wiener指数"""
        try:
            return rdMolDescriptors.CalcWienerIndex(mol)
        except:
            return 0.0

    def _count_branches(self, mol) -> int:
        """计算支链数量"""
        branches = 0
        for atom in mol.GetAtoms():
            if atom.GetDegree() > 2:
                branches += atom.GetDegree() - 2
        return branches

    def _get_max_branch_length(self, mol) -> int:
        """获取最大支链长度"""
        max_length = 0
        for atom in mol.GetAtoms():
            if atom.GetDegree() > 2:
                for neighbor in atom.GetNeighbors():
                    length = self._get_path_length(mol, atom.GetIdx(), neighbor.GetIdx())
                    if length > max_length:
                        max_length = length
        return max_length

    def _get_path_length(self, mol, start_idx, end_idx, visited=None) -> int:
        """获取两个原子之间的路径长度"""
        if visited is None:
            visited = set()

        if start_idx == end_idx:
            return 0

        if start_idx in visited:
            return 0

        visited.add(start_idx)
        max_path = 0

        for neighbor in mol.GetAtomWithIdx(start_idx).GetNeighbors():
            path = self._get_path_length(mol, neighbor.GetIdx(), end_idx, visited.copy())
            if path > max_path:
                max_path = path

        return max_path + 1 if max_path > 0 else 0

    def _get_morgan_fp(self, mol, n_bits=1024) -> list:
        """获取Morgan指纹"""
        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)
            fp_array = np.zeros(n_bits)
            DataStructs.ConvertToNumpyArray(fp, fp_array)
            return fp_array.tolist()
        except:
            return [0] * n_bits

    def extract_feature_matrix(self, smiles_list: list) -> tuple:
        """
        提取SMILES列表的特征矩阵

        Args:
            smiles_list: SMILES字符串列表

        Returns:
            (特征矩阵, 特征名称列表)
        """
        feature_dicts = []
        valid_smiles = []

        for smiles in smiles_list:
            features = self.extract_molecular_features(smiles)
            if features is not None:
                feature_dicts.append(features)
                valid_smiles.append(smiles)

        if not feature_dicts:
            return None, None

        # 扁平化特征字典
        all_features = []
        for fd in feature_dicts:
            features_vector = []
            for key, value in sorted(fd.items()):
                if key != 'morgan_fp':  # 指纹单独处理
                    features_vector.append(float(value))
            # 添加指纹特征
            features_vector.extend(fd['morgan_fp'][:10])  # 只取前10位作为示例
            all_features.append(features_vector)

        feature_matrix = np.array(all_features)
        feature_names = [key for key in sorted(feature_dicts[0].keys()) if key != 'morgan_fp']
        feature_names.extend([f'morgan_fp_{i}' for i in range(10)])

        return feature_matrix, feature_names

    def perform_pca_analysis(self, smiles_list: list, n_components=2) -> dict:
        """
        执行PCA降维分析

        Args:
            smiles_list: SMILES字符串列表
            n_components: 主成分数量

        Returns:
            分析结果字典
        """
        if not HAS_SKLEARN:
            return {'error': 'scikit-learn未安装'}

        feature_matrix, feature_names = self.extract_feature_matrix(smiles_list)
        if feature_matrix is None:
            return {'error': '无法提取特征矩阵'}

        # 标准化
        scaled_features = self.scaler.fit_transform(feature_matrix)

        # PCA
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(scaled_features)

        return {
            'pca_result': pca_result,
            'explained_variance_ratio': pca.explained_variance_ratio_.tolist(),
            'cumulative_variance': np.cumsum(pca.explained_variance_ratio_).tolist(),
            'feature_importance': self._get_pca_feature_importance(pca, feature_names),
            'valid_smiles': [smiles for smiles, features in zip(smiles_list,
                       [self.extract_molecular_features(s) for s in smiles_list])
                       if features is not None]
        }

    def _get_pca_feature_importance(self, pca, feature_names) -> dict:
        """获取PCA特征重要性"""
        importance = {}
        for i, component in enumerate(pca.components_):
            component_importance = dict(zip(feature_names, np.abs(component)))
            sorted_importance = dict(sorted(component_importance.items(),
                                           key=lambda x: x[1], reverse=True)[:10])
            importance[f'PC{i+1}'] = sorted_importance
        return importance

    def perform_clustering(self, smiles_list: list, method='kmeans',
                          n_clusters=3, eps=0.5, min_samples=5) -> dict:
        """
        执行聚类分析

        Args:
            smiles_list: SMILES字符串列表
            method: 聚类方法 ('kmeans' 或 'dbscan')
            n_clusters: K-means聚类数量
            eps: DBSCAN的eps参数
            min_samples: DBSCAN的min_samples参数

        Returns:
            聚类结果字典
        """
        if not HAS_SKLEARN:
            return {'error': 'scikit-learn未安装'}

        feature_matrix, feature_names = self.extract_feature_matrix(smiles_list)
        if feature_matrix is None:
            return {'error': '无法提取特征矩阵'}

        # 标准化
        scaled_features = self.scaler.fit_transform(feature_matrix)

        valid_indices = [i for i, smiles in enumerate(smiles_list)
                         if self.extract_molecular_features(smiles) is not None]
        valid_smiles = [smiles_list[i] for i in valid_indices]

        # 聚类
        if method == 'kmeans':
            clusterer = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            labels = clusterer.fit_predict(scaled_features)
            method_name = 'K-Means'
        elif method == 'dbscan':
            clusterer = DBSCAN(eps=eps, min_samples=min_samples)
            labels = clusterer.fit_predict(scaled_features)
            method_name = 'DBSCAN'
        else:
            return {'error': f'不支持的聚类方法: {method}'}

        # 计算轮廓系数（仅用于K-means且有多个簇时）
        silhouette = None
        if method == 'kmeans' and len(set(labels)) > 1:
            try:
                silhouette = silhouette_score(scaled_features, labels)
            except:
                pass

        # 按聚类分组
        clusters = {}
        for label, smiles in zip(labels, valid_smiles):
            cluster_key = f'Cluster_{label}'
            if cluster_key not in clusters:
                clusters[cluster_key] = []
            clusters[cluster_key].append(smiles)

        return {
            'method': method_name,
            'labels': labels.tolist(),
            'clusters': clusters,
            'n_clusters': len(set(labels)),
            'silhouette_score': silhouette,
            'valid_smiles': valid_smiles
        }

    def detect_outliers(self, smiles_list: list, contamination=0.1) -> dict:
        """
        检测异常值分子

        Args:
            smiles_list: SMILES字符串列表
            contamination: 异常值比例

        Returns:
            异常值检测结果
        """
        if not HAS_SKLEARN:
            return {'error': 'scikit-learn未安装'}

        feature_matrix, _ = self.extract_feature_matrix(smiles_list)
        if feature_matrix is None:
            return {'error': '无法提取特征矩阵'}

        # 标准化
        scaled_features = self.scaler.fit_transform(feature_matrix)

        # 异常检测
        iso_forest = IsolationForest(contamination=contamination, random_state=42)
        outliers = iso_forest.fit_predict(scaled_features)

        valid_smiles = [smiles for smiles in smiles_list
                        if self.extract_molecular_features(smiles) is not None]

        outlier_molecules = [smiles for smiles, pred in zip(valid_smiles, outliers) if pred == -1]
        normal_molecules = [smiles for smiles, pred in zip(valid_smiles, outliers) if pred == 1]

        return {
            'outliers': outlier_molecules,
            'normal': normal_molecules,
            'n_outliers': len(outlier_molecules),
            'outlier_ratio': len(outlier_molecules) / len(valid_smiles) if valid_smiles else 0
        }

    def predict_properties(self, smiles_list: list, property_name: str = 'logp') -> dict:
        """
        预测分子属性（基于简单回归）

        Args:
            smiles_list: SMILES字符串列表
            property_name: 属性名称 ('logp', 'tpsa', 'mol_weight')

        Returns:
            预测结果
        """
        feature_matrix, _ = self.extract_feature_matrix(smiles_list)
        if feature_matrix is None:
            return {'error': '无法提取特征矩阵'}

        valid_smiles = [smiles for smiles in smiles_list
                        if self.extract_molecular_features(smiles) is not None]

        # 使用实际的描述符值（实际应用中可使用训练好的模型）
        predictions = []
        for smiles in valid_smiles:
            features = self.extract_molecular_features(smiles)
            if property_name == 'logp':
                predictions.append(features['logp'])
            elif property_name == 'tpsa':
                predictions.append(features['tpsa'])
            elif property_name == 'mol_weight':
                predictions.append(features['mol_weight'])
            else:
                predictions.append(0.0)

        return {
            'property': property_name,
            'smiles': valid_smiles,
            'predictions': predictions,
            'mean': np.mean(predictions),
            'std': np.std(predictions)
        }

    def analyze_isomer_diversity(self, smiles_list: list) -> dict:
        """
        分析同分异构体的多样性

        Args:
            smiles_list: SMILES字符串列表

        Returns:
            多样性分析结果
        """
        feature_matrix, feature_names = self.extract_feature_matrix(smiles_list)
        if feature_matrix is None:
            return {'error': '无法提取特征矩阵'}

        # 标准化
        scaled_features = self.scaler.fit_transform(feature_matrix)

        # 计算多样性指标
        n_molecules = len(smiles_list)
        n_features = feature_matrix.shape[1]

        # 特征统计
        feature_means = np.mean(scaled_features, axis=0).tolist()
        feature_stds = np.std(scaled_features, axis=0).tolist()

        # 距离矩阵
        from scipy.spatial.distance import pdist, squareform
        distances = pdist(scaled_features, 'euclidean')
        mean_distance = np.mean(distances)
        max_distance = np.max(distances)
        min_distance = np.min(distances)

        return {
            'n_molecules': n_molecules,
            'n_features': n_features,
            'feature_means': feature_means[:10],  # 只显示前10个
            'feature_stds': feature_stds[:10],
            'mean_pairwise_distance': float(mean_distance),
            'max_pairwise_distance': float(max_distance),
            'min_pairwise_distance': float(min_distance),
            'diversity_score': float(mean_distance)  # 简单的多样性指标
        }

    def visualize_clusters(self, smiles_list: list, method='kmeans', n_clusters=3) -> str:
        """
        可视化聚类结果（返回图像路径）

        Args:
            smiles_list: SMILES字符串列表
            method: 聚类方法
            n_clusters: 聚类数量

        Returns:
            图像文件路径
        """
        if not HAS_MATPLOTLIB or not HAS_SKLEARN:
            return None

        # 执行PCA和聚类
        pca_result = self.perform_pca_analysis(smiles_list, n_components=2)
        clustering_result = self.perform_clustering(smiles_list, method=method, n_clusters=n_clusters)

        if 'error' in pca_result or 'error' in clustering_result:
            return None

        # 创建图像
        fig, ax = plt.subplots(figsize=(10, 8))
        colors = plt.cm.tab10(np.linspace(0, 1, clustering_result['n_clusters']))

        labels = clustering_result['labels']
        valid_smiles = clustering_result['valid_smiles']

        for i, (x, y) in enumerate(pca_result['pca_result']):
            label = labels[i]
            ax.scatter(x, y, c=[colors[label]], s=100, alpha=0.6, edgecolors='black', linewidths=0.5)
            ax.annotate(valid_smiles[i], (x, y), fontsize=8, alpha=0.7)

        ax.set_xlabel(f'PC1 (Explained Variance: {pca_result["explained_variance_ratio"][0]:.2%})')
        ax.set_ylabel(f'PC2 (Explained Variance: {pca_result["explained_variance_ratio"][1]:.2%})')
        ax.set_title(f'Molecular Clustering - {clustering_result["method"]}\nSilhouette Score: {clustering_result.get("silhouette_score", "N/A")}')
        ax.legend([f'Cluster {i}' for i in range(clustering_result['n_clusters'])])
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        # 保存图像
        import os
        output_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ml_results')
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f'cluster_visualization_{method}.png')
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        return output_path

    def generate_analysis_report(self, smiles_list: list) -> dict:
        """
        生成综合分析报告

        Args:
            smiles_list: SMILES字符串列表

        Returns:
            综合分析报告
        """
        report = {
            'summary': {},
            'features': {},
            'clustering': {},
            'diversity': {},
            'outliers': {}
        }

        # 基础统计
        valid_smiles = [smiles for smiles in smiles_list
                        if self.extract_molecular_features(smiles) is not None]
        report['summary'] = {
            'total_input': len(smiles_list),
            'valid_molecules': len(valid_smiles),
            'invalid_molecules': len(smiles_list) - len(valid_smiles)
        }

        if not valid_smiles:
            return report

        # 特征分析
        feature_matrix, feature_names = self.extract_feature_matrix(valid_smiles)
        if feature_matrix is not None:
            report['features'] = {
                'n_features': feature_matrix.shape[1],
                'feature_names': feature_names[:10],  # 只显示前10个
                'feature_matrix_shape': feature_matrix.shape
            }

        # 聚类分析
        clustering = self.perform_clustering(valid_smiles, method='kmeans', n_clusters=3)
        if 'error' not in clustering:
            report['clustering'] = {
                'method': clustering['method'],
                'n_clusters': clustering['n_clusters'],
                'silhouette_score': clustering.get('silhouette_score'),
                'cluster_sizes': {k: len(v) for k, v in clustering['clusters'].items()}
            }

        # 多样性分析
        diversity = self.analyze_isomer_diversity(valid_smiles)
        if 'error' not in diversity:
            report['diversity'] = {
                'diversity_score': diversity['diversity_score'],
                'mean_pairwise_distance': diversity['mean_pairwise_distance']
            }

        # 异常检测
        outliers = self.detect_outliers(valid_smiles, contamination=0.1)
        if 'error' not in outliers:
            report['outliers'] = {
                'n_outliers': outliers['n_outliers'],
                'outlier_ratio': outliers['outlier_ratio']
            }

        return report


# 便捷函数
def analyze_isomers(smiles_list: str | list, analysis_type: str = 'all') -> dict:
    """
    便捷的分析函数

    Args:
        smiles_list: SMILES字符串或字符串列表
        analysis_type: 分析类型 ('features', 'pca', 'clustering', 'diversity', 'outliers', 'all')

    Returns:
        分析结果
    """
    analyzer = MLIsomerAnalyzer()

    if isinstance(smiles_list, str):
        smiles_list = [smiles_list]

    if analysis_type == 'all':
        return analyzer.generate_analysis_report(smiles_list)
    elif analysis_type == 'features':
        feature_matrix, feature_names = analyzer.extract_feature_matrix(smiles_list)
        return {'feature_matrix': feature_matrix, 'feature_names': feature_names}
    elif analysis_type == 'pca':
        return analyzer.perform_pca_analysis(smiles_list)
    elif analysis_type == 'clustering':
        return analyzer.perform_clustering(smiles_list)
    elif analysis_type == 'diversity':
        return analyzer.analyze_isomer_diversity(smiles_list)
    elif analysis_type == 'outliers':
        return analyzer.detect_outliers(smiles_list)
    else:
        return {'error': f'不支持的分析类型: {analysis_type}'}


# 测试代码
if __name__ == '__main__':
    # 测试数据
    test_smiles = [
        'CCCCCCCC',  # 正辛烷
        'CCCC(C)CCC',  # 2-甲基庚烷
        'CC(C)CCCCCC',  # 2-甲基辛烷
        'CCCC(C)CCCC',  # 3-甲基庚烷
        'CC(C)CC(C)CC',  # 2,4-二甲基己烷
        'CC(C)CCC(C)C',  # 2,4-二甲基己烷
    ]

    # 创建分析器
    analyzer = MLIsomerAnalyzer()

    # 提取特征
    print("=== 提取分子特征 ===")
    for smiles in test_smiles:
        features = analyzer.extract_molecular_features(smiles)
        print(f"{smiles}: {features['mol_weight']:.2f} Da, LogP: {features['logp']:.2f}")

    # 聚类分析
    print("\n=== 聚类分析 ===")
    clustering = analyzer.perform_clustering(test_smiles, method='kmeans', n_clusters=2)
    print(f"聚类数量: {clustering['n_clusters']}")
    for cluster, molecules in clustering['clusters'].items():
        print(f"{cluster}: {len(molecules)} molecules")

    # 多样性分析
    print("\n=== 多样性分析 ===")
    diversity = analyzer.analyze_isomer_diversity(test_smiles)
    print(f"多样性得分: {diversity['diversity_score']:.4f}")
    print(f"平均距离: {diversity['mean_pairwise_distance']:.4f}")

    # 生成综合报告
    print("\n=== 综合分析报告 ===")
    report = analyzer.generate_analysis_report(test_smiles)
    print(f"有效分子数: {report['summary']['valid_molecules']}")
    print(f"聚类数量: {report['clustering'].get('n_clusters', 'N/A')}")
    print(f"异常值数量: {report['outliers'].get('n_outliers', 'N/A')}")

    # 可视化聚类
    if HAS_MATPLOTLIB:
        print("\n=== 生成聚类可视化 ===")
        viz_path = analyzer.visualize_clusters(test_smiles, method='kmeans', n_clusters=2)
        print(f"可视化已保存到: {viz_path}")
    else:
        print("\n提示: 安装matplotlib可启用可视化功能")

    # 检查依赖
    print("\n=== 依赖检查 ===")
    print(f"scikit-learn: {'✅' if HAS_SKLEARN else '❌'}")
    print(f"matplotlib: {'✅' if HAS_MATPLOTLIB else '❌'}")
    print(f"PIL: {'✅' if HAS_PIL else '❌'}")
