import numpy as np
from scipy.special import logsumexp
from sklearn.cluster import KMeans

def bernoulli_mixture_model_clustering(mutation_vectors, 
                                        max_clusters=10, 
                                        num_runs=10, 
                                        max_iterations=300, 
                                        conv_thresh=1e-4,
                                        alpha=1.01, 
                                        beta=2):
    """
    基于 Bernoulli Mixture Model 的聚类算法，使用 EM 算法估计参数。

    Args:
        mutation_vectors: 突变矩阵，形状为 (N, D)，其中 N 是 reads 数量，D 是位点数量。
        max_clusters: 最大簇数。
        num_runs: 运行 EM 算法的次数。
        max_iterations: EM 算法的最大迭代次数。
        conv_thresh: EM 算法的收敛阈值。
        alpha: Dirichlet 先验参数 alpha，用于平滑 mu。
        beta: Dirichlet 先验参数 beta，用于平滑 mu。

    Returns:
        best_cluster_labels: 每个 reads 所属的簇标签。
        best_cluster_mutation_rates: 每个簇的突变概率向量。
    """

    N, D = mutation_vectors.shape
    best_log_likelihood = -np.inf
    best_cluster_labels = None
    best_cluster_mutation_rates = None

    for K in range(1, max_clusters + 1):
        for _ in range(num_runs):
            # 初始化参数
            pi, mu = initialize_parameters(mutation_vectors, K)

            # EM 算法迭代
            for i in range(max_iterations):
                # E 步：计算后验概率
                responsibilities = e_step(mutation_vectors, pi, mu)

                # M 步：更新模型参数
                pi, mu = m_step(mutation_vectors, responsibilities, alpha, beta)

                # 计算对数似然函数值
                log_likelihood_value = calculate_log_likelihood(mutation_vectors, pi, mu)

                # 检查收敛条件
                if i > 0 and abs(log_likelihood_value - prev_log_likelihood) < conv_thresh:
                    break

                prev_log_likelihood = log_likelihood_value

            # 更新最佳模型
            if log_likelihood_value > best_log_likelihood:
                best_log_likelihood = log_likelihood_value
                best_cluster_labels = np.argmax(responsibilities, axis=1)
                best_cluster_mutation_rates = mu

    return best_cluster_labels, best_cluster_mutation_rates


def initialize_parameters(mutation_vectors, K):
    """使用 KMeans 初始化模型参数。"""
    #kmeans = KMeans(n_clusters=K, random_state=0).fit(mutation_vectors)
    kmeans = KMeans(n_clusters=K, n_init=10, random_state=0).fit(mutation_vectors)
    pi = np.ones(K) / K
    mu = kmeans.cluster_centers_
    return pi, mu

def e_step(mutation_vectors, pi, mu):
    """EM 算法的 E 步：计算每个 read 属于每个簇的后验概率。"""
    N, D = mutation_vectors.shape
    K = pi.shape[0]

    responsibilities = np.zeros((N, K))
    for n in range(N):
        for k in range(K):
            # 逐个位点计算 log 概率，避免数值下溢
            for d in range(D):
                # responsibilities[n, k] += np.log(
                #     (mu[k, d] ** mutation_vectors[n, d]) * ((1 - mu[k, d]) ** (1 - mutation_vectors[n, d]))
                # )
                
                epsilon = 1e-10  # 选择一个非常小的值
                # responsibilities[n, k] += np.log(
                #     np.maximum(mu[k, d] ** mutation_vectors[n, d], epsilon) *
                #     np.maximum((1 - mu[k, d]) ** (1 - mutation_vectors[n, d]), epsilon)
                # )
                responsibilities[n, k] += mutation_vectors[n, d] * np.log(mu[k, d] + epsilon) + \
                               (1 - mutation_vectors[n, d]) * np.log(1 - mu[k, d] + epsilon)

            responsibilities[n, k] += np.log(pi[k])
        # 使用 logsumexp 进行数值稳定的计算
        responsibilities[n, :] -= logsumexp(responsibilities[n, :]) 
    responsibilities = np.exp(responsibilities)  # 将 log 概率转换回概率
    return responsibilities


def m_step(mutation_vectors, responsibilities, alpha, beta):
    """EM 算法的 M 步，更新模型参数 pi 和 mu。"""
    N, D = mutation_vectors.shape
    K = responsibilities.shape[1]

    # 更新 pi
    N_k = np.sum(responsibilities, axis=0)
    pi = N_k / N

    # 更新 mu
    mu = np.zeros((K, D))
    for k in range(K):
        for d in range(D):
            mu[k, d] = (np.sum(responsibilities[:, k] * mutation_vectors[:, d]) + alpha - 1) / (
                N_k[k] + alpha + beta - 2
            )
    return pi, mu

def calculate_log_likelihood(mutation_vectors, pi, mu):
    """计算对数似然函数值"""
    N, D = mutation_vectors.shape
    K = pi.shape[0]

    log_likelihood_value = 0
    for n in range(N):
        temp_likelihood = 0
        for k in range(K):
            temp_likelihood += pi[k] * np.prod((mu[k] ** mutation_vectors[n]) * ((1 - mu[k]) ** (1 - mutation_vectors[n])))
        log_likelihood_value += np.log(temp_likelihood)

    return log_likelihood_value

import sys
#filename = "/disk1/work/zehui/gpu_public/public_data_result/split_structure/data/filter_model_reads_fake_wide.txt"
filename = sys.argv[1]
out_filename = sys.argv[2]

mutation_vectors = []
name_vectors = []
head = 1
for line in open(filename):
    if head > 0 :
        head -= 1
        continue
    lineL = line.strip().split("\t")
    vectors = [int(n) for n in lineL[1:]]
    mutation_vectors.append(vectors)
    name = lineL[0]
    name_vectors.append(name)
#
mutation_vectors = np.array(mutation_vectors)
# 运行聚类算法
cluster_labels, cluster_mutation_rates = bernoulli_mixture_model_clustering(
    mutation_vectors, 
    max_clusters=2,
    num_runs=10,
    max_iterations=75,
)

#print("Cluster Labels:", cluster_labels)
#print("Cluster Mutation Rates:\n", cluster_mutation_rates)
fh = open(out_filename,"w")
for n in range(len(cluster_labels)):
    try:
        fh.write(f"{name_vectors[n]}\t{cluster_labels[n]}\n")
        #print(f"{name_vectors[n]}\t{cluster_labels[n]}")
    except:
        continue
fh.close()

