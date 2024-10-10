from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from numba import jit
from tqdm import tqdm


class CalcTRAUC:
    def __init__(self, args, RP_Matrix_TR_Sample, w):
        self.args = args
        self.matrix = RP_Matrix_TR_Sample
        self.matrix.X *= w
        self.genes = self.matrix.var_names
        self.geneset = pd.read_csv(self.args.input, header=None)[0]

    @staticmethod
    @jit(nopython=True, nogil=True)
    def get_auc(params):
        i, j, labels, vec = params
        l_p = labels > 0.5
        p_c = np.sum(l_p)
        n_c = len(l_p) - p_c
        if p_c > 1 and n_c > 1:
            index = np.argsort(vec)
            l_p_rank = np.where(l_p[index])[0]
            rank_sum = np.sum(l_p_rank)
            auc = (rank_sum - p_c * (1 + p_c) / 2) / (p_c * n_c)
        else:
            auc = 0
        return i, j, auc

    def iter_params(self, gene_vec, trunk):
        start = trunk * self.args.trunk_size
        end = (trunk + 1) * self.args.trunk_size
        tr_rp = self.matrix[start:end].to_df().values
        for i in range(tr_rp.shape[0]):
            tr_vec = tr_rp[i]
            yield start + i, 0, gene_vec, tr_vec

    def run(self):
        gene_vec = np.in1d(self.genes, self.geneset)

        auc = np.zeros((self.matrix.shape[0], 1))
        with ThreadPoolExecutor(self.args.threads) as pool:
            trunk_count = int(self.matrix.shape[0] / self.args.trunk_size) + 1
            for trunk in tqdm(range(trunk_count)):
                data = self.iter_params(gene_vec, trunk)
                tasks = [pool.submit(self.get_auc, params) for params in data]
                for task in tasks:
                    i, j, score = task.result()
                    auc[i, j] = score

        auc = pd.DataFrame(auc, index=self.matrix.obs_names)
        return auc
