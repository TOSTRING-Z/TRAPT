from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from numba import jit
from tqdm import tqdm

class CalcTRAUC:    
    r"""Calculate the area under the curve (AUC) for each TR curve.

    Attributes
    ----------
    args : TRAPT.Tools.Args
        Global parameters for TRAPT.
    RP_Matrix_TR_Sample : anndata.AnnData
        The sum of TR-RP scores and D-RP scores.
    w : np.array 
        U-RP scores.
    
    Notes
    -----
    The input is the RP matrix, and the calculation is performed as follows:
    
    .. math:: IRP = (TRRP + DRP) \times URP
    """
    def __init__(self, args, RP_Matrix_TR_Sample, w):
        self.args = args
        self.matrix = RP_Matrix_TR_Sample
        self.matrix.X *= w
        self.genes = self.matrix.var_names
        self.geneset = pd.read_csv(self.args.input, header=None)[0]

    @staticmethod
    @jit(nopython=True, nogil=True)
    def get_auc(params):
        r"""Parallel computing module.

        Parameters:
        i : int
            The i-th row of the I-RP matrix.
        j : int
            Default is 0.
        labels : np.array
            Gene vector.
        vec : np.array
            I-RP vector.

        Returns:
            i, j, and the AUC score of the i-th TR.
        """
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
        r"""Parallel parameter module.

        Parameters:
        gene_vec : np.array
            Gene vector.
        trunk : int
            Number of blocks.

        Returns:
            An iterator.
        """
        start = trunk * self.args.trunk_size
        end = (trunk + 1) * self.args.trunk_size
        tr_rp = self.matrix[start:end].to_df().values
        for i in range(tr_rp.shape[0]):
            tr_vec = tr_rp[i]
            yield start + i, 0, gene_vec, tr_vec

    def run(self):
        r"""TR auc calculation module execution entry point.

        Returns:
            A pd.DataFrame of AUC scores for TRs.
        """
        gene_vec = np.in1d(self.genes, self.geneset)
        print("Calculate the AUC...")
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
