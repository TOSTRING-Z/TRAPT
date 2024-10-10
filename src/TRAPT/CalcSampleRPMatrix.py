import os
import re
import sys

from glob import glob
from multiprocessing import Pool

import anndata
import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm


class CalcSampleRPMatrix:
    def __init__(self, library='library', output='library', type='H3K27ac'):
        self.library = library
        self.output = output
        self.type = type
        self.hg38_refseq = os.path.join(library, 'hg38_refseq.ucsc')
        self.dhs_hg38_rose = os.path.join(library, 'dhs_hg38_rose.bed')
        self.DHS_TO_GENE = os.path.join(library, 'dhs_hg38_rose_DHS_TO_GENE.txt')

    def dhs2gene(self, sample):
        try:
            r = 100000
            d_ = 10000
            m = 0.01
            alpha = (r - d_) / np.log(2 / m - 1)
            sample_name = re.findall(r'dhs@(.*?).csv', sample)[0]
            vec = pd.read_csv(sample, sep='\t', header=None)[0].values
            rp_vec = []
            for gene in self.genes:
                index = self.dhs_gene_g.loc[gene, 'index']
                s = vec[index]
                d = self.dhs_gene_g.loc[gene, 'DISTANCE']
                w = np.ones(d.shape)
                w[d > d_] = 2 / (np.exp((d[d > d_] - d_) / alpha) + 1)
                w[d > r] = 0
                rp = np.mean(np.multiply(w, s))
                rp_vec.append(rp)
            rp_vec = np.log(np.array(rp_vec) + 1)
            rp_vec = sparse.csr_matrix(rp_vec)
            return sample_name, rp_vec
        except:
            print('Error %s !' % sample)
            return None

    def run(self):
        ucsc_gene = pd.read_csv(self.hg38_refseq, sep='\t')
        ucsc_gene = ucsc_gene[['name', 'name2']]
        ucsc_gene.columns = ['GENES', 'SYMBOLS']
        ucsc_gene = ucsc_gene.drop_duplicates()
        dhs_hg38 = pd.read_csv(self.dhs_hg38_rose, sep='\t', header=None)
        dhs_hg38 = dhs_hg38.reset_index()[['index', 0]]
        dhs_hg38.columns = ['index', 'DHS']
        dhs_gene = pd.read_csv(self.DHS_TO_GENE, sep='\t')
        dhs_gene = dhs_gene.iloc[:, [0, 4, 5]]
        dhs_gene.columns = ['DHS', 'GENES', 'DISTANCE']
        dhs_gene_merge = dhs_gene.merge(dhs_hg38, on='DHS')
        dhs_gene_merge = dhs_gene_merge.merge(ucsc_gene, on='GENES')
        dhs_gene_merge = dhs_gene_merge.drop_duplicates(['DHS', 'SYMBOLS'])
        self.dhs_gene_g = dhs_gene_merge.groupby('SYMBOLS').agg(
            {'DISTANCE': list, 'index': list}
        )
        self.dhs_gene_g['DISTANCE'] = self.dhs_gene_g['DISTANCE'].apply(
            lambda x: np.array(x)
        )
        self.genes = self.dhs_gene_g.index
        rp_matrix = []
        samples = sorted(
            glob(os.path.join(self.library, f'{self.type}_matrix', 'dhs@*_hg38.csv'))
        )
        args = list(map(lambda sample: (sample), samples))
        rp_matrix = []
        sample_list = []
        with Pool(16) as pool:
            for row in tqdm(pool.imap(self.dhs2gene, args), total=len(args)):
                if row:
                    sample_name, rp_vec = row
                    sample_list.append(sample_name)
                    rp_matrix.append(rp_vec)
        rp_matrix = sparse.vstack(rp_matrix, dtype='float32')
        rp_matrix_ad = anndata.AnnData(rp_matrix)
        rp_matrix_ad.var_names = self.genes
        rp_matrix_ad.obs_names = sample_list
        rp_matrix_ad.write_h5ad(
            os.path.join(self.output, f'RP_Matrix_{self.type}.h5ad')
        )


if __name__ == '__main__':
    type = sys.argv[1]
    library = sys.argv[2]
    assert type in ['H3K27ac', 'ATAC']
    CalcSampleRPMatrix(library=library, output=library, type=type).run()
