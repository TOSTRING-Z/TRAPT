import os
import re
import sys

from multiprocessing import Pool

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm


class CalcTRRPMatrix():
    def __init__(self,library='library',output='library'):
        self.library = library
        self.output = output
        self.hg38_refseq = os.path.join(library,'hg38_refseq.ucsc')
        self.dhs_hg38_rose = os.path.join(library,'dhs_hg38_rose.bed')
        self.DHS_TO_GENE = os.path.join(library,'dhs_hg38_rose_DHS_TO_GENE.txt')
    
    def dhs2gene(self,sample):
        id = re.findall('dhs@(.*?).csv',sample)[0]
        try:
            r = 100000
            d_ = 10000
            e = self.TR_info.loc[id,'Distal Intergenic Percentage']
            if e > 0:
                m = e
            else:
                m = 0.01
            alpha = (r-d_)/np.log(2/m-0.99)
            sample_name = re.findall(r'dhs@(.*?).csv',sample)[0]
            vec = pd.read_csv(sample,sep='\t',header=None)[0].values
            rp_vec = []
            for gene in self.genes:
                index = self.dhs_gene_g.loc[gene,'index']
                s = vec[index]
                d = self.dhs_gene_g.loc[gene,'DISTANCE']
                w = np.ones(d.shape)
                w[d > d_] = 2/(np.exp((d[d > d_]-d_)/alpha)+1)
                w[d > r] = 0
                rp = np.mean(np.multiply(w,s))
                rp_vec.append(rp)
            rp_vec = sparse.csr_matrix(rp_vec)
            return sample_name,rp_vec
        except:
            print('Error %s !' % sample)
            return None

    def run(self):
        ucsc_gene = pd.read_csv(self.hg38_refseq,sep='\t')
        ucsc_gene = ucsc_gene[['name','name2']]
        ucsc_gene.columns = ['GENES','SYMBOLS']
        ucsc_gene = ucsc_gene.drop_duplicates()
        dhs_hg38 = pd.read_csv(self.dhs_hg38_rose,sep='\t',header=None)
        dhs_hg38 = dhs_hg38.reset_index()[['index',0]]
        dhs_hg38.columns = ['index','DHS']
        dhs_gene = pd.read_csv(self.DHS_TO_GENE,sep='\t')
        dhs_gene = dhs_gene.iloc[:,[0,4,5]]
        dhs_gene.columns = ['DHS','GENES','DISTANCE']
        dhs_gene_merge = dhs_gene.merge(dhs_hg38,on='DHS')
        dhs_gene_merge = dhs_gene_merge.merge(ucsc_gene,on='GENES')
        dhs_gene_merge = dhs_gene_merge.drop_duplicates(['DHS','SYMBOLS'])
        self.dhs_gene_g = dhs_gene_merge.groupby('SYMBOLS').agg({
            'DISTANCE':list,
            'index': list
        })
        self.dhs_gene_g['DISTANCE'] = self.dhs_gene_g['DISTANCE'].apply(lambda x:np.array(x))
        self.genes = self.dhs_gene_g.index
        self.TR_info = pd.read_csv(os.path.join(self.library,'TRs_info.txt'),sep='\t',index_col=0)
        # New data
        # from glob import glob
        # samples = golb(os.path.join(self.library,'TR_matrix','dhs@*.csv'))
        tr_dhs_ad = h5py.File(os.path.join(self.library,'TR_DHS.h5ad'))
        samples = np.array(tr_dhs_ad['obs']['tr'],dtype=str)
        samples = list(map(lambda x:os.path.join(self.library,'TR_matrix',f'dhs@{x}.csv'),samples))

        args = list(map(lambda sample:(sample),samples))
        rp_matrix = []
        sample_list = []
        with Pool(16) as pool:
            for row in tqdm(pool.imap(self.dhs2gene, args),total=len(args)):
                if row:
                    sample_name,rp_vec = row
                    sample_list.append(sample_name)
                    rp_matrix.append(rp_vec)

        rp_matrix = sparse.vstack(rp_matrix,dtype='float32')
        rp_matrix_ad = ad.AnnData(rp_matrix)
        rp_matrix_ad.var_names = self.genes
        rp_matrix_ad.obs_names = sample_list

        obs = rp_matrix_ad.obs
        obs.index.name = 'tr'
        obs = obs.join(self.TR_info,how='left')
        obs['index'] = range(len(obs))

        rp_matrix_ad.obs = obs
        rp_matrix_ad.write_h5ad(os.path.join(self.output,f'RP_Matrix_TR.h5ad'))

if __name__ == '__main__':
    library = sys.argv[1]
    CalcTRRPMatrix(library=library,output=library).run()
