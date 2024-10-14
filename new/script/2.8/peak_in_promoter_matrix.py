import os

from multiprocessing import Pool

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm


def dhs2gene(params):
    args,sample,vec = params
    try:
        r = 10000
        rp_vec = []
        vec = vec.toarray()[:,0]
        for gene in args.genes:
            index = args.dhs_gene_g.loc[gene,'index']
            s = vec[index]
            d = args.dhs_gene_g.loc[gene,'DISTANCE']
            w = np.ones(d.shape)
            w[d > r] = 0
            rp = int(np.mean(np.multiply(w,s)) > 0)
            rp_vec.append(rp)
        rp_vec = sparse.csr_matrix(rp_vec)
        return sample,rp_vec
    except:
        print('Error %s !' % sample)
        return None

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=str, default="library")
    parser.add_argument("--output", type=str, default="new/result/2.8")
    parser.add_argument("--hg38_refseq", type=str, default="library/hg38_refseq.ucsc")
    parser.add_argument("--dhs_hg38_rose", type=str, default="library/dhs_hg38_rose.bed")
    parser.add_argument("--DHS_TO_GENE", type=str, default="library/dhs_hg38_rose_DHS_TO_GENE.txt")
    args = parser.parse_args()

    ucsc_gene = pd.read_csv(args.hg38_refseq,sep='\t')
    ucsc_gene = ucsc_gene[['name','name2']]
    ucsc_gene.columns = ['GENES','SYMBOLS']
    ucsc_gene = ucsc_gene.drop_duplicates()
    dhs_hg38 = pd.read_csv(args.dhs_hg38_rose,sep='\t',header=None)
    dhs_hg38 = dhs_hg38.reset_index()[['index',0]]
    dhs_hg38.columns = ['index','DHS']
    dhs_gene = pd.read_csv(args.DHS_TO_GENE,sep='\t')
    dhs_gene = dhs_gene.iloc[:,[0,4,5]]
    dhs_gene.columns = ['DHS','GENES','DISTANCE']
    dhs_gene_merge = dhs_gene.merge(dhs_hg38,on='DHS')
    dhs_gene_merge = dhs_gene_merge.merge(ucsc_gene,on='GENES')
    dhs_gene_merge = dhs_gene_merge.drop_duplicates(['DHS','SYMBOLS'])
    args.dhs_gene_g = dhs_gene_merge.groupby('SYMBOLS').agg({
        'DISTANCE':list,
        'index': list
    })
    args.dhs_gene_g['DISTANCE'] = args.dhs_gene_g['DISTANCE'].apply(lambda x:np.array(x))
    args.genes = args.dhs_gene_g.index
    args.TR_info = pd.read_csv(os.path.join(args.library,'TRs_info.txt'),sep='\t',index_col=0)
    tr_dhs_ad = ad.read_h5ad(os.path.join(args.library,'TR_DHS.h5ad'))

    params = ((args,tr_dhs_ad.obs.index[i],tr_dhs_ad[i].X.T) for i in np.arange(tr_dhs_ad.shape[0]))
    rp_matrix = []
    sample_list = []
    with Pool(16) as pool:
        for row in tqdm(pool.imap(dhs2gene, params),total=tr_dhs_ad.shape[0]):
            if row:
                sample_name,rp_vec = row
                sample_list.append(sample_name)
                rp_matrix.append(rp_vec)

    rp_matrix = sparse.vstack(rp_matrix,dtype='float32')
    rp_matrix_ad = ad.AnnData(rp_matrix)
    rp_matrix_ad.var_names = args.genes
    rp_matrix_ad.obs_names = sample_list

    obs = rp_matrix_ad.obs
    obs.index.name = 'tr'
    obs = obs.join(args.TR_info,how='left')
    obs['index'] = range(len(obs))

    rp_matrix_ad.obs = obs
    rp_matrix_ad.write_h5ad(os.path.join(args.output,f'Peak_Matrix_TR.h5ad'))