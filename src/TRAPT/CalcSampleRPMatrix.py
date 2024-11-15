import os
import re
from glob import glob
from multiprocessing import Pool
import anndata
import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm

def dhs2gene(args,sample):
    r"""Calculate the Epi regulatory potential score.

    Parameters:
        args : argparse.Namespace
            Global parameters.
        sample : str
            Epi sample name.
        vec : np.array
            Epi PRE score.

    Returns:
        Epi sample name, and Epi-RP score.
    """
    try:
        r = 100000
        m = 0.01
        alpha = np.log(2 / m - 0.99) / (r - 10000)
        sample_name = re.findall(r'dhs@(.*?).csv', sample)[0]
        vec = pd.read_csv(sample, sep='\t', header=None)[0].values
        rp_vec = []
        for gene in args.genes:
            index = args.dhs_gene_g.loc[gene, 'index']
            s = vec[index]
            d = args.dhs_gene_g.loc[gene, 'DISTANCE']
            w = np.ones(d.shape)
            w[d > 10000] = 3 / (2 * np.exp((d[d > 10000] - 10000) * alpha) + 1)
            w[d > r] = 0
            rp = np.mean(np.multiply(w, s))
            rp_vec.append(rp)
        rp_vec = np.log(np.array(rp_vec) + 1)
        rp_vec = sparse.csr_matrix(rp_vec)
        return sample_name, rp_vec
    except:
        print('Error %s !' % sample)
        return None

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=str, default="library")
    parser.add_argument("--output", type=str, default="library")
    parser.add_argument("--type", type=str, default="H3K27ac",help="H3K27ac/ATAC")
    parser.add_argument("--hg38_refseq", type=str, default="library/hg38_refseq.ucsc")
    parser.add_argument("--dhs_hg38_rose", type=str, default="library/dhs_hg38_rose.bed")
    parser.add_argument("--DHS_TO_GENE", type=str, default="library/dhs_hg38_rose_DHS_TO_GENE.txt")
    args = parser.parse_args()
    assert args.type in ['H3K27ac', 'ATAC']
    
    ucsc_gene = pd.read_csv(args.hg38_refseq, sep='\t')
    ucsc_gene = ucsc_gene[['name', 'name2']]
    ucsc_gene.columns = ['GENES', 'SYMBOLS']
    ucsc_gene = ucsc_gene.drop_duplicates()
    dhs_hg38 = pd.read_csv(args.dhs_hg38_rose, sep='\t', header=None)
    dhs_hg38 = dhs_hg38.reset_index()[['index', 0]]
    dhs_hg38.columns = ['index', 'DHS']
    dhs_gene = pd.read_csv(args.DHS_TO_GENE, sep='\t')
    dhs_gene = dhs_gene.iloc[:, [0, 4, 5]]
    dhs_gene.columns = ['DHS', 'GENES', 'DISTANCE']
    dhs_gene_merge = dhs_gene.merge(dhs_hg38, on='DHS')
    dhs_gene_merge = dhs_gene_merge.merge(ucsc_gene, on='GENES')
    dhs_gene_merge = dhs_gene_merge.drop_duplicates(['DHS', 'SYMBOLS'])
    args.dhs_gene_g = dhs_gene_merge.groupby('SYMBOLS').agg(
        {'DISTANCE': list, 'index': list}
    )
    args.dhs_gene_g['DISTANCE'] = args.dhs_gene_g['DISTANCE'].apply(
        lambda x: np.array(x)
    )
    args.genes = args.dhs_gene_g.index
    rp_matrix = []
    samples = sorted(
        glob(os.path.join(args.library, f'{args.type}_matrix', 'dhs@*_hg38.csv'))
    )
    params = list(map(lambda sample: (sample), samples))
    rp_matrix = []
    sample_list = []
    with Pool(16) as pool:
        for row in tqdm(pool.imap(args.dhs2gene, params), total=len(params)):
            if row:
                sample_name, rp_vec = row
                sample_list.append(sample_name)
                rp_matrix.append(rp_vec)
    rp_matrix = sparse.vstack(rp_matrix, dtype='float32')
    rp_matrix_ad = anndata.AnnData(rp_matrix)
    rp_matrix_ad.var_names = args.genes
    rp_matrix_ad.obs_names = sample_list
    rp_matrix_ad.write_h5ad(
        os.path.join(args.output, f'RP_Matrix_{args.type}.h5ad')
    )
