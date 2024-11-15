import os
from multiprocessing import Pool
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm

def dhs2gene(params):
    r"""Calculate the TR regulatory potential score.

    Parameters:
        args : argparse.Namespace
            Global parameters.
        sample : str
            TR sample name.
        vec : np.array
            TR PRE score.

    Returns:
        TR sample name, and TR-RP score.
    """
    
    args,sample,vec = params
    try:
        r = args.range
        if args.mean:
            sample_name = sample.split("@")[0]
        else:
            sample_name = sample
        e = args.TR_info.loc[sample_name,'Distal Intergenic Percentage']
        if e > 0.01:
            m = e
        else:
            m = 0.01
        alpha = np.log(2 / m - 0.99) / (r - 10000)
        print("e: %s, alpha: %s, sample: %s" % (e,alpha,sample))
        rp_vec = []
        for gene in args.genes:
            index = args.dhs_gene_g.loc[gene,'index']
            s = vec[index]
            d = args.dhs_gene_g.loc[gene,'DISTANCE']
            w = np.ones(d.shape)
            w[d > 10000] = 3 / (2 * np.exp((d[d > 10000] - 10000) * alpha) + 1)
            w[d > r] = 0
            if args.norm:
                w_e = 3 / (2 * np.exp((r - 10000) * alpha) + 1)
                w = (w - w_e) / (1 - w_e)
            rp = np.mean(np.multiply(w,s))
            rp_vec.append(rp)
        rp_vec = sparse.csr_matrix(rp_vec)
        return sample,rp_vec
    except:
        print('Error %s !' % sample)
        return None

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", type=str, default="library")
    parser.add_argument("--output", type=str, default="library")
    parser.add_argument("--hg38_refseq", type=str, default="library/hg38_refseq.ucsc")
    parser.add_argument("--dhs_hg38_rose", type=str, default="library/dhs_hg38_rose.bed")
    parser.add_argument("--DHS_TO_GENE", type=str, default="library/dhs_hg38_rose_DHS_TO_GENE.txt")
    parser.add_argument("--range", type=int, default=100000)
    parser.add_argument("--mean", type=str2bool, default=True)
    parser.add_argument("--norm", type=str2bool, default=False)
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
    TR_info = pd.read_csv(os.path.join(args.library,'TRs_info.txt'),sep='\t',index_col=0)
    if args.mean:
        args.TR_info = TR_info.groupby("tr_base").agg({"Distal Intergenic Percentage":np.mean}).copy()
    tr_dhs_ad = ad.read_h5ad(os.path.join(args.library,'TR_DHS.h5ad'))

    def iter_params():
        trunk_size = 100
        trunk_count = int(tr_dhs_ad.shape[0] / trunk_size) + 1
        for trunk in tqdm(range(trunk_count)):
            start = trunk * trunk_size
            end = (trunk + 1) * trunk_size
            tr_rp = tr_dhs_ad[start:end].to_df().values
            obs = tr_dhs_ad.obs.iloc[start:end]
            for i in range(tr_rp.shape[0]):
                tr_vec = tr_rp[i]
                sample_name = obs.index[i]
                yield args, sample_name, tr_vec
    
    rp_matrix = []
    sample_list = []
    with Pool(32) as pool:
        for row in tqdm(pool.imap(dhs2gene, iter_params()),total=tr_dhs_ad.shape[0]):
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
    obs = obs.join(TR_info,how='left')
    obs['index'] = range(len(obs))

    rp_matrix_ad.obs = obs
    rp_matrix_ad.write_h5ad(os.path.join(args.output,'RP_Matrix_TR.h5ad'))