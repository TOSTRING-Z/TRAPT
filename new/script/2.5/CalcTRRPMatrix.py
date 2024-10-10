import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--library", type=str, default="library")
parser.add_argument("--output", type=str, default="library")
args = parser.parse_args()

import re
from multiprocessing import Pool
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm

class CalcTRRPMatrix():
    def __init__(self,library='library',output='library'):
        self.library = library
        self.output = output

    def dhs2gene(self,sample):
        try:
            r = 100000
            d = 10000
            e = self.dip[sample]
            if e > 0:
                m = e
            else:
                m = 0.01
            alpha = np.log(2/m-1)/(r-d)
            return sample,sample.split("@")[0],alpha
        except:
            print('Error %s !' % sample)
            return None

    def run(self):
        TR_info = pd.read_csv(os.path.join(self.library,'TRs_info.txt'),sep='\t',index_col=0)
        dip = TR_info['Distal Intergenic Percentage']
        self.dip = dict(zip(dip.index,dip))
        tr_dhs_ad = h5py.File(os.path.join(self.library,'TR_DHS.h5ad'))
        samples = np.array(tr_dhs_ad['obs']['tr'],dtype=str)
        data = [self.dhs2gene(sample) for sample in tqdm(samples)]
        obs = pd.DataFrame(data,columns=["sample","TR","decay rate"]).set_index("sample",drop=True)
        obs.to_csv(f"{self.output}/decay_rate.txt",sep="\t")

if __name__ == '__main__':
    CalcTRRPMatrix(library=args.library,output=args.output).run()
