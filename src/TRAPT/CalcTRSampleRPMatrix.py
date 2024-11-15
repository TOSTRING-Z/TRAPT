import os
import sys

import anndata as ad
import numpy as np


class CalcTRSampleRPMatrix:
    r"""The D-RP model network aggregation module.

    Attributes
    ----------
    library : str
        Path to the background library.
    output : str
        Output path, default is 'library'.
    type : str
        H3K27ac/ATAC.
    """
    def __init__(self, library='library', output='library', type='H3K27ac'):
        self.data_ad_sample = ad.read_h5ad(
            os.path.join(library, f'RP_Matrix_{type}.h5ad')
        )
        self.data_ad_TR = ad.read_h5ad(os.path.join(library, 'RP_Matrix_TR.h5ad'))
        self.A = ad.read_h5ad(os.path.join(library, f'A_pred_{type}.h5ad'))
        self.A.X[np.diag_indices_from(self.A.X)] = 0
        i = self.data_ad_TR.shape[0]
        self.A = self.A[:i, i:]
        self.library = library
        self.output = output
        self.type = type

    def run(self):
        r"""D-RP model network aggregation module execution entry point.
        """
        self.A = self.A.to_df().values
        self.A[self.A < 0.0] = 0
        self.X = self.data_ad_sample.to_df().values
        self.A = self.A / (np.sum(self.A, axis=1, keepdims=True) + 1e-17)
        self.X = self.X / (np.linalg.norm(self.X, axis=1, keepdims=True) + 1e-17)
        matrix = self.A.dot(self.X)
        data_ad = ad.AnnData(
            matrix, obs=self.data_ad_TR.obs, var=self.data_ad_TR.var, dtype='float32'
        )
        data_ad.write_h5ad(os.path.join(self.library, f'RP_Matrix_TR_{self.type}.h5ad'))


if __name__ == '__main__':
    type = sys.argv[1]
    library = sys.argv[2]
    assert type in ['H3K27ac', 'ATAC']
    print(f'\033[1;31m # # # CalcTRSampleRPMatrix # # #\033[0m')
    CalcTRSampleRPMatrix(library=library, output=library, type=type).run()
