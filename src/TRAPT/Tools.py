import os

import anndata as ad
import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler

class Args:
    r""" TRAPT Global Parameters.

    Attributes
    ----------
    input : str
        Input path for the gene set.
    output : str
        Output path for TRAPT results.
    library : str, optional
        Path to the background library, default is the 'library' path in the current directory.
    threads : int, optional
        Number of processes used for TRAPT inference.
    trunk_size : int, optional
        Size of the chunks.
    background_genes : str, optional
        Number of background genes selected.
    use_kd : str, optional
        Use knowledge distillation.
    tr_type : str, optional
        all/tf/tcof/cr.
    source : str, optional
        all/cistrome/chip_altas/gtrd/remap/chip-atlas/remap/encode/geo.
    """
    def __init__(self, input, output, library="library", 
                 threads=16, trunk_size=2048 * 16, background_genes=6000, 
                 use_kd=True, tr_type="all", source="all") -> None:
        self.input = input
        self.output = output
        self.library = library
        self.threads = threads
        self.trunk_size = trunk_size
        self.background_genes = background_genes
        self.use_kd = use_kd
        self.tr_type = tr_type
        self.source = source
        self.liner = True
        if not os.path.exists(output):
            os.mkdir(output)

class Type:
    H3K27ac = 'H3K27ac'
    ATAC = 'ATAC'


class RPMatrix:
    def __init__(self, library, name, to_array=True):
        self.data = ad.read_h5ad(os.path.join(library, name))
        if to_array:
            self.data.X = self.data.to_df().values

    def norm(self, type="l2", axis=1):
        assert type in ["l1", "l2"]
        self.data.X = self.data.X - self.data.X.min()
        if type == "l1":
            self.data.X /= np.sum(self.data.X, axis=axis, keepdims=True).clip(min=1e-17)
        if type == "l2":
            self.data.X /= np.linalg.norm(self.data.X, axis=axis, keepdims=True).clip(
                min=1e-17
            )
        return self

    def standard_scale(self, axis=1):
        ss = StandardScaler()
        if axis == 0:
            self.data.X = ss.fit_transform(self.data.X)
        if axis == 1:
            self.data.X = ss.fit_transform(self.data.X.T).T
        return self

    def binarization(self):
        self.data.X = (self.data.X > 1e-17).astype(np.float32)
        return self

    def minmax_scale(self, axis=1):
        ss = MinMaxScaler()
        if axis == 0:
            self.data.X = ss.fit_transform(self.data.X)
        if axis == 1:
            self.data.X = ss.fit_transform(self.data.X.T).T
        return self

    def add(self, data:ad.AnnData):
        self.data.X += data.X
        return self

    def get_data(self) -> ad.AnnData:
        return self.data


class RP_Matrix:
    def __init__(self, library) -> None:
        self.TR = (
            RPMatrix(library, 'RP_Matrix_TR.h5ad')
            .norm().get_data()
        )
        self.TR_H3K27ac = (
            RPMatrix(library, 'RP_Matrix_TR_H3K27ac.h5ad')
            .norm().add(self.TR).get_data()
        )
        self.TR_ATAC = (
            RPMatrix(library, 'RP_Matrix_TR_ATAC.h5ad')
            .norm().add(self.TR).get_data()
        )
        self.H3K27ac = (
            RPMatrix(library, 'RP_Matrix_H3K27ac.h5ad')
            .standard_scale().get_data()
        )
        self.ATAC = (
            RPMatrix(library, 'RP_Matrix_ATAC.h5ad')
            .standard_scale().get_data()
        )

