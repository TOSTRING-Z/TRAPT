import os
import random
import sys

import anndata as ad
import numpy as np
import pandas as pd
import torch
import torch.nn.modules.loss
from scipy import sparse
from sklearn.preprocessing import normalize
from torch import Tensor, nn
from torch.nn.functional import binary_cross_entropy, mse_loss
from torch.utils.data import DataLoader
import torch_geometric.transforms as T
from torch_geometric.data import Data
from torch_geometric.nn import VGAE, GCNConv, InnerProductDecoder
from tqdm import tqdm

from TRAPT.Tools import RPMatrix


def seed_torch(seed=2023):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)


class InnerProductDecoderWeight(InnerProductDecoder):
    def __init__(self, A_e, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.A_e = A_e

    def forward(
        self, z: Tensor, edge_index: Tensor = None, sigmoid: bool = True
    ) -> Tensor:
        if edge_index != None:
            adj = (z[edge_index[0]] * z[edge_index[1]]).sum(dim=1)
        else:
            adj = torch.matmul(z, z.t())
        return torch.sigmoid(adj) if sigmoid else adj


class VariationalGCNEncoder(torch.nn.Module):
    def __init__(self, in_channels, h_dim, z_dim):
        super(VariationalGCNEncoder, self).__init__()
        self.conv_h = GCNConv(in_channels, h_dim, cached=True)
        self.conv_mu = GCNConv(h_dim, z_dim, cached=True)
        self.conv_logstd = GCNConv(h_dim, z_dim, cached=True)

    def predict_h(self, x, edge_index):
        h = self.conv_h(x, edge_index).relu()
        return h

    def forward(self, x, edge_index):
        h = self.conv_h(x, edge_index).relu()
        return self.conv_mu(h, edge_index), self.conv_logstd(h, edge_index)


class CVAE(nn.Module):
    def __init__(self, input_dim, condition_dim, h_dim, z_dim):
        super(CVAE, self).__init__()
        self.h_encoder = torch.nn.Linear(input_dim + condition_dim, h_dim)
        self.act = torch.nn.ReLU()
        self.dropout = torch.nn.Dropout()
        self.layer_mu = torch.nn.Linear(h_dim, z_dim)
        self.layer_logstd = torch.nn.Linear(h_dim, z_dim)
        self.decoder = torch.nn.Linear(z_dim + condition_dim, input_dim)

    def reparametrize(self, mu: Tensor, logstd: Tensor) -> Tensor:
        if self.training:
            return mu + torch.randn_like(logstd) * torch.exp(logstd)
        else:
            return mu

    def predict_h(self, x):
        h = self.act(self.h_encoder(x))
        return h

    def kl_div(self):
        return torch.mean(
            0.5
            * torch.sum(
                torch.exp(self.logstd) + self.mu**2 - 1.0 - self.logstd, dim=1
            )
        )

    def forward(self, x):
        condition = x[:, -2:]
        h = self.act(self.h_encoder(x))
        self.mu = self.layer_mu(h)
        self.logstd = self.layer_logstd(h)
        z = self.reparametrize(self.mu, self.logstd)
        decoded = self.decoder(torch.concat([z, condition], dim=1))
        return decoded

class CalcSTM:
    r"""D-RP model network reconstruction module.

    Parameters:
    RP_Matrix : TRAPT.Tools.RP_Matrix
        TR-RP matrix and Epi-RP matrix.
    type : str
        Epi-RP type.
    checkpoint_path : str 
        Model save path.
    device : str, optional
        cpu/cuda.
    """
    def __init__(self, RP_Matrix, type, checkpoint_path, device='cuda'):
        self.type = type
        assert type in ['H3K27ac', 'ATAC']
        self.checkpoint_path = checkpoint_path
        self.device = device
        print(f"Using {self.device}")
        self.RP_Matrix = RP_Matrix
        self.data = torch.concat(
            [
                torch.tensor(self.RP_Matrix.TR.to_df().values),
                torch.tensor(self.RP_Matrix.Sample.to_df().values),
            ],
            dim=0,
        ).to(self.device)

    @staticmethod
    def get_cos_similar_matrix(m1, m2):
        r"""Matrix cosine similarity calculation.
        """
        num = m1.dot(m2.T)
        denom = np.linalg.norm(m1, axis=1).reshape(-1, 1) * np.linalg.norm(m2, axis=1)
        res = num / denom
        res[np.isneginf(res)] = 0
        return res

    def get_edge_index(self, A, B, n=10):
        r"""Construct a heterogeneous network.

        Parameters:
        A : anndata.AnnData
            TR-RP matrix.
        B : anndata.AnnData
            Epi-RP matrix.
        n : int 
            Number of nearest neighbors for TR.

        Returns:
            TR-Epi heterogeneous network.
        """
        self.A_s = 0
        self.A_e = A.shape[0]
        self.B_s = self.A_e
        self.B_e = self.B_s + B.shape[0]
        adj = np.zeros((self.B_e, self.B_e))
        AB_cos_sim = self.get_cos_similar_matrix(A.to_df().values, B.to_df().values)
        AB = (AB_cos_sim >= np.sort(AB_cos_sim, axis=1)[:, -n].reshape(-1, 1)).astype(
            int
        )
        adj[self.A_s : self.A_e, self.B_s : self.B_e] = AB
        adj[self.B_s : self.B_e, self.A_s : self.A_e] = AB.T
        row, col, value = sparse.find(adj)
        return torch.tensor(np.array([row, col]), dtype=torch.long)

    @staticmethod
    def sparse_to_tensor(data, type='sparse'):
        if type == 'dense':
            return torch.tensor(data.to_df().values)
        else:
            row, col, value = sparse.find(data)
            indices = torch.tensor(np.array([row, col]), dtype=torch.long)
            values = torch.tensor(value, dtype=torch.float)
            size = torch.Size(data.shape)
            return torch.sparse_coo_tensor(indices, values, size)

    def save_graph(self):
        x, edge_index = self.data.x.to(self.device), self.data.edge_index.to(
            self.device
        )
        z = self.model_vgae.encode(x, edge_index).cpu().detach().numpy()
        A_pred = z.dot(z.T)
        TR_info = pd.DataFrame(list(self.RP_Matrix.TR.obs_names))
        TR_info['type'] = 'TR'
        Sample_info = pd.DataFrame(list(self.RP_Matrix.Sample.obs_names))
        Sample_info['type'] = self.type
        info = pd.concat([TR_info, Sample_info]).astype(str)
        info.columns = ['index', 'type']
        info = info.set_index('index')
        data = ad.AnnData(A_pred, obs=info, var=info)
        data.write_h5ad(f'{self.checkpoint_path}/A_pred_{self.type}.h5ad')

    def recon_loss(self, z, data, norm, weight):
        r"""VGAE reconstruction loss.
        """
        y_pred = self.model_vgae.decoder(z, data.edge_label_index, sigmoid=True)
        loss = (
            norm * binary_cross_entropy(y_pred, data.edge_label, weight=weight).mean()
        )
        return loss

    def init_vgae(self, h, use_kd):
        r"""D-RP student model training function.

        Parameters:
        h : torch.Tensor
            Potential representation of the D-RP teacher model.
        use_kd : bool
            Utilize knowledge distillation.
        """
        seed_torch()
        self.h = h
        self.edge_index = self.get_edge_index(self.RP_Matrix.TR, self.RP_Matrix.Sample)
        print("Edge sum:", self.edge_index.shape[1])
        self.data = Data(x=self.data, edge_index=self.edge_index)
        transform = T.Compose([
            T.NormalizeFeatures(),
            T.ToDevice(self.device),
            T.RandomLinkSplit(num_val=0, num_test=0, is_undirected=True,
                            neg_sampling_ratio=0, add_negative_train_samples=False),
        ])
        train_data, _, _ = transform(self.data)
        x = train_data.x
        edge_index = train_data.edge_label_index
        self.num_features = self.data.num_features
        def get_params(data):
            num_nodes = data.num_nodes
            y = data.edge_label
            pos_weight = (num_nodes**2 - y.sum()) / y.sum()
            weight = torch.ones_like(y,device=self.device)
            weight[y.to(bool)] = pos_weight
            norm = num_nodes**2 / ((num_nodes**2 - y.sum()) * 2)
            return num_nodes, norm, weight
        
        num_nodes, norm, weight = get_params(train_data)

        encoder = VariationalGCNEncoder(self.data.num_features, self.h_dim, self.z_dim)
        decoder = InnerProductDecoderWeight(self.A_e)
        self.model_vgae = VGAE(encoder, decoder).to(self.device)
        self.optimizer_vgae = torch.optim.Adam(self.model_vgae.parameters(), lr=0.01)

        if os.path.exists(f'{self.checkpoint_path}/model_student_{self.type}.pt'):
            checkpoint = torch.load(
                f'{self.checkpoint_path}/model_student_{self.type}.pt',
                map_location=self.device,
            )
            self.model_vgae.load_state_dict(checkpoint['model_state_dict'])
            self.optimizer_vgae.load_state_dict(checkpoint['optimizer_state_dict'])
        else:
            print('# Training ...')
            for epoch in range(self.epochs_vgae):
                self.model_vgae.train()
                self.optimizer_vgae.zero_grad()
                h = self.model_vgae.encoder.predict_h(x, edge_index)
                z = self.model_vgae.encode(x, edge_index)
                if use_kd:
                    kd_loss = mse_loss(h, self.h)
                else:
                    kd_loss = 0
                re_loss = self.recon_loss(z, train_data, norm, weight)
                kl_loss = (1 / num_nodes) * self.model_vgae.kl_loss()
                loss = re_loss + kl_loss + kd_loss
                loss.backward()
                self.optimizer_vgae.step()
                if epoch % 1 == 0:
                    print(
                        f'epoch: {epoch}, loss: {loss}, kl_loss: {kl_loss}, re_loss: {re_loss}, kd_loss: {kd_loss}'
                    )
                    checkpoint = {
                        'model_state_dict': self.model_vgae.state_dict(),
                        'optimizer_state_dict': self.optimizer_vgae.state_dict(),
                    }
                    torch.save(
                        checkpoint,
                        f'{self.checkpoint_path}/model_student_{self.type}.pt',
                    )
            torch.save(
                checkpoint,
                f'{self.checkpoint_path}/model_student_{self.type}.pt',
            )
        self.save_graph()

    def init_cvae(self):
        r"""D-RP teacher model training function.
        """
        seed_torch()
        condition = np.ones(self.data.shape[0])
        condition[: self.RP_Matrix.TR.shape[0]] = 0
        condition = np.eye(2)[condition.astype(int)]
        condition = normalize(condition, 'l1', axis=0)
        condition = torch.tensor(condition, dtype=torch.float32).to(self.device)
        dataset = torch.concat([self.data, condition], dim=1)
        batch_size = 32
        train_loader = DataLoader(
            dataset=dataset,
            batch_size=batch_size,
            shuffle=True,
        )
        self.model_cvae = CVAE(
            self.data.shape[1], condition.shape[1], self.h_dim, self.z_dim
        ).to(self.device)
        self.optimizer_cvae = torch.optim.Adam(self.model_cvae.parameters(), lr=0.01)
        loss_func = nn.MSELoss()
        if os.path.exists(f'{self.checkpoint_path}/model_teacher_{self.type}.pt'):
            checkpoint = torch.load(
                f'{self.checkpoint_path}/model_teacher_{self.type}.pt',
                map_location=self.device,
            )
            self.model_cvae.load_state_dict(checkpoint['model_state_dict'])
            self.optimizer_cvae.load_state_dict(checkpoint['optimizer_state_dict'])
        else:
            print('# Training ...')
            for epoch in range(self.epochs_cvae):
                for batch in tqdm(train_loader):
                    self.model_cvae.train()
                    self.optimizer_cvae.zero_grad()
                    decoded = self.model_cvae(batch)
                    kl_loss = self.model_cvae.kl_div() / batch_size
                    re_loss = (
                        loss_func(decoded.view(-1), batch[:, :-2].reshape(-1))
                        / batch_size
                    )
                    loss = re_loss + kl_loss
                    loss.backward()
                    self.optimizer_cvae.step()
                print(
                    f'epoch: {epoch}, loss: {loss}, kl_loss: {kl_loss}, re_loss: {re_loss}'
                )
                checkpoint = {
                    'model_state_dict': self.model_cvae.state_dict(),
                    'optimizer_state_dict': self.optimizer_cvae.state_dict(),
                }
                torch.save(
                    checkpoint,
                    f'{self.checkpoint_path}/model_teacher_{self.type}.pt',
                )

        return self.model_cvae.predict_h(dataset).detach()

    def run(self,use_kd=True):
        r"""D-RP model network reconstruction module execution entry point.
        """
        self.epochs_cvae = 100
        self.epochs_vgae = 1000
        self.h_dim = 48
        self.z_dim = 24
        # Teacher
        z = self.init_cvae()
        # Student
        self.init_vgae(z,use_kd)

if __name__ == '__main__':
    type = sys.argv[1]
    library = sys.argv[2]
    use_kd = True
    checkpoint_path = library
    if len(sys.argv) > 3:
        use_kd = False if sys.argv[3].lower() == "false" else True
    if len(sys.argv) > 4:
        checkpoint_path = sys.argv[4]

    assert type in ['H3K27ac', 'ATAC']

    class RP_Matrix:
        TR = RPMatrix(library, 'RP_Matrix_TR.h5ad').norm('l1').get_data()
        Sample = RPMatrix(library, f'RP_Matrix_{type}.h5ad').norm('l1').get_data()

    CalcSTM(RP_Matrix, type, checkpoint_path).run(use_kd)
