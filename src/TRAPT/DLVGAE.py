import os
import random
import sys

import anndata as ad
import numpy as np
import pandas as pd
import torch
import torch.nn.modules.loss
import torch_geometric.transforms as T
from scipy import sparse
from sklearn.preprocessing import normalize
from torch import Tensor, nn
from torch.nn.functional import binary_cross_entropy, kl_div, mse_loss
from torch.utils.data import DataLoader
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
    def __init__(self, RP_Matrix, type, checkpoint_path, device='cpu'):
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
        self.is_val = 0

    @staticmethod
    def get_cos_similar_matrix(m1, m2):
        num = m1.dot(m2.T)
        denom = np.linalg.norm(m1, axis=1).reshape(-1, 1) * np.linalg.norm(m2, axis=1)
        res = num / denom
        res[np.isneginf(res)] = 0
        return res

    def get_edge_index(self, A, B, n=10):
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

    def test(self, data):
        self.model_vgae.eval()
        with torch.no_grad():
            from sklearn.metrics import average_precision_score, roc_auc_score

            z = self.model_vgae.encode(
                data.x.to(self.device), data.edge_index.to(self.device)
            )
            y_pred = self.model_vgae.decoder(z, sigmoid=True)
            y_pred = y_pred.detach().cpu().numpy().flatten()
            index = np.random.choice(np.arange(len(self.y)), len(self.y) // 1000)
            return roc_auc_score(self.y[index], y_pred[index]), average_precision_score(
                self.y, y_pred
            )

    def save_graph(self):
        print(self.test(self.data))
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

    def recon_loss(self, z):
        y_pred = self.model_vgae.decoder(z, sigmoid=True)
        loss = (
            self.norm
            * binary_cross_entropy(y_pred.view(-1), self.y, weight=self.weight).mean()
        )
        return loss

    def init_vgae(self, h, use_dl):
        seed_torch()
        self.h = h
        self.edge_index = self.get_edge_index(self.RP_Matrix.TR, self.RP_Matrix.Sample)
        print("Edge sum:", self.edge_index.shape[1])
        self.data = Data(x=self.data, edge_index=self.edge_index)
        if self.is_val:
            self.transform = T.RandomLinkSplit(
                num_val=0.2,
                num_test=0,
                neg_sampling_ratio=0,
                add_negative_train_samples=False,
            )
            self.train_data, self.val_data, self.test_data = self.transform(self.data)
            x = self.train_data.x.to(self.device)
            edge_index = self.train_data.edge_index.to(self.device)
        else:
            x = self.data.x.to(self.device)
            edge_index = self.data.edge_index.to(self.device)
        self.num_features = self.data.num_features
        self.num_nodes = self.data.num_nodes
        values = torch.ones(self.edge_index.shape[1])
        size = torch.Size([self.num_nodes] * 2)
        y = torch.sparse_coo_tensor(self.edge_index, values, size)
        y = y.to_dense().view(-1)
        pos_weight = (self.num_nodes**2 - y.sum()) / y.sum()
        weight = torch.ones_like(y)
        weight[y.numpy().astype(bool)] = pos_weight
        self.weight = weight.to(self.device)
        self.norm = self.num_nodes**2 / ((self.num_nodes**2 - y.sum()) * 2)
        self.y = y.to(self.device)

        encoder = VariationalGCNEncoder(self.data.num_features, self.h_dim, self.z_dim)
        decoder = InnerProductDecoderWeight(self.A_e)
        self.model_vgae = VGAE(encoder, decoder).to(self.device)
        self.optimizer_vgae = torch.optim.Adam(self.model_vgae.parameters(), lr=0.01)
        if (
            os.path.exists(f'{self.checkpoint_path}/model_student_{self.type}.pt')
            and not self.is_val
        ):
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
                if use_dl:
                    dl_loss = mse_loss(h, self.h)
                else:
                    dl_loss = 0
                re_loss = self.recon_loss(z)
                kl_loss = (1 / self.num_nodes) * self.model_vgae.kl_loss()
                loss = re_loss + kl_loss + dl_loss
                loss.backward()
                self.optimizer_vgae.step()
                if epoch % 1 == 0:
                    if self.is_val:
                        auc, ap = self.test(self.val_data)
                        print(f"batch:", auc, ap)
                    print(
                        f"epoch: {epoch}, loss: {loss}, dl_loss: {dl_loss}, re_loss: {re_loss}, kl_loss: {kl_loss}"
                    )
                    if not self.is_val:
                        checkpoint = {
                            'model_state_dict': self.model_vgae.state_dict(),
                            'optimizer_state_dict': self.optimizer_vgae.state_dict(),
                        }
                        torch.save(
                            checkpoint,
                            f'{self.checkpoint_path}/model_student_{self.type}.pt',
                        )
            if not self.is_val:
                torch.save(
                    checkpoint,
                    f'{self.checkpoint_path}/model_student_{self.type}.pt',
                )
        if self.is_val:
            auc, ap = self.test(self.val_data)
            print(f"batch:", auc, ap)
        else:
            self.save_graph()

    def init_cvae(self):
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
        if (
            os.path.exists(f'{self.checkpoint_path}/model_teacher_{self.type}.pt')
            and not self.is_val
        ):
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
                if not self.is_val:
                    checkpoint = {
                        'model_state_dict': self.model_cvae.state_dict(),
                        'optimizer_state_dict': self.optimizer_cvae.state_dict(),
                    }
                    torch.save(
                        checkpoint,
                        f'{self.checkpoint_path}/model_teacher_{self.type}.pt',
                    )

        return self.model_cvae.predict_h(dataset).detach()

    def run(self,use_dl=True):
        self.epochs_cvae = 10
        self.epochs_vgae = 100
        # self.epochs_cvae = 1
        # self.epochs_vgae = 1
        self.h_dim = 48
        self.z_dim = 24
        # Teacher
        z = self.init_cvae()
        # Student
        self.init_vgae(z,use_dl)


if __name__ == '__main__':
    type = sys.argv[1]
    library = sys.argv[2]
    use_dl = True
    checkpoint_path = library
    if len(sys.argv) > 3:
        use_dl = False if sys.argv[3].lower() == "false" else True
    if len(sys.argv) > 4:
        checkpoint_path = sys.argv[4]

    assert type in ['H3K27ac', 'ATAC']

    class RP_Matrix:
        TR = RPMatrix(library, 'RP_Matrix_TR.h5ad').norm('l1').get_data()
        Sample = RPMatrix(library, f'RP_Matrix_{type}.h5ad').norm('l1').get_data()

    CalcSTM(RP_Matrix, type, checkpoint_path).run(use_dl)
