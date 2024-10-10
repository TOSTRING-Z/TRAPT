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
from torch.utils.data import random_split
from torch.nn.functional import binary_cross_entropy, kl_div, mse_loss
from torch.utils.data import DataLoader
from torch_geometric.data import Data
from torch_geometric.nn import VGAE, GCNConv, InnerProductDecoder
from tqdm import tqdm

from TRAPT.Tools import RPMatrix

import matplotlib.pyplot as plt
from sklearn.metrics import auc, precision_recall_curve, roc_curve

def plot_auc_auprc(y_true, y_scores):
    fpr, tpr, _ = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    pr_auc = auc(recall, precision)
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC Curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend(loc="lower right")
    plt.savefig("new/result/1.1/figure/AUC-D-RP模型.pdf")
    plt.clf()
    plt.close()
    plt.plot(recall, precision, color='blue', lw=2, label=f'Precision-Recall Curve (AUC = {pr_auc:.2f})')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend(loc="lower left")
    plt.savefig("new/result/1.1/figure/AUPRC-D-RP模型.pdf")
    plt.clf()
    plt.close()


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
        AB = (AB_cos_sim >= np.sort(AB_cos_sim, axis=1)[:, -n].reshape(-1, 1)).astype(int)
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

    def recon_loss(self, z):
        y_pred = self.model_vgae.decoder(z, sigmoid=True)
        loss = (
            self.norm
            * binary_cross_entropy(y_pred.view(-1), self.y, weight=self.weight).mean()
        )
        return loss

    def init_vgae(self, h, num_test, num_val, dl_weight = 1):
        seed_torch()
        self.h = h
        self.edge_index = self.get_edge_index(self.RP_Matrix.TR, self.RP_Matrix.Sample)
        print("Edge sum:", self.edge_index.shape[1])
        self.data = Data(x=self.data, edge_index=self.edge_index)
        self.transform = T.RandomLinkSplit(
            num_val=num_val,
            num_test=num_test,
            neg_sampling_ratio=0,
            add_negative_train_samples=False,
        )
        self.train_data, self.val_data, self.test_data = self.transform(self.data)
        x = self.train_data.x.to(self.device)
        edge_index = self.train_data.edge_index.to(self.device)
        x_val = self.val_data.x.to(self.device)
        edge_index_val = self.val_data.edge_index.to(self.device)
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

        print('# Training ...')
        losses_val = []
        losses_train = []
        for epoch in range(self.epochs_vgae):
            self.model_vgae.train()
            self.optimizer_vgae.zero_grad()
            h = self.model_vgae.encoder.predict_h(x, edge_index)
            self.z = self.model_vgae.encode(x, edge_index)
            dl_loss = mse_loss(h, self.h)
            re_loss = self.recon_loss(self.z)
            kl_loss = (1 / self.num_nodes) * self.model_vgae.kl_loss()
            loss = re_loss + kl_loss + dl_weight * dl_loss
            losses_train.append(loss.detach())
            loss.backward()
            self.optimizer_vgae.step()
            if epoch % 1 == 0:
                self.model_vgae.eval()
                with torch.no_grad():
                    h = self.model_vgae.encoder.predict_h(x_val, edge_index_val)
                    z = self.model_vgae.encode(x_val, edge_index_val)
                    dl_loss_val = mse_loss(h, self.h)
                    re_loss_val = self.recon_loss(z)
                    kl_loss_val = (1 / self.num_nodes) * self.model_vgae.kl_loss()
                    loss_val = re_loss_val + kl_loss_val + dl_weight * dl_loss_val
                    losses_val.append(loss_val)
                print(
                    f"epoch: {epoch}, loss_val: {loss_val}, dl_loss_val: {dl_loss_val}, re_loss_val: {re_loss_val}, kl_loss_val: {kl_loss_val}"
                )
        return losses_val,losses_train

    def init_cvae(self,num_val=0.1):
        seed_torch()
        condition = np.ones(self.data.shape[0])
        condition[: self.RP_Matrix.TR.shape[0]] = 0
        condition = np.eye(2)[condition.astype(int)]
        condition = normalize(condition, 'l1', axis=0)
        condition = torch.tensor(condition, dtype=torch.float32).to(self.device)
        dataset = torch.concat([self.data, condition], dim=1)
        batch_size = 32
        train_dataset, val_dataset = random_split(dataset, [1 - num_val, num_val])
        train_loader = DataLoader(
            dataset=train_dataset,
            batch_size=batch_size,
            shuffle=True,
        )
        self.model_cvae = CVAE(
            self.data.shape[1], condition.shape[1], self.h_dim, self.z_dim
        ).to(self.device)
        self.optimizer_cvae = torch.optim.Adam(self.model_cvae.parameters(), lr=0.01)
        loss_func = nn.MSELoss()
        losses = []
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

            self.model_cvae.eval()
            with torch.no_grad():
                batch = val_dataset.dataset[val_dataset.indices]
                decoded = self.model_cvae(batch)
                kl_loss_val = self.model_cvae.kl_div() / len(val_dataset.indices)
                re_loss_val = (
                    loss_func(decoded.view(-1), batch[:, :-2].reshape(-1))
                    / len(val_dataset.indices)
                )
                loss_val = re_loss_val + kl_loss_val
                losses.append(loss_val)
                print(
                    f'epoch: {epoch}, loss_val: {loss_val}, kl_loss_val: {kl_loss_val}, re_loss_val: {re_loss_val}'
                )

        return self.model_cvae.predict_h(dataset).detach(),losses

    def run(self,num_test=0.2,num_val=0.1,dl_weight=0.1,epochs_cvae=100,epochs_vgae=10):
        self.epochs_cvae = epochs_cvae
        self.epochs_vgae = epochs_vgae
        self.h_dim = 48
        self.z_dim = 24
        # Teacher
        z,teacher_losses = self.init_cvae(num_val=num_val)
        # Student
        student_losses_val,student_losses_train = self.init_vgae(z,num_test=num_test, num_val=num_val, dl_weight=dl_weight)
        return student_losses_val,student_losses_train,teacher_losses

def plot_losses(losses,label,color,save_path,save=True):
    plt.plot(np.arange(len(losses)), losses, color, label=label, linewidth=2)
    plt.title('Validation Set Loss', fontsize=14)
    plt.xlabel('Epochs', fontsize=12)
    plt.ylabel('Loss', fontsize=12)
    plt.legend(fontsize=10)
    plt.grid(True)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    if save:
        plt.savefig(save_path)
        plt.clf()
        plt.close()

if __name__ == '__main__':
    type = "H3K27ac"
    checkpoint_path = "new/result/1.1/checkpoint_path_dl"
    library = "library"
    assert type in ['H3K27ac', 'ATAC']
    random.seed(2024)
    np.random.seed(2024)
    TR_index = random.choices(np.arange(17227),k=100)
    Sample_index = random.choices(np.arange(1465),k=100)
    class RP_Matrix_:
        TR = RPMatrix(library, 'RP_Matrix_TR.h5ad').norm('l1').get_data()[TR_index,:]
        Sample = RPMatrix(library, f'RP_Matrix_{type}.h5ad').norm('l1').get_data()[Sample_index,:]

    color_palette = ['#3c5488', '#f39b7f', '#8491b4', '#91d1c2', '#fccde5']

    student_losses_val,student_losses_train,teacher_losses = CalcSTM(RP_Matrix_, type, checkpoint_path).run(0.1,0.1,1,epochs_cvae=10,epochs_vgae=100)
    plot_losses(student_losses_val,"D-RP student model (val)",color_palette[2],"new/result/1.1/figure/loss-D-RP-student[DL].svg",save=False)
    plot_losses(student_losses_train,"D-RP student model (train)",color_palette[3],"new/result/1.1/figure/loss-D-RP-student[DL].svg")
    plot_losses(teacher_losses,"D-RP teacher model",color_palette[0],"new/result/1.1/figure/loss-D-RP-teacher.svg")

    student_losses_val_not_dl,student_losses_train_not_dl,_ =  CalcSTM(RP_Matrix_, type, "new/result/1.1/checkpoint_path_not_dl").run(0.1,0.1,0,epochs_cvae=0,epochs_vgae=100)
    plot_losses(student_losses_val_not_dl,"D-RP student model (val)",color_palette[2],"new/result/1.1/figure/loss-D-RP-student[not-DL].svg",save=False)
    plot_losses(student_losses_train_not_dl,"D-RP student model (train)",color_palette[3],"new/result/1.1/figure/loss-D-RP-student[not-DL].svg")
