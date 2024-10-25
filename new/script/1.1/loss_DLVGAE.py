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
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score, roc_curve, average_precision_score

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
        self.__mu__ = self.conv_mu(h, edge_index)
        self.__logstd__ = self.conv_logstd(h, edge_index)
        return self.__mu__, self.__logstd__


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

    def recon_loss(self, z, data, norm, weight):
        y_pred = self.model_vgae.decoder(z, data.edge_label_index, sigmoid=True)
        loss = (
            norm * binary_cross_entropy(y_pred, data.edge_label, weight=weight).mean()
        )
        return loss
    
    @staticmethod
    def test(model, data):
        model.eval()
        with torch.no_grad():
            z = model.encode(data.x, data.edge_index)
            out = model.decoder(z, data.edge_label_index)
            model.train()
        auroc = roc_auc_score(data.edge_label.cpu().numpy().tolist()+[0.], out.cpu().numpy().tolist()+[0.])
        auprc = average_precision_score(data.edge_label.cpu().numpy().tolist()+[0.], out.cpu().numpy().tolist()+[0.])
        return auroc,auprc

    def init_vgae(self, h, num_test, num_val, dl_weight = 1):
        seed_torch()
        self.h = h
        self.edge_index = self.get_edge_index(self.RP_Matrix.TR, self.RP_Matrix.Sample)
        print("Edge sum:", self.edge_index.shape[1])
        self.data = Data(x=self.data, edge_index=self.edge_index)
        transform = T.Compose([
            T.NormalizeFeatures(),
            T.ToDevice(self.device),
            T.RandomLinkSplit(num_val=num_val, num_test=num_test, is_undirected=True,
                            neg_sampling_ratio=1, add_negative_train_samples=False),
        ])
        train_data, val_data, _ = transform(self.data)
        x = train_data.x
        edge_index = train_data.edge_label_index
        x_val = val_data.x
        edge_index_val = val_data.edge_label_index
        self.num_features = self.data.num_features
        def get_params(data):
            num_nodes = data.num_nodes
            y = data.edge_label
            pos_weight = (num_nodes**2 - y.sum()) / y.sum()
            weight = torch.ones_like(y)
            weight[y.numpy().astype(bool)] = pos_weight
            weight = weight.to(self.device)
            norm = num_nodes**2 / ((num_nodes**2 - y.sum()) * 2)
            return num_nodes, norm, weight
        
        num_nodes, norm, weight = get_params(train_data)

        encoder = VariationalGCNEncoder(self.data.num_features, self.h_dim, self.z_dim)
        decoder = InnerProductDecoderWeight(self.A_e)
        self.model_vgae = VGAE(encoder, decoder).to(self.device)
        self.optimizer_vgae = torch.optim.Adam(self.model_vgae.parameters(), lr=0.01)

        print('# Training ...')
        metric_val = []
        metric_train = []
        for epoch in range(self.epochs_vgae):
            self.model_vgae.train()
            self.optimizer_vgae.zero_grad()
            h = self.model_vgae.encoder.predict_h(x, edge_index)
            z = self.model_vgae.encode(x, edge_index)
            dl_loss = mse_loss(h, self.h)
            re_loss = self.recon_loss(z, train_data, norm, weight)
            kl_loss = (1 / num_nodes) * self.model_vgae.kl_loss()
            loss = re_loss + kl_loss + dl_weight * dl_loss
            train_auroc,train_auprc = self.test(self.model_vgae, train_data)
            metric_train.append([loss.detach(),train_auroc,train_auprc])
            loss.backward()
            self.optimizer_vgae.step()
            if epoch % 1 == 0:
                self.model_vgae.eval()
                with torch.no_grad():
                    val_num_nodes, val_norm, val_weight = get_params(val_data)
                    h = self.model_vgae.encoder.predict_h(x_val, edge_index_val)
                    z = self.model_vgae.encode(x_val, edge_index_val)
                    dl_loss_val = mse_loss(h, self.h)
                    re_loss_val = self.recon_loss(z, val_data, val_norm, val_weight)
                    kl_loss_val = (1 / val_num_nodes) * self.model_vgae.kl_loss()
                    loss_val = re_loss_val + kl_loss_val + dl_weight * dl_loss_val
                    val_auroc,val_auprc = self.test(self.model_vgae, val_data)
                    metric_val.append([loss_val.detach(),val_auroc,val_auprc])
                print(
                    f"epoch: {epoch}, loss_val: {loss_val}, dl_loss_val: {dl_loss_val}, re_loss_val: {re_loss_val}, kl_loss_val: {kl_loss_val}"
                )
                print(
                    f"epoch: {epoch}, loss: {loss}, dl_loss: {dl_loss}, re_loss: {re_loss}, kl_loss: {kl_loss}"
                )
        return metric_val,metric_train

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

    def run(self,num_test=0.,num_val=0.5,dl_weight=0.1,epochs_cvae=100,epochs_vgae=10):
        self.epochs_cvae = epochs_cvae
        self.epochs_vgae = epochs_vgae
        self.h_dim = 48
        self.z_dim = 24
        # Teacher
        z,teacher_losses = self.init_cvae(num_val=0.1)
        # Student
        student_metric_val,student_metric_train = self.init_vgae(z,num_test=num_test, num_val=num_val, dl_weight=dl_weight)
        return np.asanyarray(student_metric_val),np.asanyarray(student_metric_train),np.asanyarray(teacher_losses)

def plot_lines(losses,label,color,save_path=None,title='Validation set loss',xlabel="Epochs",ylabel="Loss"):
    plt.plot(np.arange(len(losses)), losses, color, label=label, linewidth=2)
    plt.title(title, fontsize=14)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.legend(fontsize=10)
    plt.grid(True)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
        plt.clf()
        plt.close()

if __name__ == '__main__':
    type = "H3K27ac"
    checkpoint_path = "new/result/1.1/checkpoint_path_kd"
    library = "library"
    assert type in ['H3K27ac', 'ATAC']
    random.seed(33)
    np.random.seed(33)
    TR_index = random.choices(np.arange(17227),k=200)
    Sample_index = random.choices(np.arange(1465),k=200)
    class RP_Matrix_:
        TR = RPMatrix(library, 'RP_Matrix_TR.h5ad').norm('l1').get_data()[TR_index,:]
        Sample = RPMatrix(library, f'RP_Matrix_{type}.h5ad').norm('l1').get_data()[Sample_index,:]

    color_palette = ['#3c5488', '#f39b7f', '#8491b4', '#91d1c2', '#fccde5']

    student_metric_val,student_metric_train,teacher_losses = CalcSTM(RP_Matrix_, type, checkpoint_path).run(0.,0.1,1,epochs_cvae=10,epochs_vgae=500)
    plot_lines(teacher_losses,"D-RP teacher model",color_palette[0],"new/result/1.1/figure/loss-D-RP-teacher.svg", title="Validation set loss of the D-RP teacher model")

    student_metric_val_not_kd,student_metric_train_not_kd,_ =  CalcSTM(RP_Matrix_, type, "new/result/1.1/checkpoint_path_not_kd").run(0.0,0.1,0,epochs_cvae=0,epochs_vgae=500)

    plot_lines(student_metric_val[:,0],"Validation set (KD)",color_palette[1])
    plot_lines(student_metric_train[:,0],"Train set (KD)",color_palette[2])
    plot_lines(student_metric_val_not_kd[:,0],"Validation set (KD)",color_palette[3])
    plot_lines(student_metric_train_not_kd[:,0],"Train set (NKD)",color_palette[4], "new/result/1.1/figure/loss-D-RP-student.svg", title="Loss of the D-RP student model")

    plot_lines(student_metric_val[:,1],"D-RP student model (KD)",color_palette[2])
    plot_lines(student_metric_val_not_kd[:,1],"D-RP student model (NKD)",color_palette[3],"new/result/1.1/figure/auroc-D-RP-student.svg",title="Validation set AUROC for the D-RP model",ylabel="AURPC")

    plot_lines(student_metric_val[:,2],"D-RP student model (KD)",color_palette[2])
    plot_lines(student_metric_val_not_kd[:,2],"D-RP student model (NKD)",color_palette[3],"new/result/1.1/figure/auprc-D-RP-student.svg",title="Validation set AUPRC for the D-RP model",ylabel="AUPRC")