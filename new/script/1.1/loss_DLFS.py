import random
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf
from keras import Input, Model
from keras import backend as K
from keras.callbacks import EarlyStopping
from keras.layers import Dense, Dropout
from keras.models import Model
from keras.optimizers import Adam
from keras.regularizers import Regularizer
import matplotlib.pyplot as plt
from sklearn.metrics import auc, precision_recall_curve, roc_curve
from TRAPT.Tools import RP_Matrix, Args, Type
from TRAPT.Run import str2bool
from TRAPT.DLFS import *

class FeatureSelectionLoss(FeatureSelection):
    def __init__(self, args, data_ad, type):
        super().__init__(args, data_ad, type)
        self.epochs = 100

    def TSFS(self, X, T):
        seed_tensorflow()
        pos_weight = float(len(T) - T.sum()) / T.sum()
        sample_weight = np.ones(len(T))
        sample_weight[T.astype(bool)] = pos_weight
        self.history = {}
        if self.args.use_kd:
            input = Input(X.shape[1:])
            y = Dense(
                2, 
                activation="relu", 
                kernel_regularizer=SparseGroupLasso(groups=self.groups)
            )(input)
            feature_extraction = Model(input, y)
            t = feature_extraction(input)
            output = Dense(1)(t)
            output = CustomSigmoid(t=1)(output)
            self.teacher = Model(input, output)
            self.teacher.compile(optimizer=Adam(self.learning_rate), loss="mse", weighted_metrics=[])
            self.history["teacher"] = self.teacher.fit(
                X,
                T,
                epochs=self.epochs,
                batch_size=self.batch_size,
                callbacks=[self.early_stopping],
                validation_split=0.1,
                sample_weight=sample_weight,
                shuffle=self.shuffle,
                verbose=1,
            )
            Y = feature_extraction.predict(X)
            seed_tensorflow()
            input = Input(X.shape[1:])
            output = Dense(
                Y.shape[-1] * 10,
                activation="relu",
                kernel_regularizer=SparseGroupLasso(groups=self.groups),
            )(input)
        else:
            Y = np.expand_dims(T,-1)
            input = Input((X.shape[-1],))
            output = Dense(
                Y.shape[-1] * 10 * 2,
                activation="relu",
                kernel_regularizer=SparseGroupLasso(groups=self.groups),
            )(input)
        o_1 = Dense(Y.shape[-1], activation="relu")(output)
        self.student = Model(input, o_1)
        self.student.compile(optimizer=Adam(self.learning_rate), loss="mse", weighted_metrics=[])
        self.history["student"] = self.student.fit(
            X,
            Y,
            epochs=self.epochs,
            batch_size=self.batch_size,
            validation_split=0.1,
            sample_weight=sample_weight,
            shuffle=self.shuffle,
            verbose=1,
        )
        weights = self.student.layers[1].get_weights()[0]
        weights = np.nan_to_num(weights, nan=0)
        weights = np.sum(np.square(weights), 1)
        return np.argsort(-weights),weights

def plot_lines(losses,label,color,save_path=None,title='Validation set loss',xlabel="Epochs",ylabel="Loss"):
    x = np.arange(len(losses))
    y = np.asanyarray(losses)
    plt.plot(x, y, color, label=label, linewidth=2)
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

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--library", "-l", type=str, default="library", help = "Enter the library path, default is './library'")
    parser.add_argument("--input", "-i", type=str, default="input/KnockTFv1/down/AR@DataSet_01_192_down500.txt", help = "Enter a gene list")
    parser.add_argument("--output", "-o", type=str, default="new/result/1.1/files", help = "Enter an output folder")
    parser.add_argument("--background_genes", "-b", type=int, default=6000, help = "Background gene count")
    parser.add_argument("--use_kd", "-d", type=str2bool, default=True, help = "Using knowledge distillation")
    args = Args(**dict(parser.parse_args()._get_kwargs()))
    rp_matrix = RP_Matrix(args.library)
    color_palette = ['#3c5488', '#f39b7f', '#8491b4', '#91d1c2', '#fccde5']

    args.use_kd = True
    FS_H3K27ac = FeatureSelectionLoss(args, rp_matrix.H3K27ac, Type.H3K27ac)
    H3K27ac_RP,H3K27ac_info_samples = FS_H3K27ac.run()
    teacher_loss = FS_H3K27ac.history["teacher"].history
    student_loss_kd = FS_H3K27ac.history["student"].history
    plot_lines(teacher_loss["loss"],"U-RP teacher model",color_palette[0],"new/result/1.1/figure/loss-U-RP-teacher.svg", title="Validation set loss of the U-RP teacher model")

    args.use_kd = False
    FS_H3K27ac = FeatureSelectionLoss(args, rp_matrix.H3K27ac, Type.H3K27ac)
    H3K27ac_RP,H3K27ac_info_samples = FS_H3K27ac.run()
    student_loss_not_kd = FS_H3K27ac.history["student"].history

    plot_lines(student_loss_kd["val_loss"],"Validation set (KD)",color_palette[1])
    plot_lines(student_loss_kd["loss"],"Train set (KD)",color_palette[2])
    plot_lines(student_loss_not_kd["val_loss"],"Validation set (NKD)",color_palette[3])
    plot_lines(student_loss_not_kd["loss"],"Train set (NKD)",color_palette[4], "new/result/1.1/figure/loss-U-RP-student.svg", title="Loss of the U-RP student model")