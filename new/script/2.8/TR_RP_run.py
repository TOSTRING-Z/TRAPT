import os
from CalcTRAUC import CalcTRAUC
from Tools import Args,RP_Matrix
from tqdm import tqdm

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--library", type=str, default="library")
parser.add_argument("--input_path", type=str, default="input/KnockTFv1/down")
parser.add_argument("--output_path", type=str, default="new/result/2.8/output-TR_RP")
args_ = parser.parse_args()


files = os.listdir(args_.input_path)
for file in tqdm(files):
    input = os.path.join(args_.input_path,file)
    output = os.path.join(args_.output_path,file.split(".")[0])
    args = Args(input, output, args_.library)
    rp_matrix = RP_Matrix(args.library)
    CTR_TR = CalcTRAUC(args, rp_matrix.TR, 1)
    RP_TR_auc = CTR_TR.run()
    RP_TR_auc.columns = ["RP_TR_auc"]
    RP_TR_auc["tr_base"] = RP_TR_auc.index.str.extract("(.*?)@", expand=False)
    RP_TR_auc.to_csv(f"{args.output}/TR_detail.txt", sep="\t")