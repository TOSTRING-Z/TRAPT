import os
import numpy as np
import pandas as pd
import argparse
import re

from TRAPT.Tools import RP_Matrix, Args, Type
from TRAPT.CalcTRAUC import CalcTRAUC
from TRAPT.DLFS import FeatureSelection

def runTRAPT(args: Args):
    r"""TRAPT execution entry function.

    Attributes
    ----------
    args : TRAPT.Tools.Args
        Global parameters for TRAPT.
    
    Returns:
        A pd.DataFrame of TR activity.
    """
    rp_matrix = RP_Matrix(args.library)
    obs = rp_matrix.TR.obs

    if os.path.exists(f"{args.output}/H3K27ac_RP.csv"):
        H3K27ac_RP = pd.read_csv(f"{args.output}/H3K27ac_RP.csv", header=None)[0]
    else:
        FS_H3K27ac = FeatureSelection(args, rp_matrix.H3K27ac, Type.H3K27ac)
        H3K27ac_RP,H3K27ac_info_samples = FS_H3K27ac.run()
        H3K27ac_RP.to_csv(f"{args.output}/H3K27ac_RP.csv", index=False, header=False)
        H3K27ac_info_samples.to_csv(f"{args.output}/H3K27ac_info_samples.csv", header=False)

    if os.path.exists(f"{args.output}/ATAC_RP.csv"):
        ATAC_RP = pd.read_csv(f"{args.output}/ATAC_RP.csv", header=None)[0]
    else:
        FS_ATAC = FeatureSelection(args, rp_matrix.ATAC, Type.ATAC)
        ATAC_RP,ATAC_info_samples = FS_ATAC.run()
        ATAC_RP.to_csv(f"{args.output}/ATAC_RP.csv", index=False, header=False)
        ATAC_info_samples.to_csv(f"{args.output}/ATAC_info_samples.csv", header=False)

    if os.path.exists(f"{args.output}/RP_TR_H3K27ac_auc.csv"):
        RP_TR_H3K27ac_auc = pd.read_csv(
            f"{args.output}/RP_TR_H3K27ac_auc.csv", index_col=0, header=None
        )
    else:
        H3K27ac_RP = H3K27ac_RP.values.flatten()
        CTR_TR = CalcTRAUC(args, rp_matrix.TR_H3K27ac, H3K27ac_RP)
        RP_TR_H3K27ac_auc = CTR_TR.run()
        RP_TR_H3K27ac_auc.to_csv(f"{args.output}/RP_TR_H3K27ac_auc.csv", header=False)

    if os.path.exists(f"{args.output}/RP_TR_ATAC_auc.csv"):
        RP_TR_ATAC_auc = pd.read_csv(
            f"{args.output}/RP_TR_ATAC_auc.csv", index_col=0, header=None
        )
    else:
        ATAC_RP = ATAC_RP.values.flatten()
        CTR_TR = CalcTRAUC(args, rp_matrix.TR_ATAC, ATAC_RP)
        RP_TR_ATAC_auc = CTR_TR.run()
        RP_TR_ATAC_auc.to_csv(f"{args.output}/RP_TR_ATAC_auc.csv", header=False)

    print("AUC integration...")
    data_auc = pd.concat([RP_TR_H3K27ac_auc, RP_TR_ATAC_auc], axis=1)
    data_auc /= np.linalg.norm(data_auc, axis=0, keepdims=True)
    TR_activity = pd.DataFrame(
        np.sum(data_auc.values, axis=1), index=data_auc.index, columns=[1]
    )
    TR_detail = pd.concat([TR_activity, data_auc], axis=1).reset_index()
    TR_detail.columns = ["TR", "TR activity", "RP_TR_H3K27ac_auc", "RP_TR_ATAC_auc"]
    obs.index.name = "TR"
    TR_detail = TR_detail.merge(obs.reset_index(), on="TR").sort_values(
        "TR activity", ascending=False
    )
    if args.tr_type != "all":
        TR_detail = TR_detail[TR_detail.TR.str.contains({
            "tf":"Sample",
            "tcof":"TcoF",
            "cr":"CR"
        }[args.tr_type], flags=re.IGNORECASE)]
    if args.source != "all":
        TR_detail = TR_detail.iloc[
            np.where(TR_detail.Source.str.contains(args.source, flags=re.IGNORECASE))
        ]
    TR_detail.to_csv(os.path.join(args.output, "TR_detail.txt"), index=False, sep="\t")
    print("Program execution completed.")
    return TR_detail

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    r"""TRAPT method entry function.
    """
    description = "usage: %prog [options] -l [LIBRARY] -i [INPUT] -o [OUTPUT]"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--library", "-l", type=str, default="library", help = "Enter the library path, default is './library'")
    parser.add_argument("--threads", "-p", type=int, default=16, help = "Number of threads to launch, default is 16")
    parser.add_argument("--trunk_size", "-s", type=int, default=32768, help = "Block size. If the memory is insufficient, set a smaller value. The default value is 32768")
    parser.add_argument("--input", "-i", type=str, default=None, help = "Enter a gene list")
    parser.add_argument("--output", "-o", type=str, default=None, help = "Enter an output folder")
    parser.add_argument("--background_genes", "-b", type=int, default=6000, help = "Background gene count")
    parser.add_argument("--use_kd", "-d", type=str2bool, default=True, help = "Using knowledge distillation")
    parser.add_argument("--tr_type", type=str, default="all", help = "all/tr/tcof/cr")
    parser.add_argument("--source", type=str, default="all", help = "all/cistrome/chip_altas/gtrd/remap/chip-atlas/remap/encode/geo")
    args = Args(**dict(parser.parse_args()._get_kwargs()))
    runTRAPT(args)

if __name__ == "__main__":
    main()
