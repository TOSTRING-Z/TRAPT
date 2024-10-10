import os
import argparse
import pandas as pd
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--ssh", type=bool, default=True)
parser.add_argument("--ssh_user", type=str, default="zhangguorui")
parser.add_argument("--ssh_ip", type=str, default="172.27.0.6")
parser.add_argument("--ssh_port", type=int, default="22")
parser.add_argument("--marks_path", type=str, default="/wyzdata8/SEdb1/all_bam/sample")
parser.add_argument("--input", type=str, default="new/result/2.16/output/AD_200/H3K27ac_info_samples.csv")
parser.add_argument("--output_path", type=str, default="/data/zgr/data/TRAPT/tool/new/result/2.16/marks")
args = parser.parse_args()

def rsync(file,args):
    target = os.path.join(args.marks_path,file)
    if args.ssh:
        os.system(f"rsync -i --progress --copy-links -avzhe 'ssh -p {args.ssh_port}' {args.ssh_user}@{args.ssh_ip}:{target} {args.output_path}")
    else:
        os.system(f"rsync -i --progress --copy-links -avzhe {target} {args.output_path}")

info = pd.read_csv(args.input,header=None)[0]
for sample in tqdm(info):
    name = sample.split("_hg38")[0]
    file = f"{name}.sort.*"
    print(file)
    rsync(file,args)
