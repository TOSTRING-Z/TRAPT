import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
import argparse

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", type=str, default=None)
parser.add_argument("--output_dir", type=str, default=None)
parser.add_argument("--library", type=str, default='library')
parser.add_argument("--processes", type=int, default=2)
parser.add_argument("--restart_step", type=int, default=24)
parser.add_argument("--threads", type=int, default=16)
parser.add_argument("--background_genes", type=int, default=6000)
parser.add_argument("--use_kd", type=str2bool, default=True)
args = parser.parse_args()

import sys
import multiprocessing
from glob import glob

from tqdm import tqdm

from TRAPT.Run import runTRAPT
from TRAPT.Tools import Args

def execute_from_command_line(argv):
    python = sys.executable
    os.execl(python, python, *argv)


def get_params(input_dir, output_dir, library, restart_step, threads, trunk_size, background_genes, use_kd):
    step = 0
    for input in tqdm(sorted(glob(f'{input_dir}/*'), reverse=False)):
        name = os.path.basename(input).split('.')[0]
        final_output = os.path.join(output_dir, name)
        os.system(f'mkdir -p {final_output}')
        args = Args(input, final_output, library, threads, trunk_size, background_genes, use_kd)
        if os.path.exists(f'{final_output}/TR_detail.txt'):
            continue
        print(f'\033[1;31m # # # {input} # # #\033[0m')
        if step >= restart_step:
            return args
        else:
            step += 1
            yield args


if __name__ == '__main__':
    trunk_size = 2048 * args.threads

    step = 0
    with multiprocessing.Pool(args.processes) as pool:
        for data in pool.imap(
            runTRAPT,
            get_params(
                args.input_dir, 
                args.output_dir, 
                args.library, 
                args.restart_step, 
                args.threads, 
                trunk_size, 
                args.background_genes, 
                args.use_kd
            ),
        ):
            print("setp: %d ***************************" % step)
            step += 1

    if step:
        argv = sys.argv
        print("\033[1;31m restart...\033[0m", argv)
        execute_from_command_line(argv)
