import json
import os
from multiprocessing import Pool

import pandas as pd
import requests
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_path", type=str, default="input/Lisa/down")
parser.add_argument("--output_path", type=str, default="new/result/2.11/chea3/down")
parser.add_argument("--processes", type=int, default=10)
args = parser.parse_args()

def run(args):
    
    input,output,name = args
    print(name)
    gene_set = list(pd.read_csv(f"{input}/{name}.txt",header=None).dropna()[0])[:400]
    url = "https://maayanlab.cloud/chea3/api/enrich/"
    headers = {"Content-Type": "application/json"}
    params = {
        "query_name":"myQuery", 
        "gene_set": gene_set
    }
    try:
        response = requests.get(url=url,data="json",params=params,headers=headers)
        data = json.loads(response.text)
    except Exception as e:
        print("Error: ",e)
        return
    for key in data.keys():
        sub_data = pd.DataFrame(data.get(key))
        sub_data.to_csv(f"{output}/{name}${key}.tsv",sep="\t",index=False)

if __name__ == "__main__":
    os.makedirs(args.output_path,exist_ok=True)
    all = set(map(lambda x:x.split("$")[0],os.listdir(args.output_path)))
    params = [(args.input_path,args.output_path,x.split('.')[0]) for x in os.listdir(args.input_path) if x.split('.')[0] not in all]
    pool = Pool(args.processes)
    pool.map(run,params)
    

