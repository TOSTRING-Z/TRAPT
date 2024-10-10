import os
import re
import time

import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from tqdm import tqdm


def download_wait(directory, timeout, nfiles=None):
    """
    Wait for downloads to finish with a specified timeout.

    Args
    ----
    directory : str
        The path to the folder where the files will be downloaded.
    timeout : int
        How many seconds to wait until timing out.
    nfiles : int, defaults to None
        If provided, also wait for the expected number of files.

    """
    seconds = 0
    dl_wait = True
    while dl_wait and seconds < timeout:
        time.sleep(1)
        dl_wait = False
        files = os.listdir(directory)
        if nfiles and len(files) != nfiles:
            dl_wait = True

        for fname in files:
            if fname.endswith('.crdownload'):
                dl_wait = True

        seconds += 1
    return seconds

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_path", type=str, default="input/Lisa/down")
parser.add_argument("--output_path", type=str, default="new/result/2.11/icistarget/down")
parser.add_argument("--download_dir", type=str, default="/home/tostring/下载")
parser.add_argument("--executable_path", type=str, default=None)
args = parser.parse_args()

os.makedirs(args.output_path,exist_ok=True)
os.makedirs(args.download_dir,exist_ok=True)

# 设置Chrome选项，例如禁用图片加载以提高性能
chrome_options = Options()
chrome_options.add_argument("--disable-gpu")
chrome_options.add_argument("--no-sandbox")
chrome_options.add_argument("--disable-dev-shm-usage")

if args.executable_path:
    # 指定ChromeDriver的路径
    service = Service(executable_path=args.executable_path)
    # 创建一个Chrome浏览器实例
    web = webdriver.Chrome(service=service, options=chrome_options)
else:
    web = webdriver.Chrome(options=chrome_options)

print("String .....")
web.get('https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/index.php')
if __name__=="__main__":
    download_final = set(os.listdir(args.output_path))
    file_name_list = list(filter(lambda x:re.findall(r'(^.*?)\.txt$',x)[0] not in download_final,os.listdir(args.input_path)))
    # 自动点击开头的cookie
    web.implicitly_wait(5)
    try:
        web.find_element(By.XPATH,'//*[@id="onetrust-accept-btn-handler"]').click()
    except:
        pass
    for file in tqdm(file_name_list):
        try:
            data = pd.read_csv(f"{args.input_path}/{file}", header=None)
            # 第一个输入框
            input = web.find_element(By.XPATH,'//*[@id="query_data"]')
            input.send_keys('\n'.join(data[0]))
            # type
            click1 = web.find_element(By.XPATH,'//*[@id="query_type"]/option[1]')
            click1.click()
            click2 = web.find_element(By.XPATH,'//*[@id="databases"]/option[1]')
            click2.click()
            click7 = web.find_element(By.XPATH,'//*[@id="databases"]/option[3]')
            click7.click()
            click8 = web.find_element(By.XPATH,'//*[@id="databases"]/option[4]')
            click8.click()
            # 输入工作名
            work_name = re.findall(r'(^.*?)\.txt$',file)[0]
            input = web.find_element(By.XPATH,'//*[@id="job_name"]')
            input.send_keys(work_name)
            # 点击提交
            input = web.find_element(By.XPATH,'//*[@id="execute"]')
            input.click()
            WebDriverWait(web, 1000, 0.5).until(EC.presence_of_element_located((By.LINK_TEXT, "archive")))
            click3 = web.find_element(By.XPATH,'//*[@id="content"]/form/p[2]/a')
            os.system(f'rm -rf {args.download_dir}/archive*')
            click3.click()
            download_wait(args.download_dir, 1000)
            os.system(f'mv {args.download_dir}/archive.zip {args.download_dir}/{work_name}.zip')
            os.system(f'unzip -d {args.output_path}/{work_name} {args.download_dir}/{work_name}.zip')
            web.back()
            input = web.find_element(By.XPATH,'//*[@id="query_data"]')
            input.clear()
            input = web.find_element(By.XPATH,'//*[@id="job_name"]')
            input.clear()
        except:
            web.implicitly_wait(5)
            try:
                web.find_element(By.XPATH,'//*[@id="onetrust-accept-btn-handler"]').click()
            except:
                pass
            web.get('https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/index.php')