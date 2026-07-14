#!/usr/bin/env python
import os
import glob
import requests
import subprocess

article_id="32301477"
version="1"

article = requests.get(f"https://api.figshare.com/v2/articles/{article_id}/versions/{version}").json()
for file in article["files"]:
    print(f'Downloading {file["name"]} from {file["download_url"]}')
    subprocess.run(["wget", file["download_url"], "-O", f"data/{file['name']}"], check=True)

# (cd data && for i in *.zip; do unzip $i; done)
for zip_file in glob.glob(os.path.join("data", "*.zip")):
    subprocess.run(["unzip", os.path.basename(zip_file)], cwd="data", check=True)

# find data -type f -name "*.gz" -exec gunzip {} \;
for gz_file in glob.glob("data/**/*.gz", recursive=True):
    subprocess.run(["gunzip", gz_file], check=True)

# rm data/*.zip
for zip_file in glob.glob(os.path.join("data", "*.zip")):
    os.remove(zip_file)

print("Data downloaded and extracted")
