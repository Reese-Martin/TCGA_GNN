import requests
import json
import os
from tqdm import tqdm
import pandas as pd

# Directory to store downloaded files
download_dir = "tcga_data"
os.makedirs(download_dir, exist_ok=True)

# GDC API endpoint
files_endpt = "https://api.gdc.cancer.gov/files"

# Filter for TCGA-LUAD RNA-seq data
params = {
    "filters": json.dumps({
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LUAD"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}}
        ]
    }),
    "fields": "file_id,file_name",
    "format": "JSON",
    "size": "1000"
}


# Requesting the list of files
# response = requests.get(files_endpt, params=params)
# file_list = response.json()['data']['hits']
# file_ids = [f['file_id'] for f in file_list]

r = requests.get("https://api.gdc.cancer.gov/files", params=params)
file_list = r.json()["data"]["hits"]

# Create manifest-like dataframe
file_ids = [{"id": f["file_id"]} for f in file_list]
pd.DataFrame(file_ids).to_csv("gdc_manifest.txt", sep="\t", index=False)
# Now use the line below paired with gdc-client to download manifest data
# gdc-client download -m gdc_manifest.txt
# Note: if gdc-client has not been added to your path, you need to provide /your/path/to/gdc-client rather than just
# gdc-client

