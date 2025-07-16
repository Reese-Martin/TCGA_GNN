import requests
import json
import os
import pandas as pd

# Directory to store downloaded files
download_dir = "tcga_data"
os.makedirs(download_dir, exist_ok=True)

# GDC API endpoint
files_endpt = "https://api.gdc.cancer.gov/files"

# Filter for TCGA-LUAD RNA-seq data
filters = {"op": "and", "content": [
          {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LUAD"]}},
          {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
          {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
          {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Solid Tissue Normal"]}}]}

params = {
    "filters": json.dumps(filters),
    "fields": "file_id,file_name",
    "format": "JSON",
    "size": "1000"}


# Requesting the list of files

r = requests.get("https://api.gdc.cancer.gov/files", params=params)
file_list = r.json()["data"]["hits"]

# Create manifest-like dataframe
file_ids = [{"id": f["file_id"]} for f in file_list]
pd.DataFrame(file_ids).to_csv("gdc_manifest.txt", sep="\t", index=False)
# Now use the line below paired with gdc-client to download manifest data
# gdc-client download -m gdc_manifest.txt
# Note: if gdc-client has not been added to your path, you need to provide /your/path/to/gdc-client rather than just
# gdc-client

### Create manifest to download all normal samples available in the TCGA,
# I am doing this to help overcome the sample size issue, where cancer samples often outnumber normal solid tissue
# samples by a significant margin (in LUAD as of 06/25, that is 9:1)

filters = {"op": "and", "content": [
        {"op": "in", "content": {"field": "data_category", "value": "Transcriptome Profiling"}},
        {"op": "in", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
        {"op": "in", "content": {"field": "experimental_strategy", "value": "RNA-Seq"}},
        {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Solid Tissue Normal"]}}]}

params = {
    "filters": json.dumps(filters),        # Must be a string
    "fields": "file_id,file_name,submitter_id,cases.case_id,cases.submitter_id,sample_type",
    "format": "JSON",
    "size": 10000                          # Max results per request
}


# Requesting the list of files

r = requests.get("https://api.gdc.cancer.gov/files", params=params)
file_list = r.json()["data"]["hits"]

# Create manifest-like dataframe
file_ids = [{"id": f["file_id"]} for f in file_list]
pd.DataFrame(file_ids).to_csv("gdc_manifest_SolidTissueNormal.txt", sep="\t", index=False)

