# using the gdc manifest that was created from Generate_Manifest.py, download the metadata for each patient file
import pandas as pd
import requests
# import json

manifest = pd.read_csv('star_counts/gdc_manifest.txt')
file_ids = manifest['id'].tolist()
fields = [
    "file_id",
    "file_name",
    "cases.submitter_id",
    "cases.project.project_id",
    "cases.diagnoses.age_at_diagnosis",
    "cases.diagnoses.vital_status",
    "cases.diagnoses.tumor_stage",
    "cases.diagnoses.tumor_grade",
    "cases.diagnoses.progression_or_recurrence",
    "cases.diagnoses.days_to_recurrence",
    "cases.diagnoses.days_to_death",
    "cases.demographic.gender",
    "cases.demographic.race",
    "cases.demographic.ethnicity",
    "cases.samples.sample_type",
    "cases.samples.submitter_id"
]

query = {
    "filters": {"op": "in", "content": {"field": "file_id", "value": file_ids}},
    "format": "JSON",
    "size": len(file_ids),
    "fields": ",".join(fields)
}

response = requests.post("https://api.gdc.cancer.gov/files", json=query)
response.raise_for_status()
data = response.json()["data"]["hits"]

metadata = []
for entry in data:
    case = entry.get("cases", [{}])[0]
    diagnosis = case.get("diagnoses", [{}])[0]
    demographic = case.get("demographic", {})
    sample = case.get("samples", [{}])[0]

    metadata.append({
        "file_id": entry.get("file_id"),
        "file_name": entry.get("file_name"),
        "project_id": case.get("project_id"),
        "case_barcode": case.get("submitter_id"),
        "sample_barcode": sample.get("submitter_id"),
        "sample_type": sample.get("sample_type"),
        "age_at_diagnosis": diagnosis.get("age_at_diagnosis"),
        "tumor_stage": diagnosis.get("tumor_stage"),
        "tumor_grade": diagnosis.get("tumor_grade"),
        "progression_or_recurrence": diagnosis.get("progression_or_recurrence"),
        "days_to_recurrence": diagnosis.get("days_to_recurrence"),
        "days_to_death": diagnosis.get("days_to_death"),
        "vital_status": diagnosis.get("vital_status"),
        "gender": demographic.get("gender"),
        "race": demographic.get("race"),
        "ethnicity": demographic.get("ethnicity")
    })

df = pd.DataFrame(metadata)
df.to_csv("Pt_metadata.csv", index=False)