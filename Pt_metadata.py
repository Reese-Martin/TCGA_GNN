# using the gdc manifest that was created from Generate_Manifest.py, download the metadata for each patient file
import pandas as pd
import requests


def Collect_Metadata(manif, filename):
    manifest = pd.read_csv(manif)
    file_ids = manifest['id'].tolist()
    input_file_ids = manifest['id'].tolist()
    file_fields = "file_id,file_name,cases.case_id,cases.submitter_id,cases.samples.submitter_id"
    file_query = {"filters": {"op": "in", "content": {"field": "file_id", "value": input_file_ids}},
                  "format": "JSON",
                  "fields": file_fields,
                  "size": len(input_file_ids)}

    file_resp = requests.post("https://api.gdc.cancer.gov/files", json=file_query)
    file_data = file_resp.json()["data"]["hits"]

    # Build a lookup from file_id to (case_id, sample_id, etc.)
    file_lookup = {}
    case_ids = set()
    for entry in file_data:
        file_id = entry["file_id"]
        file_name = entry["file_name"]
        for case in entry.get("cases", []):
            case_id = case["case_id"]
            case_submitter_id = case.get("submitter_id")
            for sample in case.get("samples", []):
                sample_submitter_id = sample.get("submitter_id")
                file_lookup[file_id] = {"file_name": file_name,
                                        "case_id": case_id,
                                        "case_barcode": case_submitter_id,
                                        "sample_barcode": sample_submitter_id}
                case_ids.add(case_id)

    #  Query cases for metadata using case_ids
    case_fields = ["case_id",
                   "submitter_id",
                   "project.project_id",
                   "diagnoses.age_at_diagnosis",
                   "diagnoses.vital_status",
                   "diagnoses.tumor_stage",
                   "diagnoses.tumor_grade",
                   "diagnoses.progression_or_recurrence",
                   "diagnoses.days_to_recurrence",
                   "diagnoses.days_to_death",
                   "demographic.gender",
                   "demographic.race",
                   "demographic.ethnicity"]

    case_query = {"filters": {"op": "in", "content": {"field": "case_id", "value": list(case_ids)}},
                  "format": "JSON",
                  "fields": ",".join(case_fields),
                  "size": len(case_ids)}

    case_resp = requests.post("https://api.gdc.cancer.gov/cases", json=case_query)
    case_data = case_resp.json()["data"]["hits"]

    # Merge file_lookup with case metadata
    metadata = []
    for case in case_data:
        case_id = case["case_id"]
        project_id = case.get("project", {}).get("project_id")
        case_barcode = case.get("submitter_id")
        demographic = case.get("demographic", {})
        diagnoses = case.get("diagnoses", [])
        diagnosis = diagnoses[0] if diagnoses else {}

        # Find all file_ids that correspond to this case_id
        for file_id, file_info in file_lookup.items():
            if file_info["case_id"] == case_id:
                metadata.append({"file_id": file_id, "file_name": file_info["file_name"], "project_id": project_id,
                                 "case_id": case_id, "case_barcode": case_barcode,
                                 "sample_barcode": file_info["sample_barcode"],
                                 "age_at_diagnosis": diagnosis.get("age_at_diagnosis"),
                                 "tumor_stage": diagnosis.get("tumor_stage"),
                                 "tumor_grade": diagnosis.get("tumor_grade"),
                                 "progression_or_recurrence": diagnosis.get("progression_or_recurrence"),
                                 "days_to_recurrence": diagnosis.get("days_to_recurrence"),
                                 "days_to_death": diagnosis.get("days_to_death"),
                                 "vital_status": diagnosis.get("vital_status"),
                                 "gender": demographic.get("gender"),
                                 "race": demographic.get("race"),
                                 "ethnicity": demographic.get("ethnicity")})

        # convert to DataFrame
    metadata_df = pd.DataFrame(metadata)
    metadata_df.to_csv(filename)


# LUAD metadata collection
LUADmanif = 'LUAD_star_counts/gdc_manifest.txt'
LUADfilename = "LUAD_star_counts/Pt_metadata_test.csv"
Collect_Metadata(LUADmanif, LUADfilename)

# All Solid Tissue Normal metadata collection
STNmanif = 'STN_star_counts/gdc_manifest_SolidTissueNormal.txt'
STNfilename = "STN_star_counts/STN_metadata.txt"
Collect_Metadata(STNmanif, STNfilename)
