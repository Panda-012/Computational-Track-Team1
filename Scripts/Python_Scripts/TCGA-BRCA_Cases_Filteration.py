import json

def meets_criteria(diagnosis):
    stage_criteria = {"Stage 0", "Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB"}
    t_criteria = {"T0", "T1", "T2", "T1a", "T1b", "T1c"}  # Expanded to include more specific T1 categories
    n_criteria = {"N0", "N1", "N0 (i+)", "N0 (mol+)", "N1mi"}  # Expanded to include more specific N0/N1 categories
    m_criteria = {"M0"}
    
    stage = diagnosis.get('ajcc_pathologic_stage', '')
    t = diagnosis.get('ajcc_pathologic_t', '')
    n = diagnosis.get('ajcc_pathologic_n', '')
    m = diagnosis.get('ajcc_pathologic_m', '')
    
    return (stage in stage_criteria) and (t in t_criteria) and (n in n_criteria) and (m in m_criteria)

def filter_cases(file_path):
    with open(file_path, 'r') as file:
        clinical_data = json.load(file)
        
    filtered_case_ids = [case['case_id'] for case in clinical_data if any(meets_criteria(d) for d in case.get('diagnoses', []))]
    return filtered_case_ids

# Replace with the path to your JSON file
file_path = 'TCGA-BRCA_Clinical_Data.json'
filtered_case_ids = filter_cases(file_path)

# Save filtered_case_ids to a .txt file
file_output_path = 'TCGA-BRCA_Filtered_Case_IDs.txt'
with open(file_output_path, 'w') as file_out:
    for case_id in filtered_case_ids:
        file_out.write(case_id + '\n')

print(f"Filtered Case IDs have been saved to {file_output_path}")
