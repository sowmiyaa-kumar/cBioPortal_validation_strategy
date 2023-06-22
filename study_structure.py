#!/usr/bin/env python
# coding: utf-8

# In[5]:


import os
import json
import regex as re

def get_data_filename(meta_file_path):
    with open(meta_file_path, 'r') as file:
        for line in file:
            if line.startswith('data_filename'):
                return line.split(': ')[1].strip()  # return the value after ': '
    return None  # return None if 'data_filename' is not found


def validate_directory(directory):
    # Check if directory exists
    if not os.path.exists(directory):
        raise Exception(f"Directory not found: {directory}")

    # Check for data and meta files
    meta_files = []
    data_files = []
    meta_study_file = None
    meta_clinical_sample_file = None
    for file in sorted(os.listdir(directory)):
        # Checks for meta file naming conventions
        # Check if file name contains the word "meta" as a whole word or as a part of another word, 
        # with an optional underscore or numeric character at the end.
        if re.search(r'(\b|_)meta(\b|[_0-9])', file, flags=re.IGNORECASE) and not file.startswith('ONCOKB_IMPORT_BACKUP') and not file.startswith('.') and not file.endswith('~'):
            meta_files.append(file)
            if "study" in file:
                meta_study_file = file
            if "clinical_sample" in file:
                meta_clinical_sample_file = file
        else:
            data_files.append(file)
    
    if len(meta_files) == 0:
        raise Exception(f'No meta files found in {directory}. Please make sure the directory '\
                        'is the path to the folder containing the files.')

    # Check for mandatory meta_study and meta_clinical_sample files
    if not meta_study_file or not meta_clinical_sample_file:
        raise Exception("Mandatory meta_study or meta_clinical_sample file not found")

    # Each data file should have a corresponding meta file (except for meta_study)
    for meta_file in meta_files:
        if meta_file != meta_study_file:
            data_filename = get_data_filename(os.path.join(directory, meta_file))
            if data_filename is None:
                raise Exception(f"Missing 'data_filename' in meta file: {meta_file}")
            if data_filename not in data_files:
                raise Exception(f"Missing data file for meta file: {meta_file}")
    
    # File-specific checks 
    # Check for mutation file and case list
    if 'mutation.txt' in data_files:
        # Case lists should be placed in a subdirectory called 'case_lists'
        if not os.path.exists(os.path.join(directory, 'case_list')):
            raise Exception(f"Sub-directory not found: 'case_lists'")
        else:
            if not os.path.exists(os.path.join(directory, 'case_lists', 'cases_sequenced.txt')):
                raise Exception("Mutation file found but missing _sequenced case list")

    return "Directory validation passed"


# Call the function with a directory path
print(validate_directory("cancer_studies/sample_cancer_study"))

