#!/usr/bin/env python
# coding: utf-8

# Use Cerberus for meta file and case list validation
import os
import regex as re
from collections import OrderedDict
from cerberus import Validator
from cerberus_schemas import META_SCHEMA_MAP
import logging

# Function to parse file to an ordered dictionary 
def parse_file_to_ordered_dict(file_path):
    meta_dict = OrderedDict()
    with open(file_path, 'r') as metafile:
        for line_index, line in enumerate(metafile):
            line = line.strip()
            if not line:
                continue
            if ':' not in line:
                raise Exception(f"\nInvalid file entry, no ':' found")

            key, value = map(str.strip, line.split(':', 1))
            meta_dict[key] = value
    return meta_dict

# Function to get the meta file type 
def get_meta_file_type(meta_dict):
    alt_type_datatype_to_meta = {
        # cancer type
        ("CANCER_TYPE", "CANCER_TYPE"): "CANCER_TYPE",
        # clinical and timeline
        ("CLINICAL", "PATIENT_ATTRIBUTES"): "PATIENT_ATTRIBUTES",
        ("CLINICAL_SAMPLE", "SAMPLE_ATTRIBUTES"): "SAMPLE_ATTRIBUTES",
        ("CLINICAL", "TIMELINE"): "TIMELINE",
        # rppa and mass spectrometry
        ("PROTEIN_LEVEL", "LOG2-VALUE"): "PROTEIN",
        ("PROTEIN_LEVEL", "Z-SCORE"): "PROTEIN",
        ("PROTEIN_LEVEL", "CONTINUOUS"): "PROTEIN",
        # cna
        ("COPY_NUMBER_ALTERATION", "DISCRETE"): "CNA_DISCRETE",
        ("COPY_NUMBER_ALTERATION", "DISCRETE_LONG"): "CNA_DISCRETE_LONG",
        ("COPY_NUMBER_ALTERATION", "CONTINUOUS"): "CNA_CONTINUOUS",
        ("COPY_NUMBER_ALTERATION", "LOG2-VALUE"): "CNA_LOG2",
        ("COPY_NUMBER_ALTERATION", "SEG"): "SEG",
        # expression
        ("MRNA_EXPRESSION", "CONTINUOUS"): "EXPRESSION",
        ("MRNA_EXPRESSION", "Z-SCORE"): "EXPRESSION",
        ("MRNA_EXPRESSION", "DISCRETE"): "EXPRESSION",
        # mutations
        ("MUTATION_EXTENDED", "MAF"): "MUTATION",
        ("MUTATION_UNCALLED", "MAF"): "MUTATION_UNCALLED",
        # others
        ("METHYLATION", "CONTINUOUS"): "METHYLATION",
        ("GENE_PANEL_MATRIX", "GENE_PANEL_MATRIX"): "GENE_PANEL_MATRIX",
        ("STRUCTURAL_VARIANT", "SV"): "STRUCTURAL_VARIANT",
        # cross-sample molecular statistics (for gene selection)
        ("GISTIC_GENES_AMP", "Q-VALUE"): "GISTIC_GENES",
        ("GISTIC_GENES_DEL", "Q-VALUE"): "GISTIC_GENES",
        ("MUTSIG", "Q-VALUE"): "MUTATION_SIGNIFICANCE",
        ("GENESET_SCORE", "GSVA-SCORE"): "GSVA_SCORES",
        ("GENESET_SCORE", "P-VALUE"): "GSVA_PVALUES",
        ("GENERIC_ASSAY", "LIMIT-VALUE"): "GENERIC_ASSAY_CONTINUOUS",
        ("GENERIC_ASSAY", "BINARY"): "GENERIC_ASSAY_BINARY",
        ("GENERIC_ASSAY", "CATEGORICAL"): "GENERIC_ASSAY_CATEGORICAL"
    }

    fields_to_meta = {
        # study
        ("cancer_study_identifier", "type_of_cancer"): "STUDY",
        # resource
        ("cancer_study_identifier", "resource_type"): {
            "PATIENT": "PATIENT_RESOURCES",
            "SAMPLE": "SAMPLE_RESOURCES",
            "STUDY": "STUDY_RESOURCES",
            "DEFINITION": "RESOURCES_DEFINITION"
        }
    }
    genetic_alteration_type = meta_dict.get('genetic_alteration_type', None)
    data_type = meta_dict.get('datatype', None)

    meta_file_type = alt_type_datatype_to_meta.get((genetic_alteration_type, data_type),
                    next((fields_to_meta[tuple_of_fields] for tuple_of_fields in fields_to_meta.keys() if set(tuple_of_fields).issubset(meta_dict.keys())), None))

    if meta_file_type is None:
        raise Exception("Could not determine the file type. Please check your meta files for correct configuration.")

    return meta_file_type

# Function to parse all the meta files 
def parse_metadata(study_dir, meta_files):
    meta = {} # Dictionary to store the parsed meta files based on their type
    for meta_file in meta_files:
        meta_path = os.path.join(study_dir, meta_file)
        meta_dict = parse_file_to_ordered_dict(meta_path)
        meta_file_type = get_meta_file_type(meta_dict)
        meta[meta_file_type] = meta_path
    return meta 

def validate_metadata(meta) -> None:
    for meta_file_type, meta_path in meta.items():
        logging.info(f'Starting validation of {meta_file_type}')
        meta_dict = parse_file_to_ordered_dict(meta_path)
        v = Validator()
        if v.validate(meta_dict, META_SCHEMA_MAP.get(meta_file_type)) != True:
            print("ERRORS:")
            print(v.errors)
            print()
        else:
            logging.info(f'Validation of {meta_file_type} complete without errors.\n')    

