#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Use Cerberus for meta_file and case list validation
from collections import OrderedDict
import os

file_dir = 'cancer_studies/msk_impact_2017/'

# Parsing meta_study file into an OrderedDict
filename = os.path.join(file_dir, "meta/meta_study.txt")
meta_study = OrderedDict()
with open(filename, 'r') as metafile:
    # split each line into key and value
    for line in metafile:
        key, value = line.strip().split(': ')
        # add the key-value pair to the dictionary
        meta_study[key] = value
        
# Parsing meta_mutations file into an OrderedDict
filename = os.path.join(file_dir, "meta/meta_mutations.txt")
meta_mutations = OrderedDict()
with open(filename, 'r') as metafile:
    # split each line into key and value
    for line in metafile:
        key, value = line.strip().split(': ')
        # add the key-value pair to the dictionary
        meta_mutations[key] = value        


# Parsing cases_sequenced file into an OrderedDict
filename = os.path.join(file_dir, "case_lists/cases_sequenced.txt")
cases_sequenced = OrderedDict()
with open(filename, 'r') as file:
    # split each line into key and value
    for line in file:
        key, value = line.strip().split(': ')
        # add the key-value pair to the dictionary
        cases_sequenced[key] = value 


# In[ ]:


from cerberus import Validator 
import regex as re

""" Define a schema in the form of nested dictionary 
(more information on https://docs.python-cerberus.org/en/stable/schemas.html)

Can view all the available validation rules at 
https://docs.python-cerberus.org/en/stable/validation-rules.html."""

# Custom validation rules 
# defining a custom validation rule to check consistency of cancer_study_identifiers across meta files 
# can also do this by extending the Validator class
def check_cancer_study_identifier(field, value, error):
    if value != meta_study['cancer_study_identifier']:
        error(field, "The cancer study identifier does not match to that of meta_study.txt.")
        
def check_duplicate_sample_ids(field, value, error):
    sample_ids = [x.strip() for x in cases_sequenced['case_list_ids'].split('\t')]
    if len(sample_ids) != len(set(sample_ids)):
        error(field, "Duplicate sample ID in case list.")

# Schema for meta_mutations
meta_mutations_schema = {'cancer_study_identifier':
                            {'type': 'string',
                             'maxlength': 255, 
                             'check_with': check_cancer_study_identifier,
                             'required': True},
                        'genetic_alteration_type':
                            {'type': 'string',
                             'allowed': ['MUTATION_EXTENDED'], 
                             'required': True},
                        'datatype':
                             {'type': 'string',
                              'allowed': ['MAF'],
                              'required': True},
                        'stable_id':
                             {'type': 'string',
                              'allowed': ['mutations'],
                              'regex': r'^[A-Za-z0-9_-]+$',
                              'required': True},
                        'show_profile_in_analysis_tab':
                             {'type': 'string',
                              'allowed': ['true'],
                              'required': True},
                        'profile_name':
                             {'type': 'string',
                              'required': True},
                        'profile_description':
                             {'type': 'string',
                              'required': True},
                        'data_filename':
                             {'type': 'string',
                              'required': True},
                        'normal_samples_list':
                             {'type': 'string',
                              'required': False},
                        'swissprot_identifier':
                             {'type': 'string',
                              'allowed': ['name', 'accession'],
                              'required': False},
                        'gene_panel':
                             {'type': 'string',
                              'required': False},
                        'variant_classification_filter':
                             {'type': 'string',
                              'regex': r'^[^,\s]+(?:,[^,\s]+)*$', 
                              'required': False},
                        'namespaces':
                             {'type': 'string',
                              'regex': r'^[^,\s]+(?:,[^,\s]+)*$', 
                              'required': False},
                        }

# Schema for cases_sequenced 
cases_sequenced_schema = {'cancer_study_identifier':
                             {'type': 'string',
                              'check_with': check_cancer_study_identifier,
                              'required': True},
                          'stable_id':
                             {'type': 'string',
                              # To check if stable id follows the right format
                              'regex': re.escape(meta_study['cancer_study_identifier']) + r'_sequenced',
                              'required': True},
                          'case_list_name':
                             {'type': 'string',
                              'required': True},
                          'case_list_description':
                             {'type': 'string',
                              'required': True},
                          'case_list_ids':
                             {'type': 'string',
                              # To check if it's a tab-delimited list with no embedded whitespaces
                              'regex': r'^[^\t\s]+(?:\t[^\t\s]+)*$',
                              'check_with': check_duplicate_sample_ids,
                              'required': True},
                          'case_list_category':
                             {'type': 'string',
                              'allowed': ['all_cases_in_study',
                                'all_cases_with_mutation_data',
                                'all_cases_with_cna_data',
                                'all_cases_with_log2_cna_data',
                                'all_cases_with_methylation_data',
                                'all_cases_with_mrna_array_data',
                                'all_cases_with_mrna_rnaseq_data',
                                'all_cases_with_rppa_data',
                                'all_cases_with_microrna_data',
                                'all_cases_with_mutation_and_cna_data',
                                'all_cases_with_mutation_and_cna_and_mrna_data',
                                'all_cases_with_gsva_data',
                                'all_cases_with_sv_data',
                                'other'],
                              'required': False},
                         }

# Validation of meta file and case list
v = Validator(meta_mutations_schema)
if v.validate(meta_mutations) != True:
    print(v.errors)
else:
    print('Validation of meta file complete')
    
v = Validator(cases_sequenced_schema)
if v.validate(cases_sequenced) != True:
    print(v.errors)
else:
    print('Validation of meta file complete')

