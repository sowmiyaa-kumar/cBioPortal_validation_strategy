import cerberus

# Schemas for meta files 
""" Define a schema in the form of nested dictionary 
(more information on https://docs.python-cerberus.org/en/stable/schemas.html)

Can view all the available validation rules at 
https://docs.python-cerberus.org/en/stable/validation-rules.html."""

# Custom validation rules 
# defining a custom validation rule to check consistency of cancer_study_identifiers across meta files 
# can also do this by extending the Validator class
def check_cancer_study_identifier(field, value, error):
    if value != meta['STUDY']['cancer_study_identifier']:
        error(field, "The cancer study identifier does not match to that of meta_study.txt.")
        
def check_duplicate_sample_ids(field, value, error):
    sample_ids = [x.strip() for x in cases_sequenced['case_list_ids'].split('\t')]
    if len(sample_ids) != len(set(sample_ids)):
        error(field, "Duplicate sample IDs in case list.")

META_SCHEMA_MAP = {
    "CANCER_TYPE": {'genetic_alteration_type': 
                        {'type': 'string',
                        'allowed': ['CANCER_TYPE'],
                        'required': True},
                    'datatype': 
                        {'type': 'string',
                        'allowed': ['CANCER_TYPE'],
                        'required': True},
                    'data_filename':
                        {'type': 'string',
                         'required': True},
                    },
    "STUDY": {'cancer_study_identifier':
                {'type': 'string',
                 'maxlength': 255,
                 'required': True},
            'type_of_cancer': 
                {'type': 'string',
                 'maxlength': 63,
                 'required': True},
            'name': 
                {'type': 'string',
                 'maxlength': 255,
                 'required': True},
            'description': 
                {'type': 'string',
                 'maxlength': 1024,
                 'required': True},
            'citation': 
                {'type': 'string',
                 'maxlength': 200,
                 'required': False},
            'pmid': 
                {'type': 'string',
                 'maxlength': 1024,
                 'required': False},
            'groups': 
                {'type': 'string',
                 'maxlength': 200,
                 'required': False},
            'short_name': 
                {'type': 'string',
                 'maxlength': 64,
                 'required': False},
            'add_global_case_list': 
                {'type': 'string',
                 'required': False},
            'tags_file': 
                {'type': 'string',
                 'required': False},
            'reference_genome': 
                {'type': 'string',
                 'required': False},
            },
    "SAMPLE_ATTRIBUTES": {'cancer_study_identifier':
                            {'type': 'string',
                             'maxlength': 255, 
                             'required': True},
                        'genetic_alteration_type':
                            {'type': 'string',
                             'allowed': ['CLINICAL'], 
                             'required': True},
                        'datatype':
                            {'type': 'string',
                             'allowed': ['SAMPLE_ATTRIBUTES'],
                             'required': True},
                        'data_filename':
                            {'type': 'string',
                             'required': True},
                        },
    "PATIENT_ATTRIBUTES": {'cancer_study_identifier':
                            {'type': 'string',
                             'maxlength': 255, 
                             'required': True},
                        'genetic_alteration_type':
                            {'type': 'string',
                             'allowed': ['CLINICAL'], 
                             'required': True},
                        'datatype':
                            {'type': 'string',
                             'allowed': ['PATIENT_ATTRIBUTES'],
                             'required': True},
                        'data_filename':
                            {'type': 'string',
                             'required': True},
                        },
    "CNA_DISCRETE": {'cancer_study_identifier':
                        {'type': 'string',
                         'maxlength': 255, 
                         'required': True},
                    'genetic_alteration_type':
                        {'type': 'string',
                         'allowed': ['COPY_NUMBER_ALTERATION'], 
                         'required': True},
                    'datatype':
                        {'type': 'string',
                         'allowed': ['DISCRETE'],
                         'required': True},
                    'stable_id':
                        {'type': 'string',
                         'allowed': ['gistic', 'cna', 'cna_rae', 'cna_consensus'],
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
                    'gene_panel':
                        {'type': 'string',
                         'required': False},
                    'pd_annotations_filename':
                        {'type': 'string',
                         'required': False},
                    },
    "CNA_DISCRETE_LONG": {'cancer_study_identifier':
                        {'type': 'string',
                         'maxlength': 255, 
                         'required': True},
                    'genetic_alteration_type':
                        {'type': 'string',
                         'allowed': ['COPY_NUMBER_ALTERATION'], 
                         'required': True},
                    'datatype':
                        {'type': 'string',
                         'allowed': ['DISCRETE_LONG'],
                         'required': True},
                    'stable_id':
                        {'type': 'string',
                         'allowed': ['gistic', 'cna', 'cna_rae', 'cna_consensus'],
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
                    'gene_panel':
                        {'type': 'string',
                         'required': False},
                    'namespaces':
                        {'type': 'string',
                         'regex': r'^[^,\s]+(?:,[^,\s]+)*$',
                         'required': False},
                    },
    "CNA_LOG2": {'cancer_study_identifier':
                    {'type': 'string',
                     'maxlength': 255, 
                     'required': True},
                'genetic_alteration_type':
                    {'type': 'string',
                     'allowed': ['COPY_NUMBER_ALTERATION'], 
                     'required': True},
                'datatype':
                    {'type': 'string',
                     'allowed': ['DISCRETE_LONG'],
                     'required': True},
                'stable_id':
                    {'type': 'string',
                     'allowed': ['gistic', 'cna', 'cna_rae', 'cna_consensus'],
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
                'gene_panel':
                    {'type': 'string',
                        'required': False},
    # MetaFileTypes.CNA_LOG2: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'show_profile_in_analysis_tab': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'gene_panel': False
    # },
    # MetaFileTypes.CNA_CONTINUOUS: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'show_profile_in_analysis_tab': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'gene_panel': False
    # },
    "SEG": {'cancer_study_identifier':
                {'type': 'string',
                 'maxlength': 255, 
                # 'check_with': check_cancer_study_identifier,
                 'required': True},
            'genetic_alteration_type':
                {'type': 'string',
                 'allowed': ['COPY_NUMBER_ALTERATION'], 
                 'required': True},
            'datatype':
                {'type': 'string',
                 'allowed': ['SEG'],
                 'required': True},
            'reference_genome_id': 
                {'type': 'string',
                 'required': True},
            'data_filename':
                {'type': 'string',
                 'required': True},
            'description': 
                {'type': 'string',
                 'maxlength': 1024,
                 'required': True},
    },
    "MUTATION": {'cancer_study_identifier':
                    {'type': 'string',
                     'maxlength': 255, 
                    # 'check_with': check_cancer_study_identifier,
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
                     'allowed': ['true', 'TRUE'],
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
                },
    # MetaFileTypes.MUTATION_UNCALLED: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'show_profile_in_analysis_tab': False,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'normal_samples_list': False,
    #     'swissprot_identifier': False,
    #     'gene_panel': False,
    #     'variant_classification_filter': False,
    #     'namespaces': False
    # },
    # MetaFileTypes.EXPRESSION: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'source_stable_id': False,
    #     'show_profile_in_analysis_tab': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'gene_panel': False
    # },
    # MetaFileTypes.METHYLATION: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'show_profile_in_analysis_tab': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'gene_panel': False
    # },
    # MetaFileTypes.PROTEIN: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'show_profile_in_analysis_tab': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'gene_panel': False
    # },
    # MetaFileTypes.GISTIC_GENES: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'reference_genome_id': True,
    #     'data_filename': True
    # },
    # MetaFileTypes.TIMELINE: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'data_filename': True
    # },
    # MetaFileTypes.CASE_LIST: {
    #     'cancer_study_identifier': True,
    #     'stable_id': True,
    #     'case_list_name': True,
    #     'case_list_description': True,
    #     'case_list_ids': True,
    #     'case_list_category': False # TODO this is used in org.mskcc.cbio.portal.model.AnnotatedPatientSets.getDefaultPatientList(), decide whether to keeep, see #494
    # },
    # MetaFileTypes.MUTATION_SIGNIFICANCE: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'data_filename': True
    # },
    # MetaFileTypes.GENE_PANEL_MATRIX: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'data_filename': True
    # },
    # MetaFileTypes.GSVA_PVALUES: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'source_stable_id': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'geneset_def_version': True
    # },
    # MetaFileTypes.GSVA_SCORES: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'source_stable_id': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'show_profile_in_analysis_tab': True,
    #     'geneset_def_version': True
    # },
    # MetaFileTypes.GENERIC_ASSAY_CONTINUOUS: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'generic_assay_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'show_profile_in_analysis_tab': True,
    #     'generic_entity_meta_properties': False,
    #     'pivot_threshold_value': False,
    #     'value_sort_order': False,
    #     'patient_level': False
    # },
    # MetaFileTypes.GENERIC_ASSAY_BINARY: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'generic_assay_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'show_profile_in_analysis_tab': True,
    #     'generic_entity_meta_properties': False,
    #     'patient_level': False
    # },
    # MetaFileTypes.GENERIC_ASSAY_CATEGORICAL: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'generic_assay_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'show_profile_in_analysis_tab': True,
    #     'generic_entity_meta_properties': False,
    #     'patient_level': False
    # },
    # MetaFileTypes.STRUCTURAL_VARIANT: {
    #     'cancer_study_identifier': True,
    #     'genetic_alteration_type': True,
    #     'datatype': True,
    #     'stable_id': True,
    #     'show_profile_in_analysis_tab': True,
    #     'profile_name': True,
    #     'profile_description': True,
    #     'data_filename': True,
    #     'gene_panel': False,
    #     'namespaces': False
    # },
    # MetaFileTypes.SAMPLE_RESOURCES: {
    #     'cancer_study_identifier': True,
    #     'resource_type': True,
    #     'data_filename': True
    # },
    # MetaFileTypes.PATIENT_RESOURCES: {
    #     'cancer_study_identifier': True,
    #     'resource_type': True,
    #     'data_filename': True
    # },
    # MetaFileTypes.STUDY_RESOURCES: {
    #     'cancer_study_identifier': True,
    #     'resource_type': True,
    #     'data_filename': True
    # },
    # MetaFileTypes.RESOURCES_DEFINITION: {
    #     'cancer_study_identifier': True,
    #     'resource_type': True,
    #     'data_filename': True
    # },
}


# Schema for cases_sequenced 
# cases_sequenced_schema = {'cancer_study_identifier':
#                              {'type': 'string',
#                               'check_with': check_cancer_study_identifier,
#                               'required': True},
#                           'stable_id':
#                              {'type': 'string',
#                               # To check if stable id follows the right format
#                               'regex': re.escape(meta_study['cancer_study_identifier']) + r'_sequenced',
#                               'required': True},
#                           'case_list_name':
#                              {'type': 'string',
#                               'required': True},
#                           'case_list_description':
#                              {'type': 'string',
#                               'required': True},
#                           'case_list_ids':
#                              {'type': 'string',
#                               # To check if it's a tab-delimited list with no embedded whitespaces
#                               'regex': r'^[^\t\s]+(?:\t[^\t\s]+)*$',
#                               'check_with': check_duplicate_sample_ids,
#                               'required': True},
#                           'case_list_category':
#                              {'type': 'string',
#                               'allowed': ['all_cases_in_study',
#                                 'all_cases_with_mutation_data',
#                                 'all_cases_with_cna_data',
#                                 'all_cases_with_log2_cna_data',
#                                 'all_cases_with_methylation_data',
#                                 'all_cases_with_mrna_array_data',
#                                 'all_cases_with_mrna_rnaseq_data',
#                                 'all_cases_with_rppa_data',
#                                 'all_cases_with_microrna_data',
#                                 'all_cases_with_mutation_and_cna_data',
#                                 'all_cases_with_mutation_and_cna_and_mrna_data',
#                                 'all_cases_with_gsva_data',
#                                 'all_cases_with_sv_data',
#                                 'other'],
#                               'required': False},
#                          }
