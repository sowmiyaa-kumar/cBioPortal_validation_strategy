#!/usr/bin/env python
# coding: utf-8

# import required packages 
import pandas as pd
import numpy as np
from copy import deepcopy
import re
import pandera as pa
from pandera import Check, Column, DataFrameSchema, Index, MultiIndex
from pandera.errors import SchemaError

# # Read mutations data into pandas dataframe
# mut_data = pd.read_csv(os.path.join(file_dir, "data_mutations.txt"), sep='\t', comment='#', header=0)
# mut_data_copy = deepcopy(mut_data)

# # Save basic statistics
# no_of_features = len(mut_data.columns) # contains 112 columns 
# no_of_genes = len(mut_data.index) # contains 11660 rows 

# # Pre-processing steps 
# # Replace all empty values with None/NaN (Pandera does not consider empty strings as missing data)
# # mut_data = mut_data.replace('', pd.NA)

# def detect_and_replace_missing_values(df):
#     # Defining the list of missing values (in lower case) for each datatype
#     missing_strings = ['unknown', 'n/a', 'na', 'null', '.', '', '?', '[not available]']
#     missing_numbers = [0, 0.0] # Add any dataset-specific missing numbers

#     for col in df.columns:
#         if df[col].dtype == object:  # if the column is of object type
#             df.loc[df[col].str.lower().isin(missing_strings), col] = pd.NA
#         elif df[col].dtype in [int, np.int64, float, np.float64]:  # if the column is of numeric type
#             df.loc[df[col].isin(missing_numbers), col] = pd.NA

#     return df

# # Usage:
# mut_data = detect_and_replace_missing_values(mut_data)

# Drop columns with all missing values
# mut_data.dropna(axis=1, how='all', inplace=True)

# Use Pandera for initial validation and data cleaning (table-wide and column-specific checks)
"""Initial validation may include:
    - Datatype checking
    - Built-in checks for each column 
    - Removing duplicates 
    - Defining column order 
    """

# Set up lists of columns/values that need to be checked for or validated against.
REQUIRED_HEADERS = [
    'Tumor_Sample_Barcode',
    'Hugo_Symbol',
    'Variant_Classification',
]

NULL_AA_CHANGE_VALUES = ('', 'NULL', 'NA')

EXTRA_VARIANT_CLASSIFICATION_VALUES = ['Splice_Region', 'Fusion']

SKIP_VARIANT_TYPES = [
    'Silent',
    'Intron',
    '3\'UTR',
    '3\'Flank',
    '5\'UTR',
    '5\'Flank',
    'IGR',
    'RNA'
]

VARIANT_CLASSIFICATION_VALUES = [
       'Frame_Shift_Del',
       'Frame_Shift_Ins',
       'In_Frame_Del',
       'In_Frame_Ins',
       'Missense_Mutation',
       'Nonsense_Mutation',
       'Splice_Site',
       'Translation_Start_Site',
       'Nonstop_Mutation',
       'Targeted_Region',
       'De_novo_Start_InFrame',
       'De_novo_Start_OutOfFrame'] + SKIP_VARIANT_TYPES + EXTRA_VARIANT_CLASSIFICATION_VALUES + ['Unknown']


REQUIRED_ASCN_COLUMNS = [
    'ASCN.ASCN_METHOD',
    'ASCN.ASCN_INTEGER_COPY_NUMBER',
    'ASCN.TOTAL_COPY_NUMBER',
    'ASCN.MINOR_COPY_NUMBER',
    'ASCN.CCF_EXPECTED_COPIES',
    'ASCN.CCF_EXPECTED_COPIES_UPPER',
    'ASCN.CLONAL',
    'ASCN.EXPECTED_ALT_COPIES'
]


# Define custom functions for checks that are not built-in 
# Table-wide custom checks
# Function to check atleast one gene identifier column in present
def at_least_one_gene_identifier(df):
    if 'Hugo_Symbol' not in df.columns and 'Entrez_Gene_Id' not in df.columns:
        return False 
    else:
        return True

# Function to check atleast one of 'HGVSp_Short' or 'Amino_Acid_Change' columns are present 
def at_least_one_aa_change_col(df):
    if 'HGVSp_Short' not in df.columns and 'Amino_Acid_Change' not in df.columns:
        return False 
    else:
        return True

# Function to check if 'SWISSPROT' column is present
def swissprot_in_data_and_meta(df):
    if 'SWISSPROT' not in df.columns:
        return False
    else:
        return True
'''    elif 'swissprot_identifier' in meta_mutations.keys():
        raise ValueError(f"WARNING - A SWISSPROT column was found in datafile without specifying \
        associated 'swissprot_identifier' in metafile, assuming \
        'swissprot_identifier: name'.")'''
        
# Function to check if 'ascn' namespace is defined and if yes, the required columns (defined above) are present    
def ascn_namespace_defined(df):
    if 'namespaces' in meta_mutations:
        namespaces = meta_mutations['namespaces'].split(',')
        for namespace in namespaces: 
            if 'ascn' == namespace.strip().lower():
                for required_ascn_column in REQUIRED_ASCN_COLUMNS:
                    if required_ascn_column not in df.columns:
                        return False
    return True

mut_schema = DataFrameSchema(
    columns={
        "Hugo_Symbol": Column(
            dtype="object",
            checks=[
                pa.Check(lambda x: x.startswith(tuple(map(str, range(10))))==False,
                         element_wise=True,
                         ignore_na=True,
                         error="ERROR - Hugo_Symbol should not start with a number."),
            ],
            nullable=True,
            unique=False,
            coerce=False,
            required=True,
            regex=False,
            description=None,
            title=None,
        ),
        "Entrez_Gene_Id": Column(
            dtype="float64",
            checks=[
                pa.Check(lambda x: x >= 0.0, 
                         element_wise=True,
                         ignore_na=True,
                         error="ERROR - Entrez gene is non-positive."), 
            ],
            nullable=True,
            unique=False,
            coerce=False,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "NCBI_Build": Column(
            dtype="object",
            checks=[
                Check.isin(["GRCh37", "GRCh38", "GRCm38", "37", "38"],
                         ignore_na=True,
                         error="ERROR - NCBI Build is not a valid cBioPortal build."),
            ],
            nullable=True,
            unique=False,
            coerce=False,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Chromosome": Column(
            dtype="object",
            checks=[
                pa.Check(lambda x: x in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "X"],
                         element_wise=True,
                         ignore_na=True, 
                         error="ERROR - Chromosome not found in the genome."),
            ],
            nullable=True,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Start_Position": Column(
            dtype="int64",
            checks=None,
            nullable=False,
            unique=False,
            # Originally float, convert to int(?)
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "End_Position": Column(
            dtype="int64",
            checks=None,
            nullable=False,
            unique=False,
            # Originally float, convert to int()
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Variant_Classification": Column(
            dtype="object",
            checks=[
                Check.isin(VARIANT_CLASSIFICATION_VALUES,
                           ignore_na = True, 
                           error = 'ERROR - Invalid value. Not in VARIANT_CLASSIFICATION_VALUES.'),
            ],
            nullable=True,
            unique=False,
            coerce=True,
            required=True,
            regex=False,
            description=None,
            title=None,
        ),
        "Variant_Type": Column(
            dtype="object",
            checks=None,
            nullable=True,
            unique=False,
            coerce=False,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Reference_Allele": Column(
            dtype="object",
            checks=[
                Check.str_matches(r'^[-ACTG]*$',
                                  ignore_na = True,
                                  error = "ERROR - Allele Based column Reference_Allele contains \
                                  invalid character."),
            ],
            nullable=True,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Tumor_Seq_Allele1": Column(
            dtype="object",
            checks=[
                Check.str_matches(r'^[-ACTG]*$',
                                  ignore_na = True,
                                  error = "ERROR - Allele Based column Tumor_Seq_Allele1 contains \
                                  invalid character."),
            ],
            nullable=True,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Tumor_Seq_Allele2": Column(
            dtype="object",
            checks=[
                Check.str_matches(r'^[-ACTG]*$',
                                  ignore_na = True,
                                  error = "ERROR - Allele Based column Tumor_Seq_Allele2 contains \
                                  invalid character."),
            ],
            nullable=True,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Tumor_Sample_Barcode": Column(
            dtype="object",
            checks=None,
            nullable=False,
            unique=False,
            coerce=False,
            required=True,
            regex=False,
            description=None,
            title=None,
        ),
        "Matched_Norm_Sample_Barcode": Column(
            dtype="object",
            checks=None,
            nullable=False,
            unique=False,
            coerce=False,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Tumor_Validation_Allele1": Column(
            dtype="object",
            checks=[
                Check.str_matches(r'^[-ACTG]*$',
                                  ignore_na = True,
                                  error = "ERROR - Allele Based column Tumor_Validation_Allele1 contains \
                                  invalid character."),
            ],
            nullable=False,
            unique=False,
            coerce=False,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Tumor_Validation_Allele2": Column(
            dtype="object",
            checks=[
                Check.str_matches(r'^[-ACTG]*$',
                                  ignore_na = True,
                                  error = "ERROR - Allele Based column Tumor_Validation_Allele2 contains \
                                  invalid character."),
            ],
            nullable=False,
            unique=False,
            coerce=False,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Match_Norm_Validation_Allele1": Column(
            dtype="object",
            checks=[
                Check.str_matches(r'^[-ACTG]*$',
                                  ignore_na = True,
                                  error = "ERROR - Allele Based column Match_Norm_Validation_Allele1 contains \
                                  invalid character."),
            ],
            nullable=False,
            unique=False,
            coerce=False,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Match_Norm_Validation_Allele2": Column(
            dtype="object",
            checks=[
                Check.str_matches(r'^[-ACTG]*$', 
                                  ignore_na = True,
                                  error = "ERROR - Allele Based column Match_Norm_Validation_Allele2 contains \
                                  invalid character."),
            ],
            nullable=False,
            unique=False,
            coerce=False,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Verification_Status": Column(
            dtype="object",
            checks=[
                pa.Check(lambda x: x.lower() in ["verified", "unknown", "na"],
                         element_wise = True,
                         ignore_na = True,
                         error = f"ERROR - Value in 'Verification_Status' not in MAF format."),
            ],
            nullable=True,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Validation_Status": Column(
            dtype="object",
            checks=[
                pa.Check(lambda x: x.lower() in ["untested", "inconclusive", "valid", "invalid", "na", "redacted", "unknown"],
                         element_wise = True, 
                         ignore_na = True, 
                         error = f"WARNING - Value in 'Validation_Status' not in MAF format."),
            ],
            nullable=True,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Mutation_Status": Column(
            dtype="object",
            checks=[
                pa.Check(lambda x: x.lower() in ["germline", "somatic", "post-transcriptional modification", "unknown"],
                         element_wise = True, 
                         ignore_na = True, 
                         error = 'WARNING - Mutation_Status value is not in MAF format'),
                pa.Check(lambda x: x.lower() not in ["loh", "none", "wildtype"],
                         element_wise = True,
                         ignore_na = True, 
                         error = 'INFO - Mutation will not be loaded due to value in Mutation_Status'),
            ],
            nullable=True,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Validation_Method": Column(
            dtype="object",
            checks=None,
            nullable=True,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "t_ref_count": Column(
            dtype="int64",
            checks=None,
            nullable=False,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "t_alt_count": Column(
            dtype="int64",
            checks=None,
            nullable=False,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "n_ref_count": Column(
            dtype="int64",
            checks=None,
            nullable=False,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "n_alt_count": Column(
            dtype="int64",
            checks=None,
            nullable=False,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "HGVSp_Short": Column(
            dtype="object",
            checks=[
                Check.str_length(min_value=0, 
                                 max_value=255, 
                                 ignore_na = True, 
                                 error = f'ERROR - cBioPortal does not support values longer than 255 characters in "HGVSp_Short."'),
            ],
            nullable=True,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "Exon_Number": Column(
            dtype="object",
            checks=None,
            nullable=True,
            unique=False,
            coerce=True,
            required=False,
            regex=False,
            description=None,
            title=None,
        ),
        "SWISSPROT": Column(
            dtype="object",
            checks=[
                pa.Check(lambda x: pd.notna(x),
                             element_wise = True,
                             ignore_na = False,
                             error = "WARNING - Missing value in SWISSPROT column; this column is recommended to make sure that the UniProt canonical isoform is used when drawning Pfam domains in the mutations view."),
            ],
            nullable=True,
            unique=False,
            coerce=False,
            required=True,
            regex=False,
            description=None,
            title=None,
        ),
    },
    checks=[
        pa.Check(lambda df: at_least_one_gene_identifier(df),
                 error = f"ERROR - At least one of the columns Hugo_Symbol or \
                 Entrez_Gene_Id needs to be present."),
        pa.Check(lambda df: at_least_one_aa_change_col(df),
                 error = f"ERROR - At least one of the columns HGVSp_Short or \
                 Amino_Acid_Change needs to be present."),
        pa.Check(lambda df: swissprot_in_data_and_meta(df),
                 error = f"WARNING - Including the SWISSPROT column is recommended to make sure that the UniProt canonical isoform is used when drawing Pfam domains in the mutations view."),
        pa.Check(lambda df: ascn_namespace_defined(df),
                 error = f"ERROR - ASCN namespace defined but MAF missing required ASCN columns."),
    ],
    index=Index(
        dtype="int64",
        checks=None,
        nullable=False,
        coerce=False,
        name=None,
        description=None,
        title=None,
    ),
    dtype=None,
    coerce=False,
    strict=False,
    name=None,
    # Column order not required
    ordered=False,
    unique=None,
    report_duplicates="all",
    # Column names should not be repeated
    unique_column_names=True,
    title=None,
    description=None,
)

# # Validated data against schema
# try: 
#     mut_schema.validate(mut_data, lazy=True)
# except pa.errors.SchemaErrors as err:
#     failure_cases = err.failure_cases
#     error_df = err.data

# failure_cases_sorted = failure_cases.sort_values(by=['schema_context','column', 'index']).reset_index()
# # failure_cases_grouped = failure_cases_sorted.groupby(['column','check','schema_context'])['failure_case','index'].agg(list).reset_index()
# failure_cases_sorted.to_csv("prototype/pandera/failure_cases.txt", sep = "\t")
# error_df.to_csv("prototype/pandera/errors.txt", sep = "\t")

