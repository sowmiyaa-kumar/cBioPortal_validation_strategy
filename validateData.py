import os
import pandas as pd
import pandera as pa 
from pandera import Check, Column, DataFrameSchema, Index, MultiIndex
from pandera.errors import SchemaError
import pandera_schemas
from pandera_schemas import mut_schema  
from pydantic import BaseModel, ValidationError, validator, root_validator
import pydantic_schemas
from pydantic_schemas import MutData

def parse_file_to_dataframe(file_path):
    data = pd.read_csv(file_path, sep='\t', comment='#', header=0)
    return data

def detect_and_replace_missing_values(df):
    # Defining the list of missing values (in lower case) for each datatype
    missing_strings = ['unknown', 'n/a', 'na', 'null', '.', '', '?', '[not available]','[not applicable]', '[pending]', '[discrepancy]', '[completed]', '[null]']
    missing_numbers = [0, 0.0] # Add any dataset-specific missing numbers

    for col in df.columns:
        if df[col].dtype == object:  # if the column is of object type
            df.loc[df[col].str.lower().isin(missing_strings), col] = pd.NA
        elif df[col].dtype in [int, np.int64, float, np.float64]:  # if the column is of numeric type
            df.loc[df[col].isin(missing_numbers), col] = pd.NA

    return df

# def preprocessing_data(study_dir, data_files):
#     data = {}
#     for data_file in data_files:
#         data_df = pd.read_csv(os.path.join(study_dir, data_file))
#         data_df = detect_and_replace_missing_values(df)
    

def pandera_validation(data_df) -> None:
    # Validated data against schema
    try: 
        pandera_schemas.mut_schema.validate(data_df, lazy=True)
    except pa.errors.SchemaErrors as err:
        failure_cases = err.failure_cases
        error_df = err.data

    failure_cases_sorted = failure_cases.sort_values(by=['schema_context','column', 'index']).reset_index()
    # failure_cases_grouped = failure_cases_sorted.groupby(['column','check','schema_context'])['failure_case','index'].agg(list).reset_index()
    failure_cases_sorted.to_csv("errors/pandera/failure_cases.txt", sep = "\t")
    error_df.to_csv("errors/pandera/errors.txt", sep = "\t")

def pydantic_validation(data_df) -> None:
    with open("errors/pydantic/errors.txt", "w") as file:
        for idx, row in data_df.iterrows():
            try:
                pydantic_schemas.MutData(**row.to_dict())
            except ValidationError as e:
                file.write(f'Error in row {idx}: {e}\n\n')