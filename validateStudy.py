import os
import os.path as osp
import argparse
import validateStructure
import validateMeta
#import validateData
import pandas as pd
import logging

def validate_study(input_dir: str) -> None:
    # First level of validation - validate the directory structure
    meta_files, data_files = validateStructure.validate_directory(input_dir)

    # Second level of validation - validate the meta files 
    meta = validateMeta.parse_metadata(input_dir, meta_files)
    print(validateMeta.validate_metadata(meta))
    
    # Third level of validation - validate the data files
    #data_df = pd.read_csv("sample_data/brca_jup_msk_2020/data_mutations.txt", sep='\t', comment='#', header=0)
    #validateData.pandera_validation(data_df)
    #validateData.pydantic_validation(data_df)

    

if __name__ == '__main__':
    # Usage example: python3 validateStudy.py -i data/

    parser = argparse.ArgumentParser(description="Transforms all files for all studies in input folder to cBioPortal "
                                                 "staging files")
    
    parser.add_argument("-i", "--input_dir",
                        required=True,
                        help="Directory containing input files.")

    args = parser.parse_args()
    
    # Set up a logger
    logging.basicConfig(level=logging.DEBUG)

    validate_study(input_dir=args.input_dir)
