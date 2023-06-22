#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import json 

# Read the gene table into a dataframe 
genes_api = pd.read_json('http://cbioportal.org/api/genes')

# Read the chromosome JSON
chrom_sizes = {"hg19": {"1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260, "6": 171115067, "7": 159138663, "X": 155270560, "8": 146364022, "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753, "17": 81195210, "18": 78077248, "20": 63025520, "Y": 59373566, "19": 59128983, "22": 51304566, "21": 48129895}, "hg38": {"1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259, "6": 170805979, "7": 159345973, "X": 156040895, "8": 145138636, "9": 138394717, "11": 135086622, "10": 133797422, "12": 133275309, "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345, "17": 83257441, "18": 80373285, "20": 64444167, "19": 58617616, "Y": 57227415, "22": 50818468, "21": 46709983}, "mm10": {"1": 195471971, "2": 182113224, "X": 171031299, "3": 160039680, "4": 156508116, "5": 151834684, "6": 149736546, "7": 145441459, "10": 130694993, "8": 129401213, "14": 124902244, "9": 124595110, "11": 122082543, "13": 120421639, "12": 120129022, "15": 104043685, "16": 98207768, "17": 94987271, "Y": 91744698, "18": 90702639, "19": 61431566}}


# In[ ]:


# Use Pydantic to perform in-depth validation on individual rows 
from pydantic import BaseModel, ValidationError, validator, root_validator
from typing import Optional
import pydantic
import warnings
import requests
# Objects are defined via models in Pydantic
# Models are classes that inherit from the BaseModel

# Function to check if the required columns are present
def keys_exist(values, keys):
    for key in keys:
        if key not in values or values[key] is None:
            return False
    return True


'''
# To check the alias table
def check_gene_alias(v):
    response = requests.get(f"http://www.cbioportal.org/api/genes/{v}/aliases")
    if response.status_code == 200:
        return True  # Gene exists
    elif response.status_code == 404:
        return False  # Gene does not exist'''

class MutData(BaseModel):
    Hugo_Symbol: str 
    Entrez_Gene_Id: float | None
    NCBI_Build: str | None
    Chromosome: str | None
    Start_Position: int | None
    End_Position: int | None
    Variant_Classification: str
    Variant_Type: str | None
    Reference_Allele: str | None
    Tumor_Seq_Allele1: str | None
    Tumor_Seq_Allele2: str | None
    Tumor_Sample_Barcode: str
    Matched_Norm_Sample_Barcode: str | None
    Tumor_Validation_Allele1: str | None
    Tumor_Validation_Allele2: str | None
    Match_Norm_Validation_Allele1: str | None
    Match_Norm_Validation_Allele2: str | None
    Verification_Status: str | None
    Validation_Status: str | None
    Mutation_Status: str | None
    Validation_Method: str | None
    t_ref_count: int | None
    t_alt_count: int | None
    n_ref_count: int | None
    n_alt_count: int | None
    HGVSp_Short: str | None
    SWISSPROT: str | None
        
    # Custom validation can be carried out using the "validator" decorator
    """Since a validator works as a class method, the first argument is the 
    class which is named cls as a convention. The second argument is the 
    value of the field being validated, it can be named as you please."""   
    
    @validator("*", pre = True)
    def remove_whitespaces(cls, value):
        if isinstance(value, str):
            return value.strip()
        return value
    
    # Hugo_Symbol checks
    @validator('Hugo_Symbol')
    @classmethod
    def not_start_with_int(cls, value):
        if value.startswith(tuple(map(str, range(10))))==True:
            raise ValueError(f"WARNING - Hugo_Symbol should not start with a number.")
        return value
    
    @validator('Hugo_Symbol')   
    @classmethod
    def validate_hugo_symbol(cls, value):
#        genes_api_response = requests.get(f"http://cbioportal.org/api/genes/{value.upper()}")
#        alias_api_response = requests.get(f"http://www.cbioportal.org/api/genes/{value.upper()}/aliases")
        if pd.notna(value) or value is None:
            if value.upper() not in genes_api["hugoGeneSymbol"].values:
                raise ValueError(f"WARNING - {value} is not known to the cBioPortal instance. Might be new or deprecated gene symbol.")
        else:
            raise ValueError(f"WARNING - Hugo Gene Symbol is missing for this record.")
        return value
    
    # Entrez_Gene_Id checks
    @validator('Entrez_Gene_Id') 
    @classmethod
    def validate_entrez_gene_id(cls, value):
#        genes_api_response = requests.get(f"http://cbioportal.org/api/genes/{int(value)}")
#        alias_api_response = requests.get(f"http://www.cbioportal.org/api/genes/{int(value)}/aliases")
        if pd.notna(value) or value is None:
            if int(value) not in genes_api["entrezGeneId"].values:
                raise ValueError(f"WARNING - {int(value)} is not known to the cBioPortal instance. Might be new or deprecated Entrez gene id.") 
        else:
            raise ValueError(f"WARNING - Entrez Gene Id is missing for this record.")
        return value
    
    # Tumor_Sample_Barcode checks 
    @validator('Tumor_Sample_Barcode')
    @classmethod
    def validate_tumor_sample_barcode(cls, value):
        if value not in SAMPLE_IDS:
            raise ValueError(f"ERROR - Sample ID not defined in clinical file.")
        return value
                       
    # Checks involving multiple columns
    @root_validator(pre = False)
    def resolve_symbol_entrez(cls, values):
        required_fields = ['Hugo_Symbol', 'Entrez_Gene_Id']
        
        if keys_exist(values, ['Hugo_Symbol']) is False and keys_exist(values, ['Entrez_Gene_Id']) is False:
            raise ValueError(f"WARNING - Both gene identifiers for this gene are not valid. This record will not be loaded.")
            
        elif keys_exist(values, ['Hugo_Symbol']) is False and keys_exist(values, ['Entrez_Gene_Id']) is True:
            raise ValueError(f"WARNING - Entrez gene id exists, but gene symbol specified is not known to cBioPortal. The gene symbol will be ignored.")
            
        elif keys_exist(values, ['Hugo_Symbol']) is True and keys_exist(values, ['Entrez_Gene_Id']) is True:
            hugo_symbol = values.get('Hugo_Symbol')
            entrez_id = values.get('Entrez_Gene_Id')
#            gene_api_response = requests.get(f"http://cbioportal.org/api/genes/{hugo_symbol.upper()}").json()

            matching_entrez_gene_id = genes_api.loc[genes_api['hugoGeneSymbol'] == hugo_symbol.upper()]['entrezGeneId'].iloc[0]

            if int(entrez_id) != matching_entrez_gene_id:
                raise ValueError(f"ERROR - {entrez_id} does not match any valid entrezGeneId for {hugo_symbol}.")       
        return values
    
    # Default values - pre=False
    @root_validator(pre = False)
    def skip_variant(cls, values):
        required_fields = ['Hugo_Symbol', 'Entrez_Gene_Id', 'Variant_Classification']
        
        if keys_exist(values, required_fields):
            if pd.isna(values['Hugo_Symbol']) and pd.isna(values['Entrez_Gene_Id']):
                is_silent = True 
                if values['Variant_Classification'] in ['IGR', 'Targeted_Region']:
                    raise ValueError(f"INFO - This variant (Gene Symbol NA, Entrez Gene Id NA) will be filtered out.")
                else:
                    raise ValueError(f"WARNING - Gene specification (Gene symbol NA, Entrez gene ID NA) for this variant \
                                     implies intergenic even though Variant_Classification is \
                                     not 'IGR' or 'Targeted_Region'; this variant will be filtered out.")
            elif values['Variant_Classification'] in SKIP_VARIANT_TYPES:
                raise ValueError(f"INFO - Line will not be loaded due to the variant classification filter.")

        return values
    
    @root_validator(skip_on_failure = False)
    def non_splice_sites(cls, values):
        required_fields = ['Variant_Classification', 'HGVSp_Short']
        
        # check if a non-blank amino acid change exists for non-splice sites
        if keys_exist(values, required_fields):
            if values["Variant_Classification"] not in ['Splice_Site']:
                if pd.isna(values["HGVSp_Short"]):
                    raise ValueError(f'No HGVSp_Short value. This mutation record will get \
                                     a generic "MUTATED" flag.')
        
        return values
    
    @root_validator(skip_on_failure = False)
    def maf_check_6(cls, values):
        required_fields = ['Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
        
        # Check if Allele columms only contain the following values: -,A,C,T,G        
        # Else return an error message
        if keys_exist(values, required_fields):
            if not re.match(r'^[-ACTG]*$', values["Reference_Allele"]) and \
            not re.match(r'^[-ACTG]*$', values["Tumor_Seq_Allele1"]) and \
            not re.match(r'^[-ACTG]*$', values["Tumor_Seq_Allele2"]):
                raise ValueError(f"ERROR - All Allele Based columns contain invalid character.")
        return values 
    
    @root_validator(skip_on_failure = False)
    def maf_check_10(cls, values):
        required_fields = ['Start_Position', 'End_Position']
        
        if keys_exist(values, required_fields):
            # Check if Start_Position is equal or smaller than End_Position
            if not (values["Start_Position"] <= values["End_Position"]):
                raise ValueError(f"ERROR - Start_Position should be smaller than or equal to End_Position.")
        return values
            
    @root_validator(skip_on_failure = False)
    def maf_check_11(cls, values):
        required_fields = ['Variant_Type', 'End_Position', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
        
        if keys_exist(values, required_fields):
            if values["Variant_Type"] == "INS":
            # Start and End Position should be the same as length reference allele or equal 1
            # and length of Reference_Allele should be equal or smaller than the Tumor_Seq_Allele1 or 2,
            # otherwise no insertion, but deletion
                if not (((values["End_Position"] - values["Start_Position"] + 1) == len(values["Reference_Allele"])) or ((values["End_Position"] - values["Start_Position"]) == 1)):
                    raise ValueError(f"ERROR - Variant_Type indicates insertion, but difference in Start_Position and \
                                     End_Position does not equal to 1 or the length or the Reference_Allele.")
                if not (len(values["Reference_Allele"]) <= len(values["Tumor_Seq_Allele1"])) or (len(values["Reference_Allele"]) <= len(values["Tumor_Seq_Allele2"])):
                    raise ValueError(f"ERROR - Variant_Type indicates insertion, but length of Reference_Allele is bigger than \
                                     the length of the Tumor_Seq_Allele1 and/or 2 and therefore indicates deletion.")

            if values["Variant_Type"] == "DEL":
                # The difference between Start_Position and End_Position for a DEL should be equal to the length
                # of the Reference_Allele
                if not (values["End_Position"] - values["Start_Position"] + 1) == len(values["Reference_Allele"]):
                    raise ValueError(f"ERROR - Variant_Type indicates deletion, but the difference between Start_Position and \
                                     End_Position are not equal to the length of the Reference_Allele.")
                # The length of the Reference_Allele should be bigger than of the Tumor_Seq_Alleles for a DEL
                if (len(values["Reference_Allele"]) < len(values["Tumor_Seq_Allele1"])) or (len(values["Reference_Allele"]) < len(values["Tumor_Seq_Allele2"])):
                    raise ValueError(f"ERROR - Variant_Type indicates deletion, but length of Reference_Allele is smaller than \
                                     the length of Tumor_Seq_Allele1 and/or Tumor_Seq_Allele2, indicating an insertion.")

            if values["Variant_Type"] == "SNP":
                # Expect alleles to have length 1 when variant type is SNP
                if not (len(values["Reference_Allele"]) == 1 and len(values["Tumor_Seq_Allele1"]) == 1 and len(values["Tumor_Seq_Allele2"]) == 1):
                    raise ValueError(f"ERROR - Variant_Type indicates a SNP, but length of Reference_Allele, Tumor_Seq_Allele1 \
                                     and/or Tumor_Seq_Allele2 do not equal 1.")

            if values["Variant_Type"] == "DNP":
                # Expect alleles to have length 2 when variant type is DNP
                if not (len(values["Reference_Allele"]) == 2 and len(values["Tumor_Seq_Allele1"]) == 2 and len(values["Tumor_Seq_Allele2"]) == 2):
                    raise ValueError(f"ERROR - Variant_Type indicates a DNP, but length of Reference_Allele, Tumor_Seq_Allele1 \
                                      and/or Tumor_Seq_Allele2 do not equal 2.")

            if values["Variant_Type"] == "TNP":
                # Expect alleles to have length 3 when variant type is TNP
                if not (len(values["Reference_Allele"]) == 3 and len(values["Tumor_Seq_Allele1"]) == 3 and len(values["Tumor_Seq_Allele2"]) == 3):
                    raise ValueError(f"ERROR - Variant_Type indicates a TNP, but length of Reference_Allele, Tumor_Seq_Allele1 \
                                     and/or Tumor_Seq_Allele2 do not equal 3.")

            if values["Variant_Type"] == "ONP":
                # Expect alleles to have length >3 when variant type is ONP and are of equal length
                if (len(values["Reference_Allele"]) != len(values["Tumor_Seq_Allele1"])) or (len(values["Tumor_Seq_Allele1"]) != len(values["Tumor_Seq_Allele2"])) \
                        or (len(values["Reference_Allele"]) <= 3 and len(values["Tumor_Seq_Allele1"]) <= 3 and len(values["Tumor_Seq_Allele2"]) <= 3):
                    raise ValueError(f"ERROR - Variant_Type indicates a ONP, but length of Reference_Allele, \
                                     Tumor_Seq_Allele1 and 2 are not bigger than 3 or are of unequal lengths.")

            # Following variant types cannot contain a deletion in the Allele columns
            if values["Variant_Type"] == "SNP" or values["Variant_Type"] == "DNP" or values["Variant_Type"] == "TNP" or values["Variant_Type"] == "ONP":
                if ("-" in values["Reference_Allele"]) or ("-" in values["Tumor_Seq_Allele1"]) or ("-" in values["Tumor_Seq_Allele2"]):
                    raise ValueError(f'ERROR - Variant_Type indicates a {values["Variant_Type"]}, but Reference_Allele, Tumor_Seq_Allele1 \
                                     and/or Tumor_Seq_Allele2 contain deletion (-).')
        
        return values
        
    @root_validator(skip_on_failure = False)
    def checkAlleleSpecialCases(cls, values):
        """ Check other special cases which should or should not occur in Allele Based columns
        Special cases are either from unofficial vcf2maf rules or discrepancies identified. """
        required_fields = ['Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2']
        
        if keys_exist(values, required_fields):
            # Check if Allele Based columns are not all the same
            if values["Reference_Allele"] == values["Tumor_Seq_Allele1"] and values["Tumor_Seq_Allele1"] == values["Tumor_Seq_Allele2"]:
                raise ValueError(f"ERROR - All Values in columns Reference_Allele, Tumor_Seq_Allele1 \
                                 and Tumor_Seq_Allele2 are equal.")

            # In case of deletion, check when Reference_Allele is the same length as both Tumor_Seq_Allele if at least
            # one of the Tumor_Seq_Alleles is a deletion ('-') otherwise a SNP
            if values["Variant_Type"] == "DEL" and len(values["Reference_Allele"]) == len(values["Tumor_Seq_Allele1"]) \
            and len(values["Reference_Allele"]) == len(values["Tumor_Seq_Allele2"]) and "-" not in values["Tumor_Seq_Allele1"] \
            and "-" not in values["Tumor_Seq_Allele2"]:
                raise ValueError(f"ERROR - Variant_Type indicates a deletion, Allele based columns are the same length, \
                                  but Tumor_Seq_Allele columns do not contain -, indicating a SNP.")

        return values 

    """ Perform MAF file check #7,#8, #9 and #13 for the Validation columns:
    https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification)
    """
    @root_validator(skip_on_failure = False)
    def maf_check_7_and_8(cls, values):
        required_fields = ['Validation_Status', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2']
        
        if keys_exist(values, required_fields):
            # Check #7 and #8: When validation status is valid or invalid the Validation_Allele columns cannot be null
            validation_status = values["Validation_Status"]
            tumor_allele1 = values["Tumor_Validation_Allele1"]
            tumor_allele2 = values["Tumor_Validation_Allele2"]
            norm_allele1 = values["Match_Norm_Validation_Allele1"]
            norm_allele2 = values["Match_Norm_Validation_Allele2"]

            if validation_status.lower() == "valid" or validation_status.lower() == "invalid":
                if pd.isna(tumor_allele1) or pd.isna(tumor_allele2) or pd.isna(norm_allele1) or pd.isna(norm_allele2):
                    raise ValueError("ERROR - Validation Status is %s, but Validation Allele columns are empty. \
                                      % validation_status")

                # Check when Allele Columns are not empty if the Allele Based columns contain either -,A,C,T,G
                elif not re.match(tumor_allele1, r'^[-ACTG]*$') and \
                    not re.match(tumor_allele2, r'^[-ACTG]*$') and \
                    not re.match(norm_allele1, r'^[-ACTG]*$') and \
                    not re.match(norm_allele2, r'^[-ACTG]*$'):
                    raise ValueError(f"ERROR - At least one of the Validation Allele Based columns (Tumor_Validation_Allele1, \
                                     Tumor_Validation_Allele2, Match_Norm_Validation_Allele1, \
                                     Match_Norm_Validation_Allele2) contains invalid character.")

            # And in case of invalid also checks if validation alleles from tumor and normal match
            elif validation_status.lower() == "invalid" \
                    and (tumor_allele1 != norm_allele1 or tumor_allele2 != norm_allele2):
                raise ValueError(f"ERROR - When Validation_Status is invalid, the Tumor_Validation_Allele \
                                 and Match_Norm_Validation_Allele columns should be equal.")
        
        return values

    @root_validator(skip_on_failure = False)
    def maf_check_13(cls, values):
        required_fields = ['Validation_Status', 'Validation_Method']
        
        if keys_exist(values, required_fields):
            validation_status = values["Validation_Status"]
            validation_method = values["Validation_Method"]

            # Check #13: When Validation_Status is valid or invalid the Validation_Method column cannot be None.
            if validation_status.lower() == "valid" or validation_status.lower() == "invalid":
                if validation_method.lower() == "none" or validation_method.lower() == "na":
                    raise ValueError(f"ERROR - Validation Status is %s, but Validation_Method is not defined." %validation_status)
        return values

    @root_validator(skip_on_failure = False)
    def maf_check_9(cls, values):
        required_fields = ['Mutation_Status', 'Validation_Status', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', \
                          'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2']
        
        if keys_exist(values, required_fields):
            mutation_status = values["Mutation_Status"]
            validation_status = values["Validation_Status"]
            tumor_allele1 = values["Tumor_Validation_Allele1"]
            tumor_allele2 = values["Tumor_Validation_Allele2"]
            norm_allele1 = values["Match_Norm_Validation_Allele1"]
            norm_allele2 = values["Match_Norm_Validation_Allele2"]
            ref_allele = values["Reference_Allele"]
            
            # Check #9: Check Validation_Status against Mutation_Status
            # If Mutation_Status is Germline and Validation_Status is valid then the Tumor_Validation_Allele should be
            # equal to the matched Norm_Validation_Allele
            if mutation_status.lower() == "germline" and validation_status.lower() == "valid":
                if tumor_allele1 != norm_allele1 or tumor_allele2 != norm_allele2:
                    raise ValueError(f"ERROR - When Validation_Status is valid and Mutation_Status is Germline, the \
                                   Tumor_Validation_Allele should be equal to the Match_Norm_Validation_Allele.")

            # When Mutation_Status is Somatic and Validation_Status is valid the Norm_Validation Alleles should be equal
            # to the reference alleles and the Tumor_Validation_Alleles should be different from the Reference_Allele
            elif mutation_status.lower() == "somatic" and validation_status.lower() == "valid":
                if (norm_allele1 != norm_allele2 or norm_allele2 != ref_allele) and \
                        (tumor_allele1 != ref_allele or tumor_allele2 != ref_allele):
                    raise ValueError(f"ERROR - When Validation_Status is valid and Mutation_Status is Somatic, the \
                                   Match_Norm_Validation_Allele columns should be equal to the Reference Allele and \
                                   one of the Tumor_Validation_Allele columns should not be.")
        return values
    


# In[ ]:


# print(MutData.schema_json(indent=2))
# add 1 (or no. of header lines) to the row no (idx) to get the line number 
with open("prototype/pydantic/errors.txt", "w") as file:
    for idx, row in mut_data.iterrows():
        try:
            MutData(**row.to_dict())
        except ValidationError as e:
            file.write(f'Error in row {idx}: {e}\n\n')

