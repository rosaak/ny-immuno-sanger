#!/usr/bin/env python
# coding: utf-8

# # This script takes sequenced ABI files from VH and aligned with genbank file - the script filters and copy it to a new dir
# - requires utils.py and fns.py

import configparser
from collections import Counter
from itertools import chain
from pathlib import Path
from time import perf_counter

import pandas as pd
import utils
import xlsxwriter
from fns import *

### START LOGGING 
Path('logs').mkdir(parents=True, exist_ok=True)
log = utils.make_logger("logs/single_log")
log.info("START")
log.info("CHECKING CONFIGURATION FILE CONTENTS...")


### GET PARAMS FROM CONFIG FILE
# checking config file
# if any error occurs then program terminates; then check the config file
config = configparser.ConfigParser()
config.read_file(open('config_single_vh.ini'))


abbr = get_analysis_abbr(config)
abbr_lower = abbr.lower()
abi_sequence_folder = checking_dirs(config['Paths']['abi_sequence_folder'], log, log_msg=True, create_dir=False)
vx_template_sequence_folder = checking_dirs(config['Paths']['vx_template_sequence_folder'], log, log_msg=True, create_dir=False)
results_dir = checking_dirs(config['Paths']['results_dir'], log, log_msg=True, create_dir=True)
h3_nt_data_sheet_filepath = check_files(config['Files']['h3_nt_data_sheet_filepath'], log)
df = check_tsv_file(h3_nt_data_sheet_filepath, log)
excel_path_file_name = create_results_excel_file_path(config['Paths']['results_dir'], config['Files']['output_excel_file_name'])

# patterns
pat_vx = get_patterns(config)
# patterns to remove from abi names
patrm_vx_abi = get_patterns_to_rm(config)
# patterns to remove from genbank names
patrm_vx_gb = get_patterns_to_rm_from_genbank(config)
# get alignment parameters
par_match, par_missmatch, par_open, par_extend, par_filter_thresh = get_alignment_params(config)
# creating dirs for copying matched abi files
res_dir_vx = checking_dirs(f"{results_dir}/{pat_vx}", log, log_msg=True, create_dir=True)
# abi extension
abi_extension = get_abi_extension(config)
# genbak extension
gb_extension = get_gb_extension(config)
# genbank filter
gb_filter_pat = get_gb_pattern_to_filter(config)


### CHECKING ABI FILES
log.info("CHECKING THE AB1 FILES...")
vx_abi_dict = {i.name.replace(abi_extension, '').replace(patrm_vx_abi, '') : str(i) for i in sorted([*Path(abi_sequence_folder).glob(f"*{pat_vx}*{abi_extension}")])}
log.info(f"[+] INFO: There are {len(vx_abi_dict)} abi files in the ab1 sequence folder")
sample_ids = checking_vx_ab1_files(log, vx_abi_dict)
### CHECKING GENBANK FILES
log.info("CHECKING THE GENBANK FILES...")
vx_template_gb_dict =  {i.name.replace(gb_extension, '') : str(i) for i in sorted(Path(vx_template_sequence_folder).glob(f"*{gb_extension}"))}

#### Filter genbank files with pattern
if len(gb_filter_pat) >=1:
    vx_template_gb_dict_ss = {i[0]:i[1] for i in vx_template_gb_dict.items() if gb_filter_pat in i[0]}
else:
    vx_template_gb_dict_ss = vx_template_gb_dict
log.info(f"[+] INFO: There are {len(vx_template_gb_dict_ss)} genbank files in the dir")

#### Extracting 6bp nts from start and end of genbak files
vx_nts = extract_nts_from_start_and_end_from_genbank(vx_template_gb_dict_ss, 6)
vx_seq_start_a, vx_seq_end_a = vx_nts.start_nts.value_counts().nlargest(1).index[0], vx_nts.end_nts.value_counts().nlargest(1).index[0]
# vx_seq_start_a, vx_seq_end_a, vx_nts.start_nts.value_counts().values[0], vx_nts.end_nts.value_counts().values[0]

### GET PROBE SEQS INTO A DICTIOANRY
# dictionary containing probe seq and name
h3_dict = df.set_index('name', drop=True).to_dict().get('h3_nt')

### FINDING THE PROBES THAT MATCH VH
log.info(f"ITERATING THROUGH EACH SAMPLE ID")
result_vx= []
for sample in sample_ids:
    _vxabi = get_vx_abi_file_path(key=sample, vx_abi_dict=vx_abi_dict)
    _vx_d = get_seqobj_from_abi(_vxabi)  # returns a seq record obj of VX
    
    # matching each probe on  Vx - normal and revcomp sequence
    vh_prob_search = find_match_on_all_h3probes(log, h3_dict, _vx_d, sample, vx_abi_dict)
    
    if len(vh_prob_search) >=1:
        result_vx.append(vh_prob_search)
log.info(f"FINISH ITERATING THROUGH EACH SAMPLE ID")
colnames=["Match","h3_name","sample_id",f"{abbr_lower}_abi_fp", "probe_seq", f"{abbr_lower}_init_sr", f"{abbr_lower}_sr_seq_r", f"{abbr_lower}_sr_trimmed", f"{abbr_lower}_sr_tqlst"]
df_vx = pd.DataFrame(chain.from_iterable(result_vx))
df_vx.columns = colnames
log.info(f"[+] INFO: There are {df_vx.shape[0]} matches in df_vx")
### COPY PROBE MATCHED ABI FILES TO A NEW LOC
log.info(f"COPY MATCHED ABI FILES INTO NEW LOCATION")
res_df_copy = copy_mtched_abi_files_to_resdir_vx(log, res_dir_vx, abbr_lower, df_vx=df_vx, log_msg=True)
log.info(f"FINISH COPYING FILES")
### PAIRWISE ALIGNMENT OF GENBANK FILES WITH MATCHED ABI FILES TO GET THE MATCHING SCORE 
log.info(f"FIND THE GB FILENASMES WITH H3 NAMES")
M_gb_abi_vx = find_gb_match_on_all_h3probes_single(log, vx_template_gb_dict_ss, df_vx, pattern=patrm_vx_gb, log_msg=True)
log.info(f"ALIGNEMNT BETWEEN THE MATCHED GB FILENASME AND ABI FILES")
vx_gb_abi_match_filtered = run_gb_alignment_and_filtering_single(M_gb_abi_vx, df_vx, vx_template_gb_dict_ss, vx_seq_start_a, vx_seq_end_a, 
                                                                par_match, par_missmatch, par_open, par_extend, par_filter_thresh, log, 
                                                                log_msg=True, is_data_vl=False)
log.info(f"MERGING DATAFRAMES")
vx_gb_abi_match_filtered.columns = [ f"{abbr}_"+i  if i not in ['Orient','gbid','H3_name', 'sample_id', 'GB_FP'] else i  for i in vx_gb_abi_match_filtered.columns ]
df_sid_abi = pd.DataFrame([ [i, vx_abi_dict.get(i).split("/")[-1]] for i in vx_gb_abi_match_filtered.sample_id.to_list()], columns=["sample_id", f"{abbr}"])
dfx = pd.merge(vx_gb_abi_match_filtered,df_sid_abi,on='sample_id')

log.info(f"[+] INFO: There are {vx_gb_abi_match_filtered.shape[0]} rows in the GB - {abbr} aligned datafrme `vx_gb_abi_match_filtered`")
log.info(f"[+] INFO: There are {dfx.shape[0]} rows in the final - filtered merged cleanedup results")
log.info(f"[+] STAT: basic stats of the numerical columns in the final dataframe : \n{dfx[dfx.dtypes[dfx.dtypes !='object'].index].describe()}")
log.info(f"WRTING EXCEL FILES")

excel_fp = get_excel_file_name(config, accessory=False)
with pd.ExcelWriter(excel_fp, engine='xlsxwriter') as writer:
    dfx.to_excel(writer, sheet_name='final_res_mean_error_prob', index=False)
    M_gb_abi_vx.to_excel(writer, sheet_name=f"ID_matched_gb_and_{abbr}_abi", index=False)

excel_fp2 = get_excel_file_name(config, accessory=True)
with pd.ExcelWriter(excel_fp2, engine='xlsxwriter') as writer:    
    df_vx.to_excel(writer, sheet_name=f"H3_probe_matched_with_{abbr}", index=False)
    res_df_copy.to_excel(writer, sheet_name=f"{abbr}_copied_files", index=False)
    vx_nts.to_excel(writer, sheet_name=f"{abbr_lower}_gb_start_end_nts", index=False)
log.info(f"FINISNED")

