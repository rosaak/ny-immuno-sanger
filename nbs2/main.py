#!/usr/bin/env python3


"""
This script finds the probes in VH and VL against the H3 probe sequence.
It then compares the sequences of the probes to the genbank reference sequence.
and calculates the number of mismatches between the two.
And finally it calculates some statistics
- requires utils.py and fns.py
"""

import configparser
import dataclasses
import os
import re
import shutil
import statistics
from collections import Counter
from itertools import chain
from pathlib import Path
from time import perf_counter

import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
import logging
import xlsxwriter

import utils
from fns import *

#### START LOGGING
Path('logs').mkdir(parents=True, exist_ok=True)
log = utils.make_logger("logs/log")
log.info("\nSTART")
log.info("\nCHECKING CONFIGURATION FILE CONTENTS...")
#### GET PARAMS FROM CONFIG FILE
# checking config file
# if any error occurs then program terminates; then check the config file
config = configparser.ConfigParser()
config.read_file(open('config.ini'))

abi_sequence_folder = checking_dirs(config['Paths']['abi_sequence_folder'], log, log_msg=True, create_dir=False)
vh_template_sequence_folder = checking_dirs(config['Paths']['vh_template_sequence_folder'], log, log_msg=True, create_dir=False)
vl_template_sequence_folder = checking_dirs(config['Paths']['vl_template_sequence_folder'], log, log_msg=True, create_dir=False)
results_dir = checking_dirs(config['Paths']['results_dir'], log, log_msg=True, create_dir=True)
h3_nt_data_sheet_filepath = check_files(config['Files']['h3_nt_data_sheet_filepath'], log)
df = check_tsv_file(h3_nt_data_sheet_filepath, log)
excel_path_file_name = create_results_excel_file_path(config['Paths']['results_dir'], config['Files']['output_excel_file_name'])

# patterns
pat_vh, pat_vl = get_patterns(config)
# patterns to remove from abi names
patrm_vh_abi, patrm_vl_abi = get_patterns_to_rm(config)
# patterns to remove from genbank names
patrm_vh_gb, patrm_vl_gb = get_patterns_to_rm_from_genbank(config)

# get alignment parameters
par_match, par_missmatch, par_open, par_extend, par_filter_thresh = get_alignment_params(config)

# get alignment start and end regions
#vh_seq_start_a, vh_seq_end_a, vl_seq_start_a, vl_seq_end_a = get_alignment_seq_start_and_end(config)

# creating dirs for copying matched abi files
res_dir_vh = checking_dirs(f"{results_dir}/{pat_vh}", log, log_msg=True, create_dir=True)
res_dir_vl = checking_dirs(f"{results_dir}/{pat_vl}", log, log_msg=True, create_dir=True)

#### GET ABI FILES
log.info("\nCHECKING THE AB1 FILES...")
vh_abi_dict = {i.name.replace('.abi', '').replace(patrm_vh_abi, '') : str(i) for i in sorted([*Path(abi_sequence_folder).glob(f"*{pat_vh}*.abi")])}
vl_abi_dict = {i.name.replace('.abi', '').replace(patrm_vl_abi, '') : str(i) for i in sorted([*Path(abi_sequence_folder).glob(f"*{pat_vl}*.abi")])}

if len(vh_abi_dict) == len(vl_abi_dict):
    log.info(f"[+] INFO: There are {len(vh_abi_dict)} abi files in the abi sequence folder")
else:
    log.warning(f"[+] WARN: There is a differrence in number of abi files in vh: {len(vh_abi_dict)} and vl {len(vl_abi_dict)}")
    
sample_ids = checking_ab1_files(log, vh_abi_dict, vl_abi_dict)

log.info("\nCHECKING THE GENBANK FILES...")
vh_template_gb_dict =  {i.name.replace('.gb', '') : str(i) for i in sorted(Path(vh_template_sequence_folder).glob('*.gb'))}
vl_template_gb_dict =  {i.name.replace('.gb', '') : str(i) for i in sorted(Path(vl_template_sequence_folder).glob('*.gb'))}

if len(vh_template_gb_dict) == len(vl_template_gb_dict):
    log.info(f"[+] INFO: There are {len(vh_template_gb_dict)} genbank files in each dir vh and vl")
else:
    log.warning(f"[+] WARN: There is a differrence in number of genbank files in vh: {len(vh_template_gb_dict)} and vl {len(vl_template_gb_dict)}")

# get alignment start and end regions
vh_nts = extract_nts_from_start_and_end_from_genbank(vh_template_gb_dict, 6)
vl_nts = extract_nts_from_start_and_end_from_genbank(vl_template_gb_dict, 6)
vh_seq_start_a, vh_seq_end_a = vh_nts.start_nts.value_counts().nlargest(1).index[0], vh_nts.end_nts.value_counts().nlargest(1).index[0]
vl_seq_start_a, vl_seq_end_a = vl_nts.start_nts.value_counts().nlargest(1).index[0], vl_nts.end_nts.value_counts().nlargest(1).index[0]
#### GET PROBE SEQS INTO A DICTIOANRY
# dictionary containing probe seq and name
h3_dict = df.set_index('name', drop=True).to_dict().get('h3_nt')
#### FINDING THE PROBES THAT MATCH VH AND CORRESPONDING VL FILES
log.info(f"\nITERATING THROUGH EACH SAMPLE ID")
result_vh, result_vl = [], []
for sample in sample_ids:
    # print(f">>{sample}")
    _vhabi, _vlabi = get_abi_file_path(key=sample, vh_abi_dict=vh_abi_dict, vl_abi_dict=vl_abi_dict)
    _vh_d = get_seqobj_from_abi(_vhabi)  # returns a seq record obj of VH
    _vl_d = get_seqobj_from_abi(_vlabi)  # returns a seq record obj of VL
    
    # matching each probe on  VH and VL - normal and revcomp sequence
    vh_prob_search = find_match_on_all_h3probes_v3(log, h3_dict, _vh_d, _vl_d, sample, vh_abi_dict, vl_abi_dict)
    
    if len(vh_prob_search) >=1:
        result_vh.append(vh_prob_search)
log.info(f"\nFINISH ITERATING THROUGH EACH SAMPLE ID")
colnames=["Match","h3_name","sample_id","vh_abi_fp","vl_abi_fp","probe_seq", "vh_init_sr", "vl_inti_sr","vh_sr_seq_r", "vh_sr_trimmed", "vh_sr_tqlst", "vl_sr_seq_r", "vl_sr_trimmed", "vl_sr_tqlst"]
df_vh = pd.DataFrame(chain.from_iterable(result_vh))
df_vh.columns = colnames
log.info(f"[+] INFO: There are {df_vh.shape[0]} matches in df_vh")
#### COPY PROBE MATCHED ABI FILES TO A NEW LOC
log.info(f"\nCOPY MATCHED ABI FILES INTO NEW LOCATION")
res_df_copy = copy_mtched_abi_files_to_resdir(log, res_dir_vh, res_dir_vl, df_vh, log_msg=True)
log.info(f"\nFINISH COPYING FILES")
#### PAIRWISE ALIGNMENT OF GENBANK FILES WITH MATCHED ABI FILES TO GET THE MATCHING SCORE
log.info(f"\nFIND THE GB FILENASMES WITH H3 NAMES")
M_gb_abi_vh = find_gb_match_on_all_h3probes(log, vh_template_gb_dict, df_vh, pattern=patrm_vh_gb, log_msg=True)
M_gb_abi_vl = find_gb_match_on_all_h3probes(log, vl_template_gb_dict, df_vh, pattern=patrm_vl_gb, log_msg=True)
all_gb_match = pd.concat([M_gb_abi_vh, M_gb_abi_vl]).reset_index(drop=True)
log.info(f"\nALIGNEMNT BETWEEN THE MATCHED GB FILENASME AND ABI FILES")

vh_gb_abi_match_filtered = run_gb_alignment_and_filtering_v2(M_gb_abi_vh, df_vh, vh_template_gb_dict, vh_seq_start_a, vh_seq_end_a, par_match, par_missmatch, par_open, par_extend, par_filter_thresh, log, log_msg=True, is_data_vl=False)

vl_gb_abi_match_filtered = run_gb_alignment_and_filtering_v2(M_gb_abi_vl, df_vh, vl_template_gb_dict, vl_seq_start_a, vl_seq_end_a, par_match, par_missmatch, par_open, par_extend, par_filter_thresh, log, log_msg=True, is_data_vl=True)

#### MERGING DATAFRAMES
log.info(f"\nMERGING DATAFRAMES")
vh_gb_abi_match_filtered.columns = [ 'VH_'+i  if i not in ['Orient','gbid','H3_name', 'sample_id', 'GB_FP'] else i  for i in vh_gb_abi_match_filtered.columns ]
vl_gb_abi_match_filtered.columns = [ 'VL_'+i  if i not in ['Orient','gbid','H3_name', 'sample_id', 'GB_FP'] else i  for i in vl_gb_abi_match_filtered.columns ]
df_sid_abi = pd.concat(
    [pd.DataFrame([ [i, vh_abi_dict.get(i).split("/")[-1], vl_abi_dict.get(i).split("/")[-1]] for i in vh_gb_abi_match_filtered.sample_id.to_list()], 
             columns=["sample_id", "VH", "VL"]),
    pd.DataFrame([ [i, vh_abi_dict.get(i).split("/")[-1], vl_abi_dict.get(i).split("/")[-1]] for i in vl_gb_abi_match_filtered.sample_id.to_list()], 
             columns=["sample_id", "VH", "VL"])], axis=0).drop_duplicates().reset_index(drop=True)
dfx1 = pd.merge(vh_gb_abi_match_filtered,df_sid_abi,on='sample_id')
dfx2 = pd.merge(vl_gb_abi_match_filtered,df_sid_abi,on='sample_id')
dfx1_ss = dfx1[['gbid','H3_name','sample_id', 'VH_Score', 'VH_Quality_score', 'VH_Low_quality', 'VH', 'VL']]
dfx2_ss = dfx2[['gbid','H3_name','sample_id', 'VL_Score', 'VL_Quality_score', 'VL_Low_quality', 'VH', 'VL']]
dfy = pd.merge(dfx1_ss, dfx2_ss, on="sample_id", how='inner')
dfz = dfy.copy()
dfz["mean_error_prob"] = [cal_mean_error_prob(i[0], i[1]) for i in zip(dfy.VH_Low_quality.tolist(), dfy.VL_Low_quality.to_list())]
dfz = dfz[["gbid_x", "H3_name_x", "sample_id", "VH_Score", "VH_Quality_score", "VH_Low_quality", "VL_Score","VL_Quality_score", "VL_Low_quality", "mean_error_prob", "VH_x", "VL_x"]]
dfz = dfz.drop_duplicates()
dfz.columns = ['gbid', 'H3_name', 'sample_id', 'VH_Score', 'VH_Mean_Quality_score', 'VH_Low_quality', 'VL_Score', 'VL_Men_Quality_score', 'VL_Low_quality', 'Mean_Error_Prob', 'VH', 'VL']

log.info(f"[+] INFO: There are {vh_gb_abi_match_filtered.shape[0]} rows in the GB - VH aligned datafrme `vh_gb_abi_match_filtered`")
log.info(f"[+] INFO: There are {vl_gb_abi_match_filtered.shape[0]} rows in the GB - VL aligned datafrme `vl_gb_abi_match_filtered`")
log.info(f"[+] INFO: There are {dfz.shape[0]} rows in the final - filtered merged cleanedup results")

log.info(f"[+] STAT: basic stats of the numerical columns in the final dataframe : \n{dfz[dfz.dtypes[dfz.dtypes !='object'].index].describe()}")
log.info(f"\nWRTING EXCEL FILES")

#### SAVE THE MAIN RESULTS
excel_fp = get_excel_file_name(config, accessory=False)
with pd.ExcelWriter(excel_fp, engine='xlsxwriter') as writer:
    dfz.to_excel(writer, sheet_name='final_res_mean_error_prob', index=False)
    all_gb_match.to_excel(writer, sheet_name='ID matched gb and vh abi', index=False)
#### SAVE THE ACCESSORY FILES
excel_fp2 = get_excel_file_name(config, accessory=True)
with pd.ExcelWriter(excel_fp2, engine='xlsxwriter') as writer:    
    df_vh.to_excel(writer, sheet_name='H3_probe_matched_VH_and_VL', index=False)
    res_df_copy.to_excel(writer, sheet_name='VH and VL Copied Files', index=False)
    vh_nts.to_excel(writer, sheet_name='vh_gb_start_end_nts', index=False)
    vl_nts.to_excel(writer, sheet_name='vl_gb_start_end_nts', index=False)
    vh_gb_abi_match_filtered.to_excel(writer, sheet_name='VH results', index=False)
    vl_gb_abi_match_filtered.to_excel(writer, sheet_name='VL results', index=False)
log.info(f"\nFINISNED")