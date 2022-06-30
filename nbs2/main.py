#!/usr/bin/env python3


"""
This script finds the probes in VH and VL against the H3 probe sequence.
It then compares the sequences of the probes to the genbank reference sequence.
and calculates the number of mismatches between the two.
And finally it calculates some statistics on the mismatches.
"""

import configparser
import dataclasses
import logging
import os
import re
import shutil
import statistics
from collections import Counter
from itertools import chain
from pathlib import Path
from time import perf_counter

import pandas as pd
import utils
from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from fns import *

# START LOGGING
Path('logs').mkdir(parents=True, exist_ok=True)
log = utils.make_logger("logs/rp_log")
log.info("\nSTART")
log.info("\nCHECKING CONFIGURATION FILE CONTENTS...")

# GET PARAMS FROM CONFIG FILE
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
vh_seq_start_a, vh_seq_end_a, vl_seq_start_a, vl_seq_end_a = get_alignment_seq_start_and_end(config)

# creating dirs for copying matched abi files
res_dir_vh = checking_dirs(f"{results_dir}/{pat_vh}", log, log_msg=True, create_dir=True)
res_dir_vl = checking_dirs(f"{results_dir}/{pat_vl}", log, log_msg=True, create_dir=True)

# GET ABI FILES
log.info("\nCHECKING THE AB1 FILES...")
vh_abi_dict = {i.name.replace('.abi', '').replace(patrm_vh_abi, '') : str(i) for i in sorted([*Path(abi_sequence_folder).glob(f"*{pat_vh}*.abi")])}
vl_abi_dict = {i.name.replace('.abi', '').replace(patrm_vl_abi, '') : str(i) for i in sorted([*Path(abi_sequence_folder).glob(f"*{pat_vl}*.abi")])}

if len(vh_abi_dict) == len(vl_abi_dict):
    log.info(f"[+] INFO: There are {len(vh_abi_dict)} abi files in the abi sequence folder")
else:
    log.warning(f"[+] WARN: There is a differrence in number of abi files in vh: {len(vh_abi_dict)} and vl {len(vl_abi_dict)}")

# sample ids from abi files - combined vh and vl
sample_ids = checking_ab1_files(log, vh_abi_dict, vl_abi_dict)

log.info("\nCHECKING THE GENBANK FILES...")
vh_template_gb_dict =  {i.name.replace('.gb', '') : str(i) for i in sorted(Path(vh_template_sequence_folder).glob('*.gb'))}
vl_template_gb_dict =  {i.name.replace('.gb', '') : str(i) for i in sorted(Path(vl_template_sequence_folder).glob('*.gb'))}

if len(vh_template_gb_dict) == len(vl_template_gb_dict):
    log.info(f"[+] INFO: There are {len(vh_template_gb_dict)} genbank files in each dir vh and vl")
else:
    log.warning(f"[+] WARN: There is a differrence in number of genbank files in vh: {len(vh_template_gb_dict)} and vl {len(vl_template_gb_dict)}")

# GET PROBE SEQS INTO A DICTIOANRY
# dictionary containing probe seq and name
h3_dict = df.set_index('name', drop=True).to_dict().get('h3_nt')

# FINDING THE PROBES THAT MATCH VH AND CORRESPONDING VL FILES
# this loop  iternate through all the sample ids and find_match_on_all_h3probes function search for the probe on the normal seq and the rev complemented seq
# if there is a match which will be returened as a list

log.info(f"\nITERATING THROUGH EACH SAMPLE ID")
result_vh, result_vl = [], []
for sample in sample_ids:
    # print(f">>{sample}")
    _vhabi, _vlabi = get_abi_file_path(key=sample, vh_abi_dict=vh_abi_dict, vl_abi_dict=vl_abi_dict)
    vh_d = get_seqobj_from_abi(_vhabi)  # returns a seq record obj of VH
    vl_d = get_seqobj_from_abi(_vlabi)  # returns a seq record obj of VL
    
    # matching each probe on  VH and VL - normal and revcomp sequence
    vh_prob_search = find_match_on_all_h3probes_v2(log, h3_dict, vh_d, sample, vh_abi_dict, vl_abi_dict)
    vl_prob_search = find_match_on_all_h3probes_v2(log, h3_dict, vl_d, sample, vh_abi_dict, vl_abi_dict)
    
    if len(vh_prob_search) >=1:
        result_vh.append(vh_prob_search)
    if len(vl_prob_search) >=1:
        result_vl.append(vl_prob_search)
log.info(f"\nFINISH ITERATING THROUGH EACH SAMPLE ID")
# if vh and vl ids match then copy the files to a new place
colnames=["Match","h3_name","sample_id","vh_abi_fp","vl_abi_fp","probe_seq", "seq_record", "abi_seq", "trimmed_seq", "trimmed_seq_qual"]
df_vh = pd.DataFrame(chain.from_iterable(result_vh))
df_vl = pd.DataFrame(chain.from_iterable(result_vl))
df_vh.columns, df_vl.columns = colnames, colnames

# if the sahpe of df_vh and df_vl are same and the h3_name in both dataframes are same then copy the matched data to a new dir
if(df_vl.shape == df_vh.shape) and (df_vl.h3_name == df_vh.h3_name).all():
    log.info(f"[+] INFO: There are {df_vl.shape[0]} matches in the dataframe between vh and vl")
    log.info(f"[+] INFO: Of which {df_vh[['vh_abi_fp', 'vl_abi_fp']].drop_duplicates().shape[0]} samples has to be moved to a new dir")
else:
    log.info(f"[+] INFO: There are {df_vh.shape[0]} matches in df_vh")
    log.info(f"[+] INFO: There are {df_vl.shape[0]} matches in df_vl")

# COPY PROBE MATCHED ABI FILES TO A NEW LOC
log.info(f"\nCOPY MATCHED ABI FILES INTO NEW LOCATION")
res_df_copy = copy_mtched_abi_files_to_resdir(log, res_dir_vh, res_dir_vl, df_vh, log_msg=True)
log.info(f"\nFINISH COPYING FILES")

## PAIRWISE ALIGNMENT OF GENBANK FILES WITH MATCHED ABI FILES TO GET THE MATCHING SCORE 
# - Iterate through the VH (and VL) genbank file names [the file names are in a dict vh_template_gb_dict]
#    - if the gb file name matches with h3_name column in df_vh or df_vl 
#      - then save that into a list gb_match_with_dfvx_lst
#          - GB_Matched|gb_id|h3_name|sample_id|gb_fp|vh_abi_fp|vl_abi_fp
# - gb file name could match an h3_name; also this h3_name might be associated with different sample ids

log.info(f"\nFIND THE GB FILENASMES WITH H3 NAMES")
M_gb_abi_vh = find_gb_match_on_all_h3probes_v2(log, vh_template_gb_dict, df_vh, pattern=patrm_vh_gb, log_msg=True)
M_gb_abi_vl = find_gb_match_on_all_h3probes_v2(log, vl_template_gb_dict, df_vl, pattern=patrm_vl_gb, log_msg=True)

# combining the dataframes - not imp
all_abi_probe_match  = pd.concat([df_vh, df_vl]).reset_index(drop=True)
all_gb_match = pd.concat([M_gb_abi_vh, M_gb_abi_vl]).reset_index(drop=True)

vh_gb_abi_match_filtered = run_alignment_and_filtering(M_gb_abi_vh, df_vh, vh_template_gb_dict, vh_seq_start_a, vh_seq_end_a, par_match, par_missmatch, par_open, par_extend, par_filter_thresh, log, log_msg=True)
vl_gb_abi_match_filtered = run_alignment_and_filtering(M_gb_abi_vl, df_vl, vl_template_gb_dict, vl_seq_start_a, vl_seq_end_a, par_match, par_missmatch, par_open, par_extend, par_filter_thresh, log, log_msg=True)

# save these dataframes

if isinstance(vh_gb_abi_match_filtered, type(None)):
    vh_gb_abi_match_filtered.to_excel("VH_Match_GB_ABI.xls")
else:
    print("Check the log file")


if isinstance(vl_gb_abi_match_filtered, type(None)):
    vl_gb_abi_match_filtered.to_excel("VL_Match_GB_ABI.xls")
else:
    print("Check the log file")
