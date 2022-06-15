#!/usr/bin/env python3


"""
This script requires finds the probes in VH and VL against the H3 reference.
It then compares the sequences of the probes to the H3 reference.
and calculates the number of mismatches between the two.
And finally it calculates some statistics on the mismatches.
"""

from pydoc import doc
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
from IPython.display import display
import logging

import utils


# Functions

def checking_dirs(dirname: str, log: logging.Logger, create_dir: bool=False):
    dirname = str(Path(dirname).resolve())
    res = utils.dir_exists_and_create(dirname, create_dir=create_dir)
    if res == 4 or res == 3:
        log.info(f"[+] INFO: {dirname} already exists")
        return(dirname)
    elif res == 2:
        log.info(f"[X] WARNING : {dirname} doesn't exists, check the config file")
        log.info(f"[X] ABORT")
        os.abort()
    elif res == 1:
        log.info(f"[+] INFO: {dirname} created")
        return(dirname)
    
def check_files(file_path: str, log: logging.Logger):
    file_path = str(Path(file_path).resolve())
    if utils.file_exists(file_path) == 0:
        log.info(f"[X] WARNING : {file_path} doesn't exists, check the config file")
        log.info(f"[X] ABORT")
        os.abort()
    else:
        log.info(f"[+] INFO: {file_path} exists")
        return(file_path)
    
def check_tsv_file(filepath: str, log: logging.Logger) -> bool:
    """
    Returns True if the tsv file had name and h3_nt columns and no missing values in it
    """
    df = pd.read_csv(filepath, sep="\t")
    if [i for i, j in zip(sorted(df.columns.tolist(), reverse=True), sorted(['name', 'h3_nt'], reverse=True)) if i == j] == ['name', 'h3_nt']:
        log.info(f"[+] INFO: {filepath} contains column names 'name' and 'h3_nt'")
        return df[['name', 'h3_nt']].dropna()
    else:
        log.info(f"[+] WARNING: {filepath} should contain column names 'name' and 'h3_nt'")
        log.info(f"[X] ABORT")
        os.abort()
        
def create_results_excel_file_path(res_fp:str, res_fn: str):
    return str(Path(os.path.join(res_fp, res_fn)).resolve())
    
def get_patterns(config: configparser.ConfigParser):
    pat_vh,pat_vl=config['Pattern']['pat_vh'], config['Pattern']['pat_vl']
    return pat_vh, pat_vl

def get_patterns_to_rm(config: configparser.ConfigParser):
    patrm_vh, patrm_vl = config['Pattern']['vh_pat_to_rm'], config['Pattern']['vl_pat_to_rm']
    return patrm_vh, patrm_vl

def checking_ab1_files(vh_abi_dict: dict, vl_abi_dict: dict) -> str:
    """
    Return sample ids
    """
    vh_keys = [i for i in [*vh_abi_dict.keys()]]
    vl_keys = [i for i in [*vl_abi_dict.keys()]]
    if len(vh_abi_dict) == len(vl_abi_dict):
        log.info(f"[+] INFO: Same number of samples (n={len(vh_abi_dict)}) in vh_abi_dict and vl_abi_dict")
        if vh_keys == vl_keys:
            log.info(f"[+] INFO: The sample names in vh_abi_dict abd vl_abi_dict are same")
            return vl_keys
        else:
            log.info(f"[+] WARNING: Some sample names in vh_abi_dict abd vl_abi_dict are different; check the patterns given in `vh_pat_to_rm` or `vl_pat_to_rm`")
            log.info(f"[X] ABORT")
            os.abort()
        
def get_abi_file_path(key:str, vh_abi_dict: dict, vl_abi_dict: dict) -> str:
    vh = vh_abi_dict.get(key)
    vl = vl_abi_dict.get(key)
    return vh, vl

def _abi_trim(seq_record):
    start = False   # flag for starting position of trimmed sequence
    segment = 20    # minimum sequence length
    trim_start = 0  # init start index
    cutoff = 0.01   # default cutoff value for calculating base score
    # calculate base score
    score_list = [cutoff - (10 ** (qual / -25.0)) for qual in seq_record.letter_annotations["phred_quality"]]
    # calculate cummulative score
    # if cummulative value < 0, set it to 0
    # first value is set to 0, because of the assumption that
    # the first base will always be trimmed out
    cummul_score = [0]
    for i in range(1, len(score_list)):
        score = cummul_score[-1] + score_list[i]
        if score < 0:
            cummul_score.append(0)
        else:
            cummul_score.append(score)
            if not start:
                    # trim_start = value when cummulative score is first > 0
                trim_start = i
                start = True
        # trim_finish = index of highest cummulative score,
        # marking the end of sequence segment with highest cummulative score
    trim_finish = cummul_score.index(max(cummul_score))
    return seq_record[trim_start:trim_finish]

def get_seqobj_from_abi(abi_fp: str) -> [Seq, Seq, list]:
    dna = SeqIO.read(abi_fp, 'abi')
    dna_trimmed = _abi_trim(dna)
    quality_trimmed = dna_trimmed.letter_annotations["phred_quality"]
    return dna, dna_trimmed, quality_trimmed

def get_seq_from_recorf(dna_obj: Seq, reverse=False) -> [Seq]:
    if not reverse:
        return dna_obj.seq
    else:
        return dna_obj.reverse_complement().seq

def find_match_on_all_h3probes(log: logging, h3_dict: dict, d: SeqIO.SeqRecord, sample:str, vh_abi_dict: dict, vl_abi_dict: dict) -> list:
    """
    returns a list of matched, h3_probe_name and h3_probe_seq, vh,vl abi file path etc.
    """
    results = []
    for h3key, h3val in h3_dict.items():
        seq_f = get_seq_from_recorf(d,reverse=False)
        seq_r = get_seq_from_recorf(d,reverse=True)
        vh_fn = vh_abi_dict.get(sample)
        vl_fn = vl_abi_dict.get(sample)
        if h3val in seq_f:
            log.info(f"[+] INFO: Forward_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
            # print(f"Forward_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}|{h3val}|{seq_f}")
            results.append(["Forward_Strand_Match", h3key, sample, vh_fn, vl_fn, h3val, seq_f])
        elif h3val in seq_r:
            log.info(f"[+] INFO: Reverse_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
            # print(f"Reverse_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}|{h3val}|{seq_r}")
            results.append(["Reverse_Strand_Match", h3key, sample, vh_fn, vl_fn, h3val, seq_f])
        else:
            pass
            # log.info(f"[-] INFO: No_match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
            # results.append(["No_match", h3key, sample, vh_fn, vl_fn, '', ''])
    return(results)


# Init Logging
Path('logs').mkdir(parents=True, exist_ok=True)
log = utils.make_logger("logs/rp_log")

config = configparser.ConfigParser()
config.read_file(open('config.ini'))

log.info("START")
log.info("\nCHECKING CONFIGURATION FILE CONTENTS...")
# checking config file
# if any error occurs then program terminates; then check the config file

abi_sequence_folder = checking_dirs(config['Paths']['abi_sequence_folder'], log, create_dir=False)
vh_template_sequence_folder = checking_dirs(config['Paths']['vh_template_sequence_folder'], log, create_dir=False)
vl_template_sequence_folder = checking_dirs(config['Paths']['vl_template_sequence_folder'], log, create_dir=False)
results_dir = checking_dirs(config['Paths']['results_dir'], log, create_dir=True)
h3_nt_data_sheet_filepath = check_files(config['Files']['h3_nt_data_sheet_filepath'], log)
df = check_tsv_file(h3_nt_data_sheet_filepath, log)
excel_path_file_name = create_results_excel_file_path(config['Paths']['results_dir'], config['Files']['output_excel_file_name'])

log.info("\nCHECKING THE AB1 FILES...")
pat_vh, pat_vl = get_patterns(config)
patrm_vh, patrm_vl = get_patterns_to_rm(config)

# abi files and names to dictionary
vh_abi_dict = {i.name.replace('.abi', '').replace(patrm_vh, '') : str(i) for i in sorted([*Path(abi_sequence_folder).glob(f"*{pat_vh}*.abi")])}
vl_abi_dict = {i.name.replace('.abi', '').replace(patrm_vl, '') : str(i) for i in sorted([*Path(abi_sequence_folder).glob(f"*{pat_vl}*.abi")])}

# checking if the abi files are present in the abi folder - returns a list
sample_ids = checking_ab1_files(vh_abi_dict, vl_abi_dict)

# dictionary containing probe seq and name
h3_dict = df.set_index('name', drop=True).to_dict().get('h3_nt')

# LOOP THROUGH THE SAMPLE IDs
# this loop  iternate through all the sample ids and find_match_on_all_h3probes function search for the probe on the normal seq and the rev complemented seq
# if there is a match which will be returened as a list

log.info(f"\nITERATING THROUGH EACH SAMPLE ID")
result_vh, result_vl = [], []
for sample in sample_ids:
    # print(f">>{sample}")
    _vhabi, _vlabi = get_abi_file_path(key=sample, vh_abi_dict=vh_abi_dict, vl_abi_dict=vl_abi_dict)
    vh_d, vh_dtrim, vh_dqtrimlst = get_seqobj_from_abi(_vhabi)
    vl_d, vl_dtrim, vl_dqtrimlst = get_seqobj_from_abi(_vhabi)
    
    # matching each probe on  VH and VL - normal and revcomp sequence
    vh_prob_search = find_match_on_all_h3probes(log, h3_dict, vh_d, sample, vh_abi_dict, vl_abi_dict)
    vl_prob_search = find_match_on_all_h3probes(log, h3_dict, vl_d, sample, vh_abi_dict, vl_abi_dict)
    if len(vh_prob_search) >=1:
        result_vh.append(vh_prob_search)
    if len(vl_prob_search) >=1:
        result_vl.append(vl_prob_search)
log.info(f"\nFINISH ITERATING THROUGH EACH SAMPLE ID")
# if vh and vl ids match then copy the files to a new place
colnames=["Match","h3_name","sample_id","vh_abi_fp","vl_abi_fp","probe_seq", "abi_seq"]
df_vh = pd.DataFrame(chain.from_iterable(result_vh))
df_vl = pd.DataFrame(chain.from_iterable(result_vl))
df_vh.columns, df_vl.columns = colnames, colnames

# if the sahpe of df_vh and df_vl are same and the h3_name in both dataframes are same then copy the matched data to a new dir
if( df_vl.shape == df_vh.shape) and (df_vl.h3_name == df_vh.h3_name).all():
        log.info(f"[+] INFO: There are {df_vl.shape[0]} matches in the dataframe between vh and vl")
        log.info(f"[+] INFO: Of which {df_vh[['vh_abi_fp', 'vl_abi_fp']].drop_duplicates().shape}[0] samples has to be moved to a new dir")

