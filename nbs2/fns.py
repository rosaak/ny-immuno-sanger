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

# fns to check config files
    

def checking_dirs(dirname: str, log: logging, log_msg: bool=True, create_dir: bool=False):
    dirname = str(Path(dirname).resolve())
    res = utils.dir_exists_and_create(dirname, create_dir=create_dir)
    if res == 4 or res == 3:
        log.info(f"[+] INFO: {dirname} already exists") if log_msg else None
        return(dirname)
    elif res == 2:
        log.info(f"[X] WARNING: {dirname} doesn't exists, check the config file") if log_msg else None
        log.info(f"[X] ABORT") if log_msg else None
        os.abort()
    elif res == 1:
        log.info(f"[+] INFO: {dirname} created") if log_msg else None
        return(dirname)

def check_files(file_path: str, log: logging):
    file_path = str(Path(file_path).resolve())
    if utils.file_exists(file_path) == 0:
        log.info(f"[X] WARNING: {file_path} doesn't exists, check the config file")
        log.info(f"[X] ABORT")
        os.abort()
    else:
        log.info(f"[+] INFO: {file_path} exists")
        return(file_path)

def check_tsv_file(filepath: str, log: logging) -> bool:
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

def get_alignment_params(config: configparser.ConfigParser):
    par_match=config['Alignment_parameter']['match']
    par_missmatch=config['Alignment_parameter']['mismatch']
    par_open=config['Alignment_parameter']['open']
    par_extend=config['Alignment_parameter']['extend']
    par_filter_thresh=config['Alignment_parameter']['filter_thresh']
    return par_match, par_missmatch, par_open, par_extend, par_filter_thresh

def checking_ab1_files(log:logging, vh_abi_dict: dict, vl_abi_dict: dict) -> str:
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

def get_seq_from_record(dna_obj: Seq, reverse=False) -> [Seq]:
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
        seq_f = get_seq_from_record(d, reverse=False)
        seq_r = get_seq_from_record(d, reverse=True)
        vh_fn = vh_abi_dict.get(sample)
        vl_fn = vl_abi_dict.get(sample)
        if h3val in seq_f:
            log.info(f"[+] PROBE_MATCHING: Forward_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
            # print(f"Forward_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}|{h3val}|{seq_f}")
            results.append(["Forward_Strand_Match", h3key, sample, vh_fn, vl_fn, h3val, seq_f, d])
        elif h3val in seq_r:
            log.info(f"[+] PROBE_MATCHING: Reverse_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
            # print(f"Reverse_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}|{h3val}|{seq_r}")
            results.append(["Reverse_Strand_Match", h3key, sample, vh_fn, vl_fn, h3val, seq_r, d])
        else:
            pass
            # log.info(f"[-] INFO: No_match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
            # results.append(["No_match", h3key, sample, vh_fn, vl_fn, '', ''])
    return(results)


def copy_mtched_abi_files_to_resdir_x(log: logging, df_vh:pd.DataFrame, df_vl:pd.DataFrame, res_dir_vh:str, res_dir_vl:str)->pd.DataFrame:
    """
    creates a dir within the results_dir for VH and VL, with h3_name
    and copy vh and vl matched abi files into it
    """
    resx = []
    for e,i in enumerate(df_vh.h3_name.unique()[:], start=1):
        _dfvh = df_vh[df_vh.h3_name == i]
        _dfvl = df_vl[df_vl.h3_name == i]
        for j in range(_dfvh.shape[0]):
            __dfvh = _dfvh.iloc[j]
            __dfvl = _dfvl.iloc[j]
            # creating dirs 
            res_dir_vh_msid = checking_dirs(f"{res_dir_vh}/{__dfvh.h3_name}", log, log_msg=False, create_dir=True)
            res_dir_vl_msid = checking_dirs(f"{res_dir_vl}/{__dfvl.h3_name}", log, log_msg=False, create_dir=True)
            shutil.copy(__dfvh.vh_abi_fp, res_dir_vh_msid)
            shutil.copy(__dfvl.vl_abi_fp, res_dir_vl_msid)
            log.info(f"[+] COPY: {__dfvh.vh_abi_fp} to {res_dir_vh_msid}")
            log.info(f"[+] COPY: {__dfvh.vl_abi_fp} to {res_dir_vl_msid}")
            resx.append([__dfvh.h3_name, __dfvh.sample_id, f"{res_dir_vh}/{__dfvh.sample_id}", __dfvh.vh_abi_fp])
            resx.append([__dfvl.h3_name, __dfvl.sample_id, f"{res_dir_vl}/{__dfvl.sample_id}", __dfvl.vl_abi_fp])
    return pd.DataFrame(resx, columns=["h3_name", "sample_id", "abi_out_loc", "abi_initial_filepath"])


def copy_mtched_abi_files_to_resdir(res_dir_vh:str, res_dir_vl:str, vx_abi_fp:str, log: logging, log_msg:bool= True, df_vh:pd.DataFrame = None)->pd.DataFrame:
    """
    creates a dir within the results_dir for VH and VL, with h3_name
    and copy matched vh and vl abi files into it
    """
    def _cp_files(log: logging, df:pd.DataFrame, res_dir_vx:str, vx_abi_fp:str):
        res = []
        for i in df.h3_name.unique()[:]:
            _df = df[df.h3_name == i]
            for j in range(_df.shape[0]):
                __df = _df.iloc[j]
                res_dir_msid = checking_dirs(f"{res_dir_vx}/{__df.h3_name}", log, log_msg=False, create_dir=True)
                shutil.copy(__df[vx_abi_fp], res_dir_msid)
                log.info(f"[+] COPY: {__df[vx_abi_fp]} to {res_dir_msid}") if log_msg else None
                res.append([__df.h3_name, __df.sample_id, f"{res_dir_vx}/{__df.sample_id}", __df[vx_abi_fp]])
        return pd.DataFrame(res, columns=["h3_name", "sample_id", "abi_out_loc", "abi_initial_filepath"])
    
    df1 = _cp_files(log, df_vh, res_dir_vh, 'vh_abi_fp')
    df2 = _cp_files(log, df_vh, res_dir_vl, 'vl_abi_fp')
    return pd.concat([df1, df2])
        
    
    # resx, resy = [],[]
    # for i df_vh.h3_name.unique()[:]:
    #     _dfvh = df_vh[df_vh.h3_name == i]
    #     for j in range(_dfvh.shape[0]):
    #         __dfvh = _dfvh.iloc[j]
    #         res_dir_vh_msid = checking_dirs(f"{res_dir_vh}/{__dfvh.h3_name}", log, log_msg=False, create_dir=True)
    #         shutil.copy(__dfvh.vh_abi_fp, res_dir_vh_msid)
    #         log.info(f"[+] INFO: copying {__dfvh.vh_abi_fp} to {res_dir_vh_msid}")
    #         resx.append([__dfvh.h3_name, __dfvh.sample_id, f"{res_dir_vh}/{__dfvh.sample_id}", __dfvh.vh_abi_fp])
    # for i in df_vl.h3_name.unique()[:]:
    #     _dfvl = df_vl[df_vl.h3_name == i]
    #     for j in range(_dfvl.shape[0]):
    #         __dfvl = _dfvl.iloc[j]
    #         res_dir_vl_msid = checking_dirs(f"{res_dir_vl}/{__dfvl.h3_name}", log, log_msg=False, create_dir=True)
    #         shutil.copy(__dfvl.vl_abi_fp, res_dir_vl_msid)
    #         log.info(f"[+] INFO: copying {__dfvh.vl_abi_fp} to {res_dir_vl_msid}")
    #         resy.append([__dfvl.h3_name, __dfvl.sample_id, f"{res_dir_vl}/{__dfvl.sample_id}", __dfvl.vl_abi_fp])
    # return pd.DataFrame(resx, columns=["h3_name", "sample_id", "abi_out_loc", "abi_initial_filepath"])
    
    
def main():
    pass

if __name__ == '__main__':
    main()