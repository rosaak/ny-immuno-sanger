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

def get_patterns_to_rm_from_genbank(config: configparser.ConfigParser):
    patrm_vh_from_gbname, patrm_vl_from_gbname = config['Pattern']['vh_pat_to_rm_from_gb_name'], config['Pattern']['vl_pat_to_rm_from_gb_name']
    return patrm_vh_from_gbname, patrm_vl_from_gbname

def get_alignment_params(config: configparser.ConfigParser):
    par_match=config['Alignment_parameter']['match']
    par_missmatch=config['Alignment_parameter']['mismatch']
    par_open=config['Alignment_parameter']['open']
    par_extend=config['Alignment_parameter']['extend']
    par_filter_thresh=config['Alignment_parameter']['filter_thresh']
    return int(par_match), int(par_missmatch), int(par_open), int(par_extend), int(par_filter_thresh)

def get_alignment_seq_start_and_end(config: configparser.ConfigParser):
    vh_seq_aln_start_a=config['Alignment_parameter']['vh_seq_aln_start_a']
    vh_seq_aln_end_a=config['Alignment_parameter']['vh_seq_aln_end_a']
    vl_seq_aln_start_a=config['Alignment_parameter']['vl_seq_aln_start_a']
    vl_seq_aln_end_a=config['Alignment_parameter']['vl_seq_aln_end_a']
    return vh_seq_aln_start_a, vh_seq_aln_end_a, vl_seq_aln_start_a, vl_seq_aln_end_a


def get_excel_file_name(config: configparser.ConfigParser,  accessory:bool=False)-> str:
    """
    Create the excel file name from config entry. datetime will be appened to the filename
    Input:
    config: config parser
    aceesory: A suffix of accissory to add to the excel file name
    """
    if not accessory:
        return os.path.join(config["Paths"]["results_dir"], config["Files"]["output_excel_file_name"] +"_"+str(utils.now())+".xlsx")
    else:
        return os.path.join(config["Paths"]["results_dir"], config["Files"]["output_excel_file_name"] +"_accessory_"+str(utils.now())+".xlsx")


def extract_nts_from_start_and_end_from_genbank(gb_dict:dict, nt_len:int=6) -> pd.DataFrame:
    
    def _get_nts(fp:str, nt_len:int=6):
        _ = SeqIO.read(fp, 'gb')
        return f"{fp},{_.name},{_.seq[0:nt_len]},{_.seq[-nt_len:]}".split(",")
    
    res = pd.DataFrame([_get_nts(i[1]) for i in gb_dict.items()])
    res.columns = ["gb_fp", "gb_name", "start_nts", "end_nts"]
    return res

def checking_ab1_files(log:logging, vh_abi_dict: dict, vl_abi_dict: dict) -> list[str]:
    """
    Return sample ids list if vh and vl are same
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

def get_seqobj_from_abi(abi_fp: str):
    dna = SeqIO.read(abi_fp, 'abi')
    return dna

def get_trimmed_seq_record(seq_record: SeqIO.SeqRecord, get_quality: bool=False):
    dna_trimmed = _abi_trim(seq_record)
    quality_trimmed = dna_trimmed.letter_annotations["phred_quality"]
    if get_quality:
        return dna_trimmed, quality_trimmed
    else:
        return dna_trimmed

def get_seq_from_record(dna_obj: Seq, reverse=False) -> [Seq]:
    if not reverse:
        return dna_obj.seq
    else:
        return dna_obj.reverse_complement().seq


# def find_match_on_all_h3probes(log: logging, h3_dict: dict,  d: SeqIO.SeqRecord, sample:str, vh_abi_dict: dict, vl_abi_dict: dict) -> list:
#     """
#     Returns: [orientation of match, h3_probe_name, sample_id, vh_abi_dp, vl_abi_fp, probe_seq, abi_seq, sequence_record]
#     inputs
#     log: logger 
#     h3_dict: dictionary containing h3_name and probe_name
#     d: Sequence Record of 
#     sample:
#     vh_abi_dict:
#     vl_abi_dict:
#     """
#     results = []
#     for h3key, h3val in h3_dict.items():
#         d_trimmed, d_trimquallst = get_trimmed_seq_record(d, get_quality=True)
#         seq_f = get_seq_from_record(d, reverse=False)
#         seq_r = get_seq_from_record(d, reverse=True)
#         vh_fn = vh_abi_dict.get(sample)
#         vl_fn = vl_abi_dict.get(sample)
#         if h3val in seq_f:

#             log.info(f"[+] PROBE_MATCHING: Forward_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
#             # print(f"Forward_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}|{h3val}|{seq_f}")
#             results.append(["Forward_Strand_Match", h3key, sample, vh_fn, vl_fn, h3val, seq_f, d])
#         elif h3val in seq_r:
#             log.info(f"[+] PROBE_MATCHING: Reverse_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
#             # print(f"Reverse_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}|{h3val}|{seq_r}")
#             results.append(["Reverse_Strand_Match", h3key, sample, vh_fn, vl_fn, h3val, seq_r, d])
#         else:
#             pass
#             # log.info(f"[-] INFO: No_match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
#             # results.append(["No_match", h3key, sample, vh_fn, vl_fn, '', ''])
#     return(results)


# def find_match_on_all_h3probes_v2(log: logging, h3_dict: dict,  d: SeqIO.SeqRecord, sample:str, vh_abi_dict: dict, vl_abi_dict: dict) -> list:
#     """
#     Returns: [orientation of match, h3_probe_name, sample_id, vh_abi_dp, vl_abi_fp, probe_seq, seq_record, abi_seq, trimmed_seq, trimmed_seq_qual]
#     inputs
#     log: logger 
#     h3_dict: dictionary containing h3_name and probe_name
#     d: Sequence Record from the abi file
#     sample: sample id - comes from sample_ids
#     vh_abi_dict: dictionary containing vh names and full file path
#     vl_abi_dict: dictionary containing vl names and full file path
#     """
#     results = []
#     for h3key, h3val in h3_dict.items():
#         d_trimmed, d_trimquallst = get_trimmed_seq_record(d, get_quality=True)
#         seq_f = get_seq_from_record(d, reverse=False)
#         seq_r = get_seq_from_record(d, reverse=True)
#         vh_fn = vh_abi_dict.get(sample)
#         vl_fn = vl_abi_dict.get(sample)
#         if h3val in seq_f:
#             log.info(f"[+] PROBE_MATCHING: Forward_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
#             results.append(["Forward_Strand_Match", h3key, sample, vh_fn, vl_fn, h3val, d, seq_f, d_trimmed, d_trimquallst])
#         elif h3val in seq_r:
#             log.info(f"[+] PROBE_MATCHING: Reverse_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
#             results.append(["Reverse_Strand_Match", h3key, sample, vh_fn, vl_fn, h3val, d, seq_r, d_trimmed, d_trimquallst])
#         else:
#             pass
#     return(results)



def find_match_on_all_h3probes_v3(log: logging, h3_dict: dict,  vh_sr: SeqIO.SeqRecord, vl_sr: SeqIO.SeqRecord, sample:str, vh_abi_dict: dict, vl_abi_dict: dict) -> list:
    """
    If vh matches then get the corresponding vl also
    Returns: [orientation of match, h3_probe_name, sample_id, vh_abi_dp, vl_abi_fp, probe_seq, vh_seq_record, vl_seq_record, vh_trimmed_seq, vh_trimmed_seq_qual, vl_trimmed_seq, vl_trimmed_seq_qual]
    inputs
    log: logger 
    h3_dict: dictionary containing h3_name and probe_name
    d: Sequence Record from the abi file
    sample: sample id - comes from sample_ids
    vh_abi_dict: dictionary containing vh names and full file path
    vl_abi_dict: dictionary containing vl names and full file path
    """
    results = []
    for h3key, h3val in h3_dict.items():
        vh_sr_trimmed, vh_sr_trimquallst = get_trimmed_seq_record(vh_sr, get_quality=True)
        vl_sr_trimmed, vl_sr_trimquallst = get_trimmed_seq_record(vl_sr, get_quality=True)
        
        vh_sr_seq_f = get_seq_from_record(vh_sr, reverse=False)
        vl_sr_seq_f = get_seq_from_record(vl_sr, reverse=False)
        vh_sr_seq_r = get_seq_from_record(vh_sr, reverse=True)
        vl_sr_seq_r = get_seq_from_record(vl_sr, reverse=True)

        vh_fn = vh_abi_dict.get(sample)
        vl_fn = vl_abi_dict.get(sample)
        
        if h3val in vh_sr_seq_f:
            log.info(f"[+] PROBE_MATCHING: Forward_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
            results.append(["Forward_Strand_Match", h3key, sample, vh_fn, vl_fn, h3val, vh_sr, vl_sr, vh_sr_seq_f, vh_sr_trimmed, vh_sr_trimquallst, vl_sr_seq_f, vl_sr_trimmed, vl_sr_trimquallst])
        elif h3val in vh_sr_seq_r:
            log.info(f"[+] PROBE_MATCHING: Reverse_Strand_Match|{h3key}|{sample}|{vh_fn}|{vl_fn}")
            results.append(["Reverse_Strand_Match", h3key, sample, vh_fn, vl_fn, h3val, vh_sr, vl_sr, vh_sr_seq_r, vh_sr_trimmed, vh_sr_trimquallst,  vl_sr_seq_r, vl_sr_trimmed, vl_sr_trimquallst])
                           
        else:
            pass
    return(results)



# def copy_mtched_abi_files_to_resdir_x(log: logging, df_vh:pd.DataFrame, df_vl:pd.DataFrame, res_dir_vh:str, res_dir_vl:str)->pd.DataFrame:
#     """
#     creates a dir within the results_dir for VH and VL, with h3_name
#     and copy vh and vl matched abi files into it
#     """
#     resx = []
#     for e,i in enumerate(df_vh.h3_name.unique()[:], start=1):
#         _dfvh = df_vh[df_vh.h3_name == i]
#         _dfvl = df_vl[df_vl.h3_name == i]
#         for j in range(_dfvh.shape[0]):
#             __dfvh = _dfvh.iloc[j]
#             __dfvl = _dfvl.iloc[j]
#             # creating dirs 
#             res_dir_vh_msid = checking_dirs(f"{res_dir_vh}/{__dfvh.h3_name}", log, log_msg=False, create_dir=True)
#             res_dir_vl_msid = checking_dirs(f"{res_dir_vl}/{__dfvl.h3_name}", log, log_msg=False, create_dir=True)
#             shutil.copy(__dfvh.vh_abi_fp, res_dir_vh_msid)
#             shutil.copy(__dfvl.vl_abi_fp, res_dir_vl_msid)
#             log.info(f"[+] COPY: {__dfvh.vh_abi_fp} to {res_dir_vh_msid}")
#             log.info(f"[+] COPY: {__dfvh.vl_abi_fp} to {res_dir_vl_msid}")
#             resx.append([__dfvh.h3_name, __dfvh.sample_id, f"{res_dir_vh}/{__dfvh.sample_id}", __dfvh.vh_abi_fp])
#             resx.append([__dfvl.h3_name, __dfvl.sample_id, f"{res_dir_vl}/{__dfvl.sample_id}", __dfvl.vl_abi_fp])
#     return pd.DataFrame(resx, columns=["h3_name", "sample_id", "abi_out_loc", "abi_initial_filepath"])


def copy_mtched_abi_files_to_resdir(log: logging, res_dir_vh: str, res_dir_vl:str, df_vh:pd.DataFrame = None,  log_msg:bool= True)->pd.DataFrame:
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
    return pd.concat([df1, df2]).reset_index(drop=True)
        
    
def find_gb_match_on_all_h3probes(log: logging, gb_tmplate_dict: dict, df_vx: pd.DataFrame, pattern='VH-', log_msg: bool=True) -> pd.DataFrame:
    gb_match_with_dfvx_lst = []
    for gb_id in [*gb_tmplate_dict.keys()]:
        _ = df_vx[df_vx.h3_name == gb_id.replace(pattern, '')]
        if not (_).empty:
            for j in range(_.shape[0]):
                log.info(f"[+] GENBANK_MATCHING: {gb_id}|{_.iloc[j].h3_name}|{_.iloc[j].sample_id}|{gb_tmplate_dict.get(gb_id)}|{_.iloc[j].vh_abi_fp}|{_.iloc[j].vl_abi_fp}") if log_msg else None
                gb_match_with_dfvx_lst.append(['GB_Matched',gb_id, _.iloc[j].h3_name, _.iloc[j].sample_id, gb_tmplate_dict.get(gb_id), _.iloc[j].vh_abi_fp, _.iloc[j].vl_abi_fp])
        else:
            pass
    _df = pd.DataFrame(gb_match_with_dfvx_lst)
    _df.columns = ["Match", "gb_id", "h3_name", "sample_id", "gb_fp", "vh_abi_fp", "vl_abi_fp"]
    return _df
    
    
def align(A_seq: SeqIO.SeqRecord, B_seq: SeqIO.SeqRecord, match :int, mismatch: int, gap_open: int, gap_extend: int):
    return pairwise2.align.globalms(sequenceA=A_seq, sequenceB=B_seq, match=match, mismatch=mismatch, open=gap_open, extend=gap_extend, penalize_end_gaps = False)
    
    
def run_alignment_and_filtering(M_gb_abi_vx: pd.DataFrame, 
                                df_vx:pd.DataFrame, 
                                vx_template_gb_dict: dict, 
                                vx_seq_start_a:str, 
                                vx_seq_end_a:str,
                                par_match:int,
                                par_missmatch:int,
                                par_open:int,
                                par_extend:int,
                                par_filter_thresh:int,
                                log :logging, 
                                log_msg: bool=False):
    """
    Returns a dataframe containing gbid,H3_name,GB_FP,ABI_FP,Score,Quality_score,Low_quality
    when for all the datasets where score is above 700
    """
    filtered_res = []
    log.info(f"[+] FILTERED_MATCH: outer_index|inner_index|GB_ID_MATCHED|H3_name|GB_FP|ABI_VH_FP|Score|Quality_score|Low_quality")if log_msg else None
    for e,i in enumerate(range(M_gb_abi_vx.shape[0])[:]):
        gbid = M_gb_abi_vx.iloc[i].gb_id
        gbfp = M_gb_abi_vx.iloc[i].gb_fp
        _dfvx_ss = df_vx[df_vx.h3_name == M_gb_abi_vx.iloc[i].h3_name]   
        for e2, abx in enumerate(range(_dfvx_ss.shape[0])):
            # match between gb and vx-abi
            try:
                seqA = SeqIO.read(vx_template_gb_dict.get(gbid), 'gb')
                seqB = _dfvx_ss.iloc[abx].abi_seq
                seqBinit = _dfvx_ss.iloc[abx].seq_record.seq
                quality = _dfvx_ss.iloc[abx].trimmed_seq_qual
                for a in align(seqA.seq, seqB, match=par_match, mismatch=par_missmatch, gap_open=par_open, gap_extend=par_extend):
                    score = int(a[2])
                    if score >= int(par_filter_thresh):

                        # getting the aligned seqA and seqB
                        aligned_seq_a = a[0]
                        aligned_seq_a = aligned_seq_a.replace('-', '')
                        aligned_seq_a_start = aligned_seq_a[0:6]
                        aligned_seq_a_length = len(aligned_seq_a) - 6
                        aligned_seq_a_end = aligned_seq_a[aligned_seq_a_length:len(aligned_seq_a)]
                        aligned_seq_b = a[1]
                        # futher filtering of seqB
                        if '-' not in aligned_seq_b:
                            if vh_seq_start_a in aligned_seq_a_start:
                                if vh_seq_end_a in aligned_seq_a_end: 
                                    aligned_seq = a[0]
                                    aligned_seq = aligned_seq.replace('-', '')
                                    aligned_seq = Seq(aligned_seq)
                                    aligned_seq_rev = aligned_seq.reverse_complement()
                                    dna_string = str(seqBinit)
                                    aligned_seq_rev = str(aligned_seq_rev)
                                    start = re.search(aligned_seq_rev, dna_string).start()
                                    end = re.search(aligned_seq_rev, dna_string).end()
                                    quality = quality[start:end]
                                    quality_score = statistics.mean(quality)
                                    lowest_quality = min(quality)
                                    # print(f"{e} | {abx} | {gbid} | {_dfvx_ss.iloc[abx].h3_name} | {gbfp} | {_dfvx_ss.iloc[abx].vh_abi_fp} | {score}| {quality_score} | {lowest_quality}")
                                    log.info(f"[+] FILTERED_GOODMATCH: {e}|{abx}|{gbid}|{_dfvx_ss.iloc[abx].sample_id}|{_dfvx_ss.iloc[abx].h3_name}|{gbfp}|{_dfvx_ss.iloc[abx].vh_abi_fp}|{score}|{quality_score}|{lowest_quality}")if log_msg else None
                                    filtered_res.append([gbid, _dfvx_ss.iloc[abx].sample_id, _dfvx_ss.iloc[abx].h3_name, gbfp, _dfvx_ss.iloc[abx].vh_abi_fp, score, quality_score, lowest_quality])
                else:
                    log.info(f"[+] FILTERED_BADMATCH: {e}|{abx}|{gbid}|{_dfvx_ss.iloc[abx].sample_id}|{_dfvx_ss.iloc[abx].h3_name}|{gbfp}|{_dfvx_ss.iloc[abx].vh_abi_fp}|{score}|---|---")if log_msg else None
            except:
                pass
    if not len(filtered_res) == 0:
        dfres = pd.DataFrame(filtered_res)
        dfres.columns = ["gbid", "sample_id", "H3_name", "GB_FP", "ABI_FP", "Score", "Mean_Quality_score", "Low_quality"]
        return dfres
    else:
        return None
    
# def run_gb_alignment_and_filtering_v1(M_gb_abi_vx, df_vx, vx_template_gb_dict, vx_seq_start_a, vx_seq_end_a, par_match, par_missmatch, par_open, par_extend, par_filter_thresh, log, log_msg=True, is_data_vl=False):
#     filtered_res = []
#     for e,i in enumerate(range(M_gb_abi_vx.shape[0])[0:]):
#         gbid = M_gb_abi_vx.iloc[i].gb_id
#         gbfp = M_gb_abi_vx.iloc[i].gb_fp
#         _dfvx_ss = df_vx[df_vx.h3_name == M_gb_abi_vx.iloc[i].h3_name]
#         for e2, abx in enumerate(range(_dfvx_ss.shape[0])):
#             # genbank
#             seqA = SeqIO.read(vx_template_gb_dict.get(gbid), 'gb')
#             # vl
#             flag = ''
#             if is_data_vl:
#                 flag = 'VL'
#                 seqB = get_seqobj_from_abi(_dfvx_ss.iloc[abx].vl_abi_fp)
#                 seqB_vxinit = _dfvx_ss.iloc[abx]["vl_inti_sr"].seq
#                 seqB_vx, quality = get_trimmed_seq_record(seqB, get_quality=True)
#                 seqB_vx = get_seq_from_record(seqB_vx, reverse=True)
#             else:
#                 flag = 'VH'
#                 # seqB = get_seqobj_from_abi(_dfvx_ss.iloc[abx].vh_abi_fp)
#                 # # seqB = _dfvx_ss.iloc[abx].vh_sr_seq_r
#                 # seqB_vxinit = _dfvx_ss.iloc[abx].vh_init_sr.seq
#                 # quality = _dfvx_ss.iloc[abx].vh_sr_tqlst
                
#                 seqB_vx = _dfvx_ss.iloc[abx].vh_sr_seq_r
#                 seqB_vxinit = _dfvx_ss.iloc[abx].vh_init_sr.seq
#                 quality = _dfvx_ss.iloc[abx].vh_sr_tqlst
                
#             # match between gb and vl-abi
#             try:          
#                 for a in align(seqA.seq, seqB_vx, match=par_match, mismatch=par_missmatch, gap_open=par_open, gap_extend=par_extend):  #😄
#                     score = int(a[2])
#                     # print(f"{e} | {abx} | {flag} | {_dfvx_ss.iloc[abx].vl_abi_fp} |{score}")
#                     if score >= int(par_filter_thresh):
#                         # getting the aligned seqA and seqB
#                         aligned_seq_a = a[0]
#                         aligned_seq_a = aligned_seq_a.replace('-', '')
#                         aligned_seq_a_start = aligned_seq_a[0:6]
#                         aligned_seq_a_length = len(aligned_seq_a) - 6
#                         aligned_seq_a_end = aligned_seq_a[aligned_seq_a_length:len(aligned_seq_a)]
#                         aligned_seq_b = a[1]
#                         #print(f"{e} | {abx} | VL | {_dfvx_ss.iloc[abx].vl_abi_fp} |{score} | {aligned_seq_a_start} | {aligned_seq_a_end} | {len(aligned_seq_a)} | {len(aligned_seq_b)} | {aligned_seq_a[0:6] +'..'+ aligned_seq_a[-6:]} | {aligned_seq_b[0:6] +'..'+ aligned_seq_b[-6:]}")
#                         # futher filtering of seqB
#                         if '-' not in aligned_seq_b:
#                             if vx_seq_start_a in aligned_seq_a_start:    # 😄
#                                 if vx_seq_end_a in aligned_seq_a_end:    # 😄
#                                     aligned_seq = a[0]
#                                     aligned_seq = aligned_seq.replace('-', '')
#                                     aligned_seq = Seq(aligned_seq)
#                                     aligned_seq_rev = aligned_seq.reverse_complement()
#                                     dna_string = str(seqB_vxinit)        #😄
#                                     aligned_seq_rev = str(aligned_seq_rev)
#                                     start = re.search(aligned_seq_rev, dna_string).start()
#                                     end = re.search(aligned_seq_rev, dna_string).end()
#                                     quality = quality[start:end]   #😄
#                                     quality_score = statistics.mean(quality) #😄
#                                     lowest_quality = min(quality) #😄
#                                     # print(f"{e} | {abx} | {gbid} | {_dfvx_ss.iloc[abx].h3_name} | {gbfp} | {_dfvx_ss.iloc[abx].vh_abi_fp} | {score}| {quality_score} | {lowest_quality}")
#                                     log.info(f"[+] FILTERED_GOODMATCH: {e}|{abx}|{flag}|{gbid}|{_dfvx_ss.iloc[abx].h3_name}|{gbfp}|{_dfvx_ss.iloc[abx].vh_abi_fp}|{score}|{quality_score}|{lowest_quality}")if log_msg else None
#                                     filtered_res.append([flag, gbid, _dfvx_ss.iloc[abx].h3_name, gbfp, _dfvx_ss.iloc[abx].vh_abi_fp, score, quality_score, lowest_quality])
#                     else:
#                         log.info(f"[+] FILTERED_BADMATCH: {e}|{abx}|{flag}|{gbid}|{_dfvx_ss.iloc[abx].h3_name}|{gbfp}|{_dfvx_ss.iloc[abx].vh_abi_fp}|{score}|---|---")if log_msg else None
#             except:
#                 pass
#     if not len(filtered_res) == 0:
#         dfres = pd.DataFrame(filtered_res)
#         dfres.columns = ["Orient", "gbid", "H3_name", "GB_FP", "ABI_FP", "Score", "Quality_score", "Low_quality"]
#         return dfres
#     else:
#         None
        

# def run_gb_alignment_and_filtering(M_gb_abi_vx, df_vx, vx_template_gb_dict, vx_seq_start_a, vx_seq_end_a, par_match, par_missmatch, par_open, par_extend, par_filter_thresh, log, log_msg=True, is_data_vl=False):
#     filtered_res = []
#     for e,i in enumerate(range(M_gb_abi_vx.shape[0])[0:]):
#         gbid = M_gb_abi_vx.iloc[i].gb_id
#         gbfp = M_gb_abi_vx.iloc[i].gb_fp
#         sample_id = M_gb_abi_vx.iloc[i].sample_id
#         _dfvx_ss = df_vx[df_vx.sample_id == M_gb_abi_vx.iloc[i].sample_id]
#         for e2, abx in enumerate(range(_dfvx_ss.shape[0])):
#             # genbank
#             seqA = SeqIO.read(vx_template_gb_dict.get(gbid), 'gb')
#             # vl
#             flag = ''
#             if is_data_vl:
#                 flag = 'VL'
#                 seqB = get_seqobj_from_abi(_dfvx_ss.iloc[abx].vl_abi_fp)
#                 seqB_vxinit = _dfvx_ss.iloc[abx]["vl_inti_sr"].seq
#                 seqB_vx, quality = get_trimmed_seq_record(seqB, get_quality=True)
#                 seqB_vx = get_seq_from_record(seqB_vx, reverse=True)
#             else:
#                 flag = 'VH'                
#                 seqB_vx = _dfvx_ss.iloc[abx].vh_sr_seq_r
#                 seqB_vxinit = _dfvx_ss.iloc[abx].vh_init_sr.seq
#                 quality = _dfvx_ss.iloc[abx].vh_sr_tqlst
#             # match between gb and vl-abi
#             try:          
#                 for a in align(seqA.seq, seqB_vx, match=par_match, mismatch=par_missmatch, gap_open=par_open, gap_extend=par_extend):
#                     score = int(a[2])
#                     if score >= int(par_filter_thresh):
#                         # getting the aligned seqA and seqB
#                         aligned_seq_a = a[0]
#                         aligned_seq_a = aligned_seq_a.replace('-', '')
#                         aligned_seq_a_start = aligned_seq_a[0:6]
#                         aligned_seq_a_length = len(aligned_seq_a) - 6
#                         aligned_seq_a_end = aligned_seq_a[aligned_seq_a_length:len(aligned_seq_a)]
#                         aligned_seq_b = a[1]
#                         # futher filtering of seqB
#                         if '-' not in aligned_seq_b:
#                             if vx_seq_start_a in aligned_seq_a_start:
#                                 if vx_seq_end_a in aligned_seq_a_end:
#                                     aligned_seq = a[0]
#                                     aligned_seq = aligned_seq.replace('-', '')
#                                     aligned_seq = Seq(aligned_seq)
#                                     aligned_seq_rev = aligned_seq.reverse_complement()
#                                     dna_string = str(seqB_vxinit)
#                                     aligned_seq_rev = str(aligned_seq_rev)
#                                     start = re.search(aligned_seq_rev, dna_string).start()
#                                     end = re.search(aligned_seq_rev, dna_string).end()
#                                     quality = quality[start:end]
#                                     quality_score = statistics.mean(quality)
#                                     lowest_quality = min(quality)
#                                     # print(f"{e} | {abx} | {gbid} | {_dfvx_ss.iloc[abx].h3_name} | {gbfp} | {_dfvx_ss.iloc[abx].vh_abi_fp} | {score}| {quality_score} | {lowest_quality}")
#                                     log.info(f"[+] FILTERED_GOODMATCH: {e}|{abx}|{flag}|{gbid}|{_dfvx_ss.iloc[abx].h3_name}|{sample_id}|{gbfp}|{_dfvx_ss.iloc[abx].vh_abi_fp}|{score}|{quality_score}|{lowest_quality}")if log_msg else None
#                                     filtered_res.append([flag, gbid, _dfvx_ss.iloc[abx].h3_name, sample_id, gbfp, _dfvx_ss.iloc[abx].vh_abi_fp, score, quality_score, lowest_quality])
#                     else:
#                         log.info(f"[+] FILTERED_BADMATCH: {e}|{abx}|{flag}|{gbid}|{_dfvx_ss.iloc[abx].h3_name}|{sample_id}|{gbfp}|{_dfvx_ss.iloc[abx].vh_abi_fp}|{score}|---|---")if log_msg else None
#             except:
#                 pass
#     if not len(filtered_res) == 0:
#         dfres = pd.DataFrame(filtered_res)
#         dfres.columns = ["Orient", "gbid", "H3_name", "sample_id", "GB_FP", "ABI_FP", "Score", "Quality_score", "Low_quality"]
#         return dfres
#     else:
#         None
        

def run_gb_alignment_and_filtering_v2(M_gb_abi_vx, 
                                   df_vx, 
                                   vx_template_gb_dict, 
                                   vx_seq_start_a, 
                                   vx_seq_end_a, 
                                   par_match, 
                                   par_missmatch, 
                                   par_open, 
                                   par_extend, 
                                   par_filter_thresh, 
                                   log,
                                   log_msg=True,
                                   is_data_vl=False):
    filtered_res = []
    for e,i in enumerate(range(M_gb_abi_vx.shape[0])[0:]):
        gbid = M_gb_abi_vx.iloc[i].gb_id
        gbfp = M_gb_abi_vx.iloc[i].gb_fp
        sample_id = M_gb_abi_vx.iloc[i].sample_id
        _dfvx_ss = df_vx[(df_vx.sample_id == M_gb_abi_vx.iloc[i].sample_id )]
        if (_dfvx_ss.h3_name == M_gb_abi_vx.iloc[i].h3_name).all():
            for e2, abx in enumerate(range(_dfvx_ss.shape[0])):
                # genbank
                seqA = SeqIO.read(vx_template_gb_dict.get(gbid), 'gb')
                # vl
                flag = ''
                if is_data_vl:
                    flag = 'VL'
                    seqB = get_seqobj_from_abi(_dfvx_ss.iloc[abx].vl_abi_fp)
                    seqB_vxinit = _dfvx_ss.iloc[abx]["vl_inti_sr"].seq
                    seqB_vx, quality = get_trimmed_seq_record(seqB, get_quality=True)
                    seqB_vx = get_seq_from_record(seqB_vx, reverse=True)
                else:
                    flag = 'VH'                
                    seqB_vx = _dfvx_ss.iloc[abx].vh_sr_seq_r
                    seqB_vxinit = _dfvx_ss.iloc[abx].vh_init_sr.seq
                    quality = _dfvx_ss.iloc[abx].vh_sr_tqlst
                # match between gb and vl-abi
                try:          
                    for a in align(seqA.seq, seqB_vx, match=par_match, mismatch=par_missmatch, gap_open=par_open, gap_extend=par_extend):
                        score = int(a[2])
                        if score >= int(par_filter_thresh):
                            # getting the aligned seqA and seqB
                            aligned_seq_a = a[0]
                            aligned_seq_a = aligned_seq_a.replace('-', '')
                            aligned_seq_a_start = aligned_seq_a[0:6]
                            aligned_seq_a_length = len(aligned_seq_a) - 6
                            aligned_seq_a_end = aligned_seq_a[aligned_seq_a_length:len(aligned_seq_a)]
                            aligned_seq_b = a[1]
                            # futher filtering of seqB
                            if '-' not in aligned_seq_b:
                                if vx_seq_start_a in aligned_seq_a_start:
                                    if vx_seq_end_a in aligned_seq_a_end:
                                        aligned_seq = a[0]
                                        aligned_seq = aligned_seq.replace('-', '')
                                        aligned_seq = Seq(aligned_seq)
                                        aligned_seq_rev = aligned_seq.reverse_complement()
                                        dna_string = str(seqB_vxinit)
                                        aligned_seq_rev = str(aligned_seq_rev)
                                        start = re.search(aligned_seq_rev, dna_string).start()
                                        end = re.search(aligned_seq_rev, dna_string).end()
                                        quality = quality[start:end]
                                        quality_score = statistics.mean(quality)
                                        lowest_quality = min(quality)
                                        # print(f"{e} | {abx} | {gbid} | {_dfvx_ss.iloc[abx].h3_name} | {gbfp} | {_dfvx_ss.iloc[abx].vh_abi_fp} | {score}| {quality_score} | {lowest_quality}")
                                        log.info(f"[+] FILTERED_GOODMATCH: {e}|{abx}|{flag}|{gbid}|{_dfvx_ss.iloc[abx].h3_name}|{sample_id}|{gbfp}|{_dfvx_ss.iloc[abx].vh_abi_fp}|{score}|{quality_score}|{lowest_quality}")if log_msg else None
                                        filtered_res.append([flag, gbid, _dfvx_ss.iloc[abx].h3_name, sample_id, gbfp, _dfvx_ss.iloc[abx].vh_abi_fp, score, quality_score, lowest_quality])
                        else:
                            log.info(f"[+] FILTERED_BADMATCH: {e}|{abx}|{flag}|{gbid}|{_dfvx_ss.iloc[abx].h3_name}|{sample_id}|{gbfp}|{_dfvx_ss.iloc[abx].vh_abi_fp}|{score}|---|---")if log_msg else None
                except:
                    pass
    if not len(filtered_res) == 0:
        dfres = pd.DataFrame(filtered_res)
        dfres.columns = ["Orient", "gbid", "H3_name", "sample_id", "GB_FP", "ABI_FP", "Score", "Quality_score", "Low_quality"]
        return dfres
    else:
        None
        
        
def cal_mean_error_prob(vh_low_qual:int, vl_low_qual:int):
    return (10**(-vh_low_qual/10) + 10**(-vl_low_qual/10))/2        

def main():
    pass

if __name__ == '__main__':
    main()