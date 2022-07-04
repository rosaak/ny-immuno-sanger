import configparser
import logging
import os
import re
import shutil
import statistics
from pathlib import Path

import pandas as pd
import utils
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq

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

def get_analysis_abbr(config: configparser.ConfigParser):
    return config['Analysis']['analysis_abbr']

def get_abi_extension(config: configparser.ConfigParser):
    return config['ABI']['abi_extension']

def get_gb_extension(config: configparser.ConfigParser):
    return config['GENBANK']['gb_extension']

def get_gb_pattern_to_filter(config: configparser.ConfigParser):
    return config['GENBANK']['gb_pat_to_filter']

def get_patterns(config: configparser.ConfigParser):
    return config['Pattern']['pat_vx']

def get_patterns_to_rm(config: configparser.ConfigParser):
    return config['Pattern']['vx_pat_to_rm']

def get_patterns_to_rm_from_genbank(config: configparser.ConfigParser):
    return config['Pattern']['vx_pat_to_rm_from_gb_name']

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
        
def checking_vx_ab1_files(log:logging, vx_abi_dict: dict) -> list[str]:
    """
    Return sample ids list
    """
    vx_keys = [i for i in [*vx_abi_dict.keys()]]
    
    log.info(f"[+] INFO: The number of samples in vx_abi_dict {len(vx_keys)}")
    return vx_keys

def get_abi_file_path(key:str, vh_abi_dict: dict, vl_abi_dict: dict) -> str:
    vh = vh_abi_dict.get(key)
    vl = vl_abi_dict.get(key)
    return vh, vl

def get_vx_abi_file_path(key:str, vx_abi_dict: dict) -> str:
    vx = vx_abi_dict.get(key)
    return vx

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

def find_match_on_all_h3probes(log: logging, h3_dict: dict,  vx_sr: SeqIO.SeqRecord, sample:str, vx_abi_dict: dict) -> list:
    """
    Returns: [orientation of match, h3_probe_name, sample_id, vx_abi_dp, probe_seq, vx_seq_record, vx_trimmed_seq_qual]
    inputs
    log:         logger 
    h3_dict:     dictionary containing h3_name and probe_name
    vx_sr:       Sequence Record from the abi file
    sample:      sample id - comes from sample_ids
    vx_abi_dict: dictionary containing vh names and full file path
    """
    results = []
    for h3key, h3val in h3_dict.items():
        vx_sr_trimmed, vx_sr_trimquallst = get_trimmed_seq_record(vx_sr, get_quality=True)        
        vx_sr_seq_f = get_seq_from_record(vx_sr, reverse=False)
        vx_sr_seq_r = get_seq_from_record(vx_sr, reverse=True)
        vx_fn = vx_abi_dict.get(sample)
        
        if h3val in vx_sr_seq_f:
            log.info(f"[+] PROBE_MATCHING: Forward_Strand_Match|{h3key}|{sample}|{vx_fn}")
            results.append(["Forward_Strand_Match", h3key, sample, vx_fn, h3val, vx_sr, vx_sr_seq_f, vx_sr_trimmed, vx_sr_trimquallst])
        elif h3val in vx_sr_seq_r:
            log.info(f"[+] PROBE_MATCHING: Reverse_Strand_Match|{h3key}|{sample}|{vx_fn}")
            results.append(["Reverse_Strand_Match", h3key, sample, vx_fn, h3val, vx_sr, vx_sr_seq_r, vx_sr_trimmed, vx_sr_trimquallst])
        else:
            pass
    return(results)

def copy_mtched_abi_files_to_resdir_vx(log: logging, res_dir_vx: str, abbr_lower:str, df_vx:pd.DataFrame=None, log_msg:bool= True)->pd.DataFrame:
    """
    creates a dir within the results_dir for VX, with h3_name
    and copy matched vx abi files into it
    """
    def _cp_files(log: logging, df:pd.DataFrame, res_dir_vx:str, vx_abi_fp:str):
        res = []
        lst = df.h3_name.unique()
        for i in lst:
            _df = df[df.h3_name == i]
            
            for j in range(_df.shape[0]):
                __df = _df.iloc[j]
                res_dir_msid = checking_dirs(f"{res_dir_vx}/{__df.h3_name}", log, log_msg=False, create_dir=True)
                try:
                    shutil.copy(__df[vx_abi_fp], res_dir_msid)
                    log.info(f"[+] COPY: {__df[vx_abi_fp]} to {res_dir_msid}") if log_msg else None
                    res.append([__df.h3_name, __df.sample_id, f"{res_dir_vx}/{__df.sample_id}", __df[vx_abi_fp]])
                except:
                    log.info(f"{e}")
        return pd.DataFrame(res, columns=["h3_name", "sample_id", "abi_out_loc", "abi_initial_filepath"])
    
    df1 = _cp_files(log, df_vx, res_dir_vx, f"{abbr_lower}_abi_fp")
    return df1

def find_gb_match_on_all_h3probes_single(log: logging, gb_tmplate_dict: dict, df_vx: pd.DataFrame, pattern='VH-', log_msg: bool=True) -> pd.DataFrame:
    gb_match_with_dfvx_lst = []
    flag_col_name = ''
    if 'vh' in pattern.lower():
        flag_col_name = 'vh_abi_fp'
    elif 'vl' in pattern.lower():
        flag_col_name = 'vl_abi_fp'
    for gb_id in [*gb_tmplate_dict.keys()]:
        _ = df_vx[df_vx.h3_name == gb_id.replace(pattern, '')]
        if not (_).empty:
            for j in range(_.shape[0]):
                log.info(f"[+] GENBANK_MATCHING: {gb_id}|{_.iloc[j].h3_name}|{_.iloc[j].sample_id}|{gb_tmplate_dict.get(gb_id)}|{_.iloc[j][flag_col_name]}") if log_msg else None
                gb_match_with_dfvx_lst.append(['GB_Matched',gb_id, _.iloc[j].h3_name, _.iloc[j].sample_id, gb_tmplate_dict.get(gb_id), _.iloc[j][flag_col_name]  ])
        else:
            pass
    _df = pd.DataFrame(gb_match_with_dfvx_lst)
    _df.columns = ["Match", "gb_id", "h3_name", "sample_id", "gb_fp", flag_col_name]
    return _df

def align(A_seq: SeqIO.SeqRecord, B_seq: SeqIO.SeqRecord, match :int, mismatch: int, gap_open: int, gap_extend: int):
    return pairwise2.align.globalms(sequenceA=A_seq, sequenceB=B_seq, match=match, mismatch=mismatch, open=gap_open, extend=gap_extend, penalize_end_gaps = False)

def run_gb_alignment_and_filtering_single(M_gb_abi_vx, 
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
        dfres.columns = ["Orient", "gbid", "H3_name", "sample_id", "GB_FP", "ABI_FP", "Score", "Mean_Quality_score", "Low_quality"]
        return dfres
    else:
        None

def cal_mean_error_prob(vh_low_qual:int, vl_low_qual:int):
    return (10**(-vh_low_qual/10) + 10**(-vl_low_qual/10))/2        

def main():
    pass

if __name__ == '__main__':
    main()