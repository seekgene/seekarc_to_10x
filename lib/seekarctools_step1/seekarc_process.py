#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import dnaio
import itertools
from collections import defaultdict, Counter
from seekarc_utils import AdapterFilter, logger, get_new_bc, cmd_execute

def cut_fq(fq:str, outdir:str, reads_num:int=400000):
    fq_tmp = os.path.join(outdir, os.path.basename(fq))
    cmd = f'zcat "{fq}"|head -n {reads_num}|gzip > {fq_tmp}'
    cmd_execute(cmd, check=True)
    return fq_tmp

def barcode_report(logfile):
    with open(logfile) as fh:
        summary = json.load(fh)
    d = summary["stat"]
    return float(d["valid"])/d["total"]

def summary(seq, seq_q, seq_dict, qua_dict):
    for i, (base, q) in enumerate(zip(seq, seq_q)):
        seq_dict[(i,base)] += 1
        qua_dict[(i,q)] += 1
    return seq_dict, qua_dict

def process_barcode(fq1, fq2, fq_out, fqout_multi, r1_structure, shift, shift_pattern,
                    barcode_wl_dict, linker_wl_dict, match_type_dict, adapter1=[["AAAAAAAAAAAA", "3"],],
                    adapter2=[["AAAAAAAAAAAA", "3"],], do_B_correction=True, do_L_correction=True,
                    use_multi=True, use_short_read=False, paired_out=True, use_cr=False,
                    min_len1=20, min_len2=20):
    
    barcode_list_flag = False
    linker_list_flag = False
    if len(barcode_wl_dict)>0:
        barcode_list_flag = True

    if len(linker_wl_dict)>0:
        linker_list_flag = True

    stat_Dict = defaultdict(int)
    Barcode_Counter = Counter()
    
    Barcode_GC_Counter = Counter()
    UMI_GC_Counter = Counter()
    R2_GC_Counter = Counter()
    Barcode_Q_Counter = Counter()
    UMI_Q_Counter = Counter()
    R2_Q_Counter = Counter()
    
    adapter_filter = AdapterFilter(adapter1=adapter1, adapter2=adapter2)
    
    fh = dnaio.open(fq1, fq2, fileformat="fastq", mode="r")
    if paired_out:
        outfh = dnaio.open(fq_out[0], fq_out[1], fileformat="fastq", mode="w")
    else:
        outfh = dnaio.open(fq_out[0], fileformat="fastq", mode="w")

    if use_multi:
        if paired_out:
            outfh_multi = dnaio.open(fqout_multi[0], fqout_multi[1], fileformat="fastq", mode="w")
        else:
            outfh_multi = dnaio.open(fqout_multi[0], fileformat="fastq", mode="w")
    
    for r1, r2 in fh:
        stat_Dict["total"] += 1
        
        start_pos = 0
        end_pos = 0
        sequence = r1.sequence
        qualities = r1.qualities
        
        # deal with shift
        if shift:
            shift_pos = sequence[:7].find(shift_pattern)
            if shift_pos < 0:
                stat_Dict["no_anchor"] += 1
                logger.debug(f"{r1.name},{sequence},{sequence[:7]} no anchor!")
                continue
            else:
                start_pos = shift_pos + 1
        
        # get barcode/umi/quality sequence          
        old_seqs = defaultdict(list)
        new_seqs = defaultdict(list)
        seq_quals = defaultdict(list)
        B = 0
        L = 0
        is_valid = True
        is_multi = False
        is_correct = False
        is_B_no_correction = False
        is_L_no_correction = False

        for _, (code, n) in enumerate(r1_structure):
            end_pos = start_pos + n
            seq = sequence[start_pos:end_pos]
            quals = qualities[start_pos:end_pos]

            if code == "B":
                old_seqs["B"].append(seq)
                seq_quals["B"].append(quals)

                if barcode_list_flag: # match barcode in whitelist
                    if seq in barcode_wl_dict.get(B, barcode_wl_dict[0]):
                        new_seqs["B"].append({seq})
                    else:
                        if do_B_correction:
                            bc_set = get_new_bc(seq, barcode_wl_dict.get(B, barcode_wl_dict[0]), match_type_dict.get(B, match_type_dict[0]))

                            if len(bc_set) == 0:
                                is_valid = False
                                is_B_no_correction = True
                                logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} no barcode!")
                                break
                            elif len(bc_set) == 1:
                                new_seqs["B"].append(bc_set)
                                is_correct = True
                                logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} -> {list(bc_set)} do_B_correction!")
                            else:
                                new_seqs["B"].append(bc_set)
                                logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} -> {list(bc_set)} do_B_correction!")
                                is_multi = True
                        else:
                            is_valid = False
                            break
                else:
                    new_seqs["B"].append({seq})
                B += 1

            elif code == "L":   
                if linker_list_flag: # linker correction step
                    if seq in linker_wl_dict.get(L, linker_wl_dict[0]):
                        pass
                    else:
                        if do_L_correction:
                            lk_set = get_new_bc(seq, linker_wl_dict.get(L, linker_wl_dict[0]))
                            if len(lk_set) == 0:
                                is_valid = False
                                is_L_no_correction = True
                                logger.debug(f"{r1.name},{sequence[:start_pos]}[{seq}]{sequence[end_pos:]},{seq} -> {list(lk_set)} no linker!")
                                break
                        else:
                            is_valid = False
                            break
                L += 1
                
            elif code == "U":
                new_seqs["U"].append(seq)
                seq_quals["U"].append(quals)
                
            start_pos = start_pos + n

        # check double instances
        if is_valid:

            barcode_old = "".join(old_seqs["B"])
            Barcode_Counter[barcode_old] += 1

            #get base summary for umi/r2
            umi = "".join(new_seqs["U"])
            umi_q = "".join(seq_quals["U"])
            barcode_q = "".join(seq_quals["B"])
            
            UMI_GC_Counter, UMI_Q_Counter = summary(umi, umi_q, UMI_GC_Counter, UMI_Q_Counter)
            R2_GC_Counter, R2_Q_Counter = summary(r2.sequence, r2.qualities, R2_GC_Counter, R2_Q_Counter)
            
            r1.sequence = sequence[start_pos:]
            r1.qualities = qualities[start_pos:]
                        
            if is_multi: #write r2 multi files
                if use_multi:         
                    #update barcode quality
                    Barcode_Q_Counter.update(enumerate(barcode_q))
                    bc_new_lst = []
                    for element in itertools.product(*new_seqs["B"]):          
                        barcode_new = "".join(element)
                        bc_new_lst.append(barcode_new)
                        
                    bc_new_all = ":".join(bc_new_lst)
                    r2.name = "_".join([barcode_old, bc_new_all, umi, r2.name])
                    r1.name = "_".join([barcode_old, bc_new_all, umi, r1.name])
                    r1.sequence = sequence[start_pos:]
                    r1.qualities = qualities[start_pos:]
                    outfh_multi.write(r1, r2)
            else:  #write r2 files
                stat_Dict["valid"] += 1
                flag, r1, r2 = adapter_filter.filter(r1, r2)
                if flag:
                    if (not use_short_read) or len(r1) == 0 or len(r2) == 0:
                        if len(r1) < min_len1 or len(r2) < min_len2:
                            stat_Dict["too_short"] += 1
                            continue
                    else:
                        stat_Dict["trimmed"] += 1

                barcode_new = "".join([_.pop() for _ in new_seqs["B"]])
                Barcode_GC_Counter, Barcode_Q_Counter = summary(barcode_old, barcode_q, Barcode_GC_Counter, Barcode_Q_Counter)
                
                #find alterations
                if is_correct:
                    _alt = "".join([str(i)+o for i, (o,n) in enumerate(zip(barcode_old, barcode_new)) if o != n])
                else:
                    _alt = "M"

                if use_cr:
                    r2.name = r2.name
                    r1.name = r1.name
                    r1.sequence = barcode_new+umi+r1.sequence
                    r1.qualities = qualities[:len(barcode_new+umi)]+r1.qualities
                else:
                    r2.name = "_".join([barcode_new, umi, _alt, r2.name])
                    r1.name = "_".join([barcode_new, umi, _alt, r1.name])
                outfh.write(r1, r2)
            if is_correct:
                stat_Dict["B_corrected"] += 1
        else:
            if is_B_no_correction:
                stat_Dict["B_no_correction"] += 1

            if is_L_no_correction:
                stat_Dict["L_no_correction"] += 1
    if use_multi:
        outfh_multi.close()
    outfh.close()

    return {
            "stat": Counter(stat_Dict),
            "barcode_count": Barcode_Counter,
            "barcode_gc": Barcode_GC_Counter,
            "umi_gc": UMI_GC_Counter,
            "r2_gc": R2_GC_Counter,
            "barcode_q": Barcode_Q_Counter,
            "umi_q": UMI_Q_Counter,
            "r2_q": R2_Q_Counter
        }
