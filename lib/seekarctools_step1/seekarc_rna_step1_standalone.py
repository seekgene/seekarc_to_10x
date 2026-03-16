#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Standalone script for SeekARC RNA Step 1 (Extraction)
Mimics seekarctools RNA step 1 but uses shared modules.
"""

import os
import sys
import json
import logging
import argparse
import shutil
import subprocess
from functools import partial
from collections import defaultdict

# Shared Modules
try:
    from seekarc_config import CHEMISTRY, R1_MINLEN, R2_MINLEN
    from seekarc_utils import (
        cmd_execute, check_path, parse_structure, read_file, 
        AdapterFilter, QcStat, logger
    )
    from seekarc_pipeline import Pipeline
    from seekarc_process import process_barcode, barcode_report, cut_fq
except ImportError:
    # Allow running if shared modules are in the same directory
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from seekarc_config import CHEMISTRY, R1_MINLEN, R2_MINLEN
    from seekarc_utils import (
        cmd_execute, check_path, parse_structure, read_file, 
        AdapterFilter, QcStat, logger
    )
    from seekarc_pipeline import Pipeline
    from seekarc_process import process_barcode, barcode_report, cut_fq

__version__ = "1.1.0"

# Directory of this script
__srcdir = os.path.dirname(os.path.abspath(__file__))

def barcode_main(fq1:list, fq2:list, samplename: str, outdir:str,
                 barcode:list=[], match_type:list=[], shift:str=False, shift_pattern:str="A",
                 structure:str="B8L8B8L10B8U12T15", linker: list=[],
                 core:int=4, do_B_correction=True, do_L_correction=True,
                 use_multi=True, use_short_read=False, adapter1=[["TTTTTTTTTTTT", "5"], ],
                 adapter2=[["AAAAAAAAAAAA", "3"], ], paired_out=True, use_cr=False, **kwargs):
    logger.info("extract barcode start!")
    r1_structure = parse_structure(structure)
    
    barcode_wl_dict = read_file(barcode)
    linker_wl_dict = read_file(linker)
    match_type_dict = {ind: val for ind, val in enumerate(match_type)}

    if len(barcode_wl_dict)>0 and do_B_correction:
        logger.info("barcode one base mismatch allowed.")
    else:
        logger.info("barcode mismatch NOT allowed.")

    if "L" in structure:
        if len(linker_wl_dict)>0 and do_L_correction:
            logger.info("linker one base mismatch allowed.")
        else:
            logger.info("linker mismatch NOT allowed.")

    if use_multi:
        logger.info("rescue barcode match multi barcode in whitelist.")
    else:
        logger.info("ignore barcode match multi barcode in whitelist.")
    
    # Use shared process_barcode with RNA specific min_len and use_cr
    worker_func = partial(
        process_barcode,
        r1_structure=r1_structure,
        shift=shift,
        shift_pattern=shift_pattern,
        barcode_wl_dict=barcode_wl_dict,
        linker_wl_dict=linker_wl_dict,
        match_type_dict=match_type_dict,
        do_B_correction=do_B_correction,
        do_L_correction=do_L_correction,
        use_multi=use_multi,
        use_short_read=use_short_read,
        adapter1=adapter1,
        adapter2=adapter2,
        use_cr=use_cr,
        min_len1=R1_MINLEN,
        min_len2=R2_MINLEN
    )
    
    stat = QcStat()
    
    gexname = f"{samplename}_E"
    os.makedirs(f"{outdir}/step1", exist_ok=True)
    fqout = os.path.join(outdir, f"step1/{gexname}")
    fqout_multi = os.path.join(outdir, f"step1/{gexname}_multi")
    json_multi = os.path.join(outdir, f"step1/{gexname}_multi.json")
    
    pipeline = Pipeline(
        func=worker_func,
        fq1=fq1,
        fq2=fq2,
        fqout=fqout,
        fqout_multi=fqout_multi,
        stat=stat,
        core=core,
        paired_out=paired_out
    )
    pipeline.run()

    fqout1 = f"{fqout}_1.fq.gz"
    fqout2 = f"{fqout}_2.fq.gz"
    if use_multi:
        logger.info("deal multi start!")
        fqout_multi1 = f"{fqout_multi}_1.fq.gz"
        fqout_multi2 = f"{fqout_multi}_2.fq.gz"
        adapter_filter = AdapterFilter(adapter1=adapter1, adapter2=adapter2)
        multi_stat = defaultdict(int)
        
        import dnaio
        
        with dnaio.open(fqout1, fqout2, mode="a") as f:
            fh = dnaio.open(fqout_multi1, fqout_multi2, fileformat="fastq", mode="r")
            for r1, r2 in fh:
                multi_stat["total"] += 1
                final_barcode = None
                
                bc_old, r2_candidate, umi, r2_name = r2.name.split("_", 3)
                r2_candidate = r2_candidate.split(":")
                
                read_num = 0
                for _ in sorted(r2_candidate):
                    v = stat.data["barcode_count"].get(_, 0)
                    if v > read_num:
                        read_num = v
                        final_barcode = _

                if not final_barcode:
                    multi_stat["B_no_correction"] += 1
                    stat.data["stat"]["B_no_correction"] += 1
                    continue
                    
                multi_stat["valid"] += 1
                stat.data["stat"]["valid"] += 1

                flag, r1, r2 = adapter_filter.filter(r1, r2)
                if flag:
                    if (not use_short_read) or len(r1) == 0 or len(r2) == 0:
                        if len(r1) < R1_MINLEN or len(r2) < R2_MINLEN:
                            multi_stat["too_short"] += 1
                            stat.data["stat"]["too_short"] += 1
                            continue
                    else:
                        multi_stat["trimmed"] += 1
                        stat.data["stat"]["trimmed"] += 1

                alt_l = [str(i)+o for i, (o,n) in enumerate(zip(bc_old, final_barcode)) if o != n]
                _alt = "".join([alt for alt in alt_l])
                _, _, _, r1_name = r1.name.split("_", 3)
                if use_cr:
                    r2.name = r2_name
                    r1.name = r1_name
                    r1.sequence = final_barcode+umi+r1.sequence
                    r1.qualities = 'F'*len(final_barcode+umi)+r1.qualities
                else:
                    r2.name = "_".join([final_barcode, umi, _alt, r2_name])
                    r1.name = "_".join([final_barcode, umi, _alt, r1_name])
                f.write(r1, r2)

        with open(json_multi, "w") as fh:
            json.dump(multi_stat, fp=fh, indent=4)

    if "barcode_count" in stat.data:
        del stat.data["barcode_count"]
    logger.info("deal multi done!")
    stat.data["stat"]["chemistry"] = kwargs.get("chemistry", "custom")
    stat.data["stat"]["gexname"] = gexname
    stat.save(os.path.join(outdir, f"{gexname}_summary.json"))
    logger.info("extract barcode done!")
    return fqout1, fqout2

if __name__ == "__main__":
    # Basic Argument Parser
    parser = argparse.ArgumentParser(description="SeekARC RNA Step 1 (Standalone)")
    parser.add_argument("--fq1", required=True, nargs='+', help="R1 fastq files")
    parser.add_argument("--fq2", required=True, nargs='+', help="R2 fastq files")
    parser.add_argument("--samplename", required=True, help="Sample name")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--chemistry", default="Auto", help="Chemistry version")
    parser.add_argument("--core", type=int, default=4, help="Number of cores")
    
    # Optional arguments matching function signature
    parser.add_argument("--shift", default="store_true", help="Enable shift")
    parser.add_argument("--no-shift", dest="shift", action="store_false")
    parser.set_defaults(shift=False)
    
    parser.add_argument("--structure", default="B8L8B8L10B8U12T15", help="Read structure")
    parser.add_argument("--barcode", nargs='+', default=[], help="Barcode files")
    parser.add_argument("--linker", nargs='+', default=[], help="Linker files")
    parser.add_argument("--match_type", nargs='+', default=[1,], help="Match type (e.g. 1 3)")
    parser.add_argument("--use_cr", action="store_true", help="Use CR barcode")
    
    args = parser.parse_args()

    barcode_main(
        fq1=args.fq1,
        fq2=args.fq2,
        samplename=args.samplename,
        outdir=args.outdir,
        core=args.core,
        chemistry=args.chemistry,
        shift=args.shift,
        structure=args.structure,
        barcode=args.barcode,
        linker=args.linker,
        match_type=args.match_type,
        use_cr=args.use_cr
    )