#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Standalone script for SeekARC ATAC Step 1 (Extraction)
Refactored to use shared modules.
"""

import os
import sys
import json
import logging
import argparse
import subprocess
# import multiprocessing as mp
from functools import partial
from collections import defaultdict

# Shared Modules
try:
    from seekarc_config import CHEMISTRY, ATAC_R1_MINLEN, ATAC_R2_MINLEN
    from seekarc_utils import (
        cmd_execute, check_path, parse_structure, read_file, 
        AdapterFilter, QcStat, logger
    )
    from seekarc_pipeline import Pipeline
    from seekarc_process import process_barcode, barcode_report, cut_fq
    from seekarc_atac_cutreads import cutreads
except ImportError:
    # Allow running if shared modules are in the same directory
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from seekarc_config import CHEMISTRY, ATAC_R1_MINLEN, ATAC_R2_MINLEN
    from seekarc_utils import (
        cmd_execute, check_path, parse_structure, read_file, 
        AdapterFilter, QcStat, logger
    )
    from seekarc_pipeline import Pipeline
    from seekarc_process import process_barcode, barcode_report, cut_fq
    from seekarc_atac_cutreads import cutreads

__version__ = "1.1.0"

# Directory of this script
__srcdir = os.path.dirname(os.path.abspath(__file__))

# ==========================================
# Helper Functions (Specific to this script or wrapping shared logic)
# ==========================================

def chemistry_auto(fq1, fq2, outdir, use_short_read=False, **kwargs):
    _outdir = f"{outdir}/.test"
    rawdata = f"{_outdir}/data"
    os.makedirs(rawdata, exist_ok=True)
    fq1_1M, fq2_1M = cut_fq(fq1[0], rawdata), cut_fq(fq2[0], rawdata)

    rate_dict = {}
    if kwargs["chemistry"] == "Auto":
        test_chems = ["MM", "DDV2"]
    elif kwargs["chemistry"] == "custom":
        CHEMISTRY["custom"] = {
            "shift": kwargs["shift"],
            "shift_pattern": kwargs["shift_pattern"],
            "structure": kwargs["structure"],
            "barcode": kwargs["barcode"],
            "linker": kwargs["linker"],
            "match_type": [1,],           
        }
        test_chems = ["custom",]
    else:
        test_chems = [kwargs["chemistry"],]

    for chem in test_chems:
        logger.info(f"test {chem}!")
        barcode_main([fq1_1M,], [fq2_1M,], chem, f"{_outdir}/{chem}", core=1, use_short_read=use_short_read, **CHEMISTRY[chem])
        rate = barcode_report(f"{_outdir}/{chem}/{chem}_A_summary.json")
        rate_dict[chem] = rate
    for k, v in rate_dict.items():
        logger.info(f"valid barcode rate of {k}: {v*100:.3f}%")
    if kwargs["chemistry"] == "Auto":
        if rate_dict["DDV2"] > rate_dict["MM"]:
            chemistry = "DDV2"
        else:
            chemistry = "MM"
    else:
        chemistry = kwargs["chemistry"]
    return CHEMISTRY[chemistry]

def check_atac_options(**kwargs):
    kwargs = chemistry_auto(**kwargs)
    return kwargs

def barcode_main(fq1:list, fq2:list, samplename: str, outdir:str,
                 barcode:list=[], match_type:list=[], shift:str=True, shift_pattern:str="A",
                 structure:str="B8L8B8L10B8U12T15", linker: list=[],
                 core:int=4, do_B_correction=True, do_L_correction=True,
                 use_multi=True, use_short_read=False, adapter1=[["TTTTTTTTTTTT", "5"], ],
                 adapter2=[["AAAAAAAAAAAA", "3"], ], paired_out=True, **kwargs):
    logger.info("extract barcode start!")
    r1_structure = parse_structure(structure)
    match_type_dict = {ind: val for ind, val in enumerate(match_type)}
    
    barcode_wl_dict = read_file(barcode)
    linker_wl_dict = read_file(linker)

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
    
    # Use shared process_barcode with ATAC specific min_len
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
        min_len1=ATAC_R1_MINLEN,
        min_len2=ATAC_R2_MINLEN
    )
    
    stat = QcStat()
    
    atacname = f"{samplename}_A"
    os.makedirs(f"{outdir}/step1", exist_ok=True)
    fqout = os.path.join(outdir, f"step1/{atacname}")
    fqout_multi = os.path.join(outdir, f"step1/{atacname}_multi")
    json_multi = os.path.join(outdir, f"step1/{atacname}_multi.json")
    
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
        
        # Note: dnaio.open is handled inside seekarc_utils/process logic, 
        # but here we are doing post-processing on the output files.
        # We need dnaio imported for this part.
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
                        if len(r1) < ATAC_R1_MINLEN or len(r2) < ATAC_R2_MINLEN:
                            multi_stat["too_short"] += 1
                            stat.data["stat"]["too_short"] += 1
                            continue
                    else:
                        multi_stat["trimmed"] += 1
                        stat.data["stat"]["trimmed"] += 1

                alt_l = [str(i)+o for i, (o,n) in enumerate(zip(bc_old, final_barcode)) if o != n]
                _alt = "".join([alt for alt in alt_l])
                r2.name = "_".join([final_barcode, umi, _alt, r2_name])

                r1.name = r2.name
                f.write(r1, r2)

        with open(json_multi, "w") as fh:
            json.dump(multi_stat, fp=fh, indent=4)

    if "barcode_count" in stat.data:
        del stat.data["barcode_count"]
    logger.info("deal multi done!")
    stat.data["stat"]["chemistry"] = kwargs.get("chemistry", "custom")
    stat.data["stat"]["atacname"] = atacname
    stat.save(os.path.join(outdir, f"{atacname}_summary.json"))
    logger.info("extract barcode done!")
    return fqout1, fqout2

def cutreads(step1afq1, step1afq2, atacjson, outdir, atacname, **kwargs):
    # Delegate to the standalone module
    from seekarc_atac_cutreads import cutreads as _cutreads
    return _cutreads(step1afq1, step1afq2, atacjson, outdir, atacname, **kwargs)


# ==========================================
# Main CLI
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="SeekARC ATAC Step 1 (Standalone)")
    parser.add_argument("--fq1", required=True, nargs='+', help="Read 1 FASTQ file(s)")
    parser.add_argument("--fq2", required=True, nargs='+', help="Read 2 FASTQ file(s)")
    parser.add_argument("--samplename", required=True, help="Sample name")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--barcode", help="Path to barcode whitelist file", nargs='+')
    parser.add_argument("--linker", help="Path to linker whitelist file", nargs='+')
    parser.add_argument("--match_type", nargs='+', default=[1,], help="Match type (e.g. 1 3)")
    parser.add_argument("--core", type=int, default=4, help="Number of cores (default: 4)")
    parser.add_argument("--structure", help="Read structure (default: B16U12)")
    
    # Optional arguments matching seekarctools
    parser.add_argument("--chemistry", default="Auto", help="Chemistry (default: Auto)")
    parser.add_argument("--shift", action="store_true", help="Shift reads")
    parser.add_argument("--shift-pattern", default="A", help="Shift pattern")
    
    # Additional flags
    parser.add_argument("--do-B-correction", action="store_true", default=True)
    parser.add_argument("--no-do-B-correction", dest="do_B_correction", action="store_false")
    parser.add_argument("--do-L-correction", action="store_true", default=True)
    parser.add_argument("--no-do-L-correction", dest="do_L_correction", action="store_false")
    parser.add_argument("--use-multi", action="store_true", default=True)
    parser.add_argument("--no-use-multi", dest="use_multi", action="store_false")
    parser.add_argument("--use-short-read", action="store_true", default=False)
    parser.add_argument("--paired-out", action="store_true", default=True)

    args = parser.parse_args()
    
    # # Resolve chemistry
    # chem_config = {}
    # if args.chemistry in CHEMISTRY:
    #     chem_config = CHEMISTRY[args.chemistry]
    # elif args.chemistry == "Auto":
    #     # Default fallback if Auto logic is complex; here we just warn or use a default
    #     # For standalone, maybe default to 'DDV1' or require explicit chemistry?
    #     # Let's use DDV1 as a safe default if Auto
    #     chem_config = CHEMISTRY.get("DDV1", {})
    #     logger.info("Chemistry set to Auto, defaulting to DDV1 config.")
    
    # # Override with command line args if provided
    # kwargs = vars(args)
    
    # # Merge chemistry config
    # for k, v in chem_config.items():
    #     if k not in kwargs or kwargs[k] is None:
    #         kwargs[k] = v
    #     # If structure is not provided in args but exists in chemistry, use it
    #     if k == 'structure' and args.structure is None:
    #          kwargs['structure'] = v
    #     # If barcode is not provided in args but exists in chemistry, use it
    #     if k == 'barcode' and args.barcode is None:
    #          kwargs['barcode'] = v
    #     if k == 'linker' and args.linker is None:
    #          kwargs['linker'] = v

    # # Ensure structure is set
    # if not kwargs.get('structure'):
    #     kwargs['structure'] = "B16U12" # Fallback
    kwargs = vars(args)
    fq1_out, fq2_out = barcode_main(**kwargs)
    
    # Execute cutreads
    # Setup paths for cutreads
    kwargs["atacname"] = f'{kwargs["samplename"]}_A'
    kwargs["atacjson"] = os.path.join(kwargs["outdir"], f"{kwargs['atacname']}_summary.json")
    kwargs["step1afq1"] = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_1.fq.gz")
    kwargs["step1afq2"] = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_2.fq.gz")

    # Run Cutreads
    logger.info("Running cutreads...")
    cutreads(**kwargs)

# if __name__ == "__main__":
#     # Basic Argument Parser
#     parser = argparse.ArgumentParser(description="SeekARC ATAC Step 1 (Standalone)")
#     parser.add_argument("--fq1", required=True, nargs='+', help="Read 1 FASTQ file(s)")
#     parser.add_argument("--fq2", required=True, nargs='+', help="Read 2 FASTQ file(s)")
#     parser.add_argument("--samplename", required=True, help="Sample name")
#     parser.add_argument("--outdir", required=True, help="Output directory")
#     parser.add_argument("--barcode", help="Path to barcode whitelist file", nargs='+')
#     parser.add_argument("--linker", help="Path to linker whitelist file", nargs='+')
#     parser.add_argument("--match_type", nargs='+', default=[], help="Match type (e.g. 1 3)")
#     parser.add_argument("--core", type=int, default=4, help="Number of cores (default: 4)")
#     parser.add_argument("--structure", help="Read structure")
    
#     # Optional arguments matching seekarctools
#     parser.add_argument("--chemistry", default="Auto", help="Chemistry (default: Auto)")
#     parser.add_argument("--shift", action="store_true", help="Shift reads")
#     parser.add_argument("--shift-pattern", default="A", help="Shift pattern")
    
#     # Additional flags
#     parser.add_argument("--do-B-correction", action="store_true", default=True)
#     parser.add_argument("--no-do-B-correction", dest="do_B_correction", action="store_false")
#     parser.add_argument("--do-L-correction", action="store_true", default=True)
#     parser.add_argument("--no-do-L-correction", dest="do_L_correction", action="store_false")
#     parser.add_argument("--use-multi", action="store_true", default=True)
#     parser.add_argument("--no-use-multi", dest="use_multi", action="store_false")
#     parser.add_argument("--use-short-read", action="store_true", default=False)
#     parser.add_argument("--paired-out", action="store_true", default=True)

#     args = parser.parse_args()
    
#     # Resolve chemistry
#     chem_config = {}
    
#     # Run check options (which runs auto-detection if needed)
#     # We construct a kwargs dict to pass to check_atac_options
#     kwargs = vars(args)
    
#     # Ensure structure is set if not provided (will be overwritten by chemistry if auto)
#     if not kwargs.get('structure'):
#         kwargs['structure'] = "B16U12" # Fallback

#     logger.info("Starting SeekARC ATAC Step 1")
#     logger.info(f"Sample: {args.samplename}")
#     logger.info(f"Outdir: {args.outdir}")

#     # Check/Auto-detect chemistry
#     # check_atac_options calls chemistry_auto which returns the config dict
#     chem_config = check_atac_options(**kwargs)
    
#     # Update kwargs with chemistry config (only if not manually overridden)
#     # Actually, chemistry_auto returns the config for the detected/selected chemistry.
#     # We should merge this into kwargs.
#     for k, v in chem_config.items():
#         # If the argument was not provided by user (None) or if we want chemistry to override
#         # Usually user args override chemistry defaults.
#         # args.structure is in kwargs. If it's None (default), we take from chemistry.
#         if kwargs.get(k) is None:
#              kwargs[k] = v
#         # Special handling for lists like barcode/linker which might be empty list [] in args
#         if k in ['barcode', 'linker'] and not kwargs.get(k):
#              kwargs[k] = v

#     # Execute barcode_main
#     fq1_out, fq2_out = barcode_main(**kwargs)
    
#     # Execute cutreads
#     atacname = f"{args.samplename}_A"
#     summary_json = os.path.join(args.outdir, f"step1/{atacname}_summary.json")
#     try:
#         cutreads(fq1_out, fq2_out, summary_json, args.outdir, atacname)
#     except Exception as e:
#         logger.error(f"Error in cutreads: {e}")
#         # Don't fail the whole pipeline if cutreads fails, but log it.


#     # Save parameters
#     os.makedirs(kwargs['outdir'], exist_ok=True)
#     with open(os.path.join(kwargs['outdir'], ".params.json"), "w") as fh:
#         json.dump(kwargs, fh, indent=4)

#     # Run Barcode Main
#     barcode_main(**kwargs)

#     # Setup paths for cutreads
#     kwargs["atacname"] = f'{kwargs["samplename"]}_A'
#     kwargs["atacjson"] = os.path.join(kwargs["outdir"], f"{kwargs['atacname']}_summary.json")
#     kwargs["step1afq1"] = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_1.fq.gz")
#     kwargs["step1afq2"] = os.path.join(kwargs["outdir"], "step1", f"{kwargs['atacname']}_2.fq.gz")

#     # Run Cutreads
#     logger.info("Running cutreads...")
#     cutreads(**kwargs)
    
#     logger.info("Done.")

if __name__ == "__main__":
    main()
