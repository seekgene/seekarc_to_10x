# -*- coding: utf-8 -*-
import os
import json
import dnaio
import logging
import re
from seekarc_config import ATAC_R1_MINLEN, ATAC_R2_MINLEN

logger = logging.getLogger("seekarc_step1")

def cutreads(step1afq1, step1afq2, atacjson, outdir, atacname, **kwargs):
    total = 0
    available = 0
    outafq1 = os.path.join(outdir, 'step1', atacname+"_cutR1.fastq.gz")
    outafq2 = os.path.join(outdir, 'step1', atacname+"_cutR2.fastq.gz")
    outafq3 = os.path.join(outdir, 'step1', atacname+"_cutR3.fastq.gz")
    pattern = r'CGTCCGTCGTTGCTCGTAGATGTGTATAAGAGACAG'
    
    # Check if files exist
    if not os.path.exists(step1afq1) or not os.path.exists(step1afq2):
         logger.error(f"Input files for cutreads not found: {step1afq1}, {step1afq2}")
         return
         
    # Ensure output directory exists
    os.makedirs(os.path.dirname(outafq1), exist_ok=True)

    with dnaio.open(file1=step1afq1, file2=step1afq2,  mode='r') as fh, \
        dnaio.open(outafq1, mode='w') as fh1, \
        dnaio.open(outafq2, mode='w') as fh2, \
        dnaio.open(outafq3, mode='w') as fh3:
        
        for r1,r2 in fh:
            total += 1
            
            # r1 name format: "barcode_umi_alt_origName"
            tmp = r1.name.split('_')
            
            r1seq = r1.sequence
            r1qua = r1.qualities
            if len(r1seq) < 50:
                continue
            
            match = re.search(pattern, r1seq)
            if match:
                start_index = match.end()
                if start_index + 50 < len(r1seq):
                    r1_sequence = r1seq[start_index:start_index + 50]
                    r1_qualities = r1qua[start_index:start_index + 50]
                else:
                    r1_sequence = r1seq[start_index:]
                    r1_qualities = r1qua[start_index:]
            else:
                if len(r1seq) > 93:
                    r1_sequence = r1seq[43:93]
                    r1_qualities = r1qua[43:93]
                else:
                    r1_sequence = r1seq[43:]
                    r1_qualities = r1qua[43:]
            
            # R1 Output (Genomic)
            re1 = dnaio.Sequence(
                    name = tmp[-1],
                    sequence = r1_sequence,
                    qualities = r1_qualities
                    ) 
            
            # R2 Output (Barcode)
            re2 = dnaio.Sequence(
                    name = tmp[-1].replace(' 1:',' 2:'),
                    sequence = tmp[0],
                    qualities = 'F'*len(tmp[0])
                    ) 

            # R3 Output (Genomic R2)
            r3name = r2.name.split('_')[-1]
            r3seq = r2.sequence
            r3qua = r2.qualities
            if len(r3seq) > 50:
                r3_sequence = r3seq[:50]
                r3_qualities = r3qua[:50]
            else:
                r3_sequence = r3seq
                r3_qualities = r3qua
            
            re3 = dnaio.Sequence(
                    name = r3name,
                    sequence = r3_sequence,
                    qualities = r3_qualities
                    )

            fh1.write(re1)
            fh2.write(re2)
            fh3.write(re3)
            available += 1

    # Update summary json
    if os.path.exists(atacjson):
        try:
            with open(atacjson, "r") as fh:
                step1_summary = json.load(fh)
        except:
            step1_summary = {"stat": {}}
            
        if "stat" not in step1_summary:
            step1_summary["stat"] = {}

        with open(atacjson, "w") as fh1:
            step1_summary["stat"]["step1_readspair"] = total
            step1_summary["stat"]["step1_available"] = available
            
            # Helper to convert numpy types if any
            def default(o):
                if hasattr(o, 'item'): return o.item()
                raise TypeError
                
            json.dump(
                step1_summary,
                fh1,
                indent=4,
                default=default
            )
    
    return outafq1, outafq2
