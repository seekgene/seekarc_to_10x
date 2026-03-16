# -*- coding: utf-8 -*-
import os
import sys
import re
import json
import logging
from subprocess import run
from collections import defaultdict
from xopen import xopen
from cutadapt.adapters import BackAdapter, RightmostFrontAdapter
from seekarc_config import __version__, __srcdir

# Logger setup
logger = logging.getLogger("seekarc_step1")
if not logger.handlers:
    handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
logger.setLevel(logging.INFO)

def check_path(path):
    if not os.path.exists(path):
        logger.info(f"Error : The path of '{path}' is not exists")
        sys.exit(1)
    else:
        return path

def parse_structure(string:str) -> tuple:
    regex = re.compile(r'([BLUXT])(\d+)')
    groups = regex.findall(string)
    return tuple([(_[0], int(_[1])) for _ in groups])

def read_file(file_list: list) -> dict:
    wl_dict = dict()
    if not file_list:
        return wl_dict
    for i, wl_file in enumerate(file_list):
        white_list = set()
        if not os.path.exists(wl_file):
            logger.warning(f"Barcode file not found: {wl_file}")
            basename = os.path.basename(wl_file)
            possible_path = os.path.join(__srcdir, 'barcode', os.path.basename(os.path.dirname(wl_file)), basename)
            if os.path.exists(possible_path):
                 wl_file = possible_path
                 logger.info(f"Found barcode file at: {wl_file}")
            else:
                 if os.path.exists(basename):
                     wl_file = basename
                     logger.info(f"Found barcode file at: {wl_file}")

        try:
            with xopen(wl_file, "r") as fh:
                for l in fh:
                    if l.startswith("#"): continue
                    la = l.strip()
                    if not la: continue
                    white_list.add(la)
        except Exception as e:
            logger.error(f"Failed to read barcode file {wl_file}: {e}")
                 
        wl_dict[i] = white_list
    return wl_dict

def get_new_bc(bc:str, white_list:set, distance:int=1)->set:
    if distance == 1:
        BASE_LIST = ["T", "C", "G", "A"]
        mm_dict = dict()
        for i, c in enumerate(bc):
            if c == "N":
                mm_dict = { bc[:i] + base + bc[i+1:]:f"{i}{base}" for base in BASE_LIST }
                break  
            else:
                mm_dict.update({ bc[:i] + base + bc[i+1:]:f"{i}{base}" for base in BASE_LIST if base!=c })
        bc_set = set(mm_dict.keys()).intersection(white_list)
    else:
        bc_dict = defaultdict(set)
        for bc_true in white_list:
            hmm = sum(ch1 != ch2 for ch1,ch2 in zip(bc_true,bc))
            if hmm <= distance:
                bc_dict[hmm].add(bc_true)
        bc_set = set()
        if len(bc_dict) != 0:
            sorted_items = sorted(bc_dict.items(), key=lambda x: x[0])
            bc_set = sorted_items[0][1]
    return bc_set

class AdapterFilter:
    def __init__(self, adapter1:list=[], adapter2:list=[],):
        self.adapter1 = [BackAdapter(sequence=_, min_overlap=10, read_wildcards=True) if p=="3" else RightmostFrontAdapter(sequence=_, min_overlap=5, read_wildcards=True) for _, p in adapter1]
        self.adapter2 = [BackAdapter(sequence=_, min_overlap=10, read_wildcards=True) if p=="3" else RightmostFrontAdapter(sequence=_, min_overlap=10, read_wildcards=True) for _, p in adapter2]
    
    def filter(self, r1=None, r2=None) -> tuple:
        flag = False
        if r1 and self.adapter1:
            for _ in self.adapter1:
                m = _.match_to(r1.sequence)
                if m:
                    flag = True
                    r1 =  m.trimmed(r1)
        if r2 and self.adapter2:
            for _ in self.adapter2:
                m = _.match_to(r2.sequence)
                if m:
                    flag = True
                    r2 =  m.trimmed(r2)
        return flag, r1, r2

class QcStat:
    def __init__(self):
        self.data = { }
    def update(self, **d):
        if not self.data:
            self.data = d
        else:
            for k, v in d.items():
                self.data[k] += v
    @staticmethod
    def _sort_gc(d):
        if not d: return {}
        idx_max = max([k[0] for k in d])
        return {b: [d.get((i, b), 0) for i in range(idx_max+1)] for b in 'ATCGN'}
    @staticmethod
    def _sort_q(d, phred=33):
        if not d: return {}
        idx_max = max([k[0] for k in d])
        q_max = max([ord(k[1])-phred for k in d])
        return {i: [d.get((i, chr(q+phred)), 0) for q in range(q_max+1)] for i in range(idx_max+1)}
    def save(self, path='summary.json'):
        tmp = {'__version__': __version__}
        for k in self.data:
            if k.endswith('_gc'):
                if not self.data[k]:
                    tmp[k] = {}
                    continue
                else:
                    tmp[k] = self._sort_gc(self.data[k])
            elif k.endswith('_q'):
                tmp[k] = self._sort_q(self.data[k])
            else:
                tmp[k] = dict(self.data[k])
        with open(path, 'w') as fh:
            json.dump(tmp, fh, indent=2)

def cmd_execute(
    args, check:bool=False, text:bool=True,
    capture_output:bool=True, env:os.environ=None
):
    if isinstance(args, list):
        args = [str(_) for _ in args]
        logger.info(" ".join(args))
        _call = run(args, check=False, text=text, capture_output=capture_output, env=env)
    elif isinstance(args, str):
        logger.info(args)
        _call = run(args, shell=True, check=False, text=text, capture_output=capture_output, env=env)
    if check:
        if _call.returncode!=0:
            logger.info(f"stderr: {_call.stderr}")
            logger.info(f"stdout: {_call.stdout}")
            sys.exit(1)
    return _call
