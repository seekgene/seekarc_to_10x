# -*- coding: utf-8 -*-
import os

# ==========================================
# Configuration & Constants
# ==========================================
__version__ = "1.1.0"

# Constants
R1_MINLEN = 20
R2_MINLEN = 60
ATAC_R1_MINLEN = 10
ATAC_R2_MINLEN = 10

# Directory of this script, used to locate resources if present
# When running as standalone, resources might be relative to this file
# or provided via arguments. Here we assume relative to a 'barcode' directory if it exists.
__srcdir = os.path.dirname(os.path.abspath(__file__))

# Chemistry configurations
CHEMISTRY = {
    '__SO01V3':{
        'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8L8B8L10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'linker': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'Linker1.txt'),
                   os.path.join(__srcdir, 'barcode', 'SO01V3', 'Linker2.txt'),),
    },
    '__nolinker':{
        'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
    },
    "__P3CBGB":{
        'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
    },
    "DDV1":{
       'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'match_type': (1,),
        'adapter1': [["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"], ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ["TTTTTTTTTTTT", "5"]],
        'adapter2': [["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ["AAAAAAAAAAAA", "3"], ["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"]],
    },
    "DDV2":{
        'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        'match_type': (1,),
        'adapter1': [["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"], ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ["TTTTTTTTTTTT", "5"] ], 
        'adapter2': [["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ["AAAAAAAAAAAA", "3"], ["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"]],
    },
    "DDVS":{
        'shift': False,
        'structure': 'B17U12X24B10',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),
                    os.path.join(__srcdir, 'barcode', 'P3CBGB', 'sample.barcode.txt.gz'),),
        'match_type': (1,3,),
        'adapter1': [["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"], ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ["TTTTTTTTTTTT", "5"] ],
        'adapter2': [["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"], ["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ["AAAAAAAAAAAA", "3"]],
    },
    "DD5V1":{
       'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        'match_type': (1,),
        'sc5p': True,
        'adapter1': [["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"] ], 
        'adapter2': [["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ["TTTTTTTTTTTT", "5"], ["CCCATATAAGAAA", "3"], ["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"] ],
    },
    "MM":{
        'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'match_type': (1,),
        'adapter2': [["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ["AAAAAAAAAAAA", "3"], ["AAAAAAAAAAAA", "5"], ["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"]], 
        'adapter1': [["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"],["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"], ["TTTTTTTTTTTT", "5"]], 
    },
    "MM-D":{
        'shift': True,
        'shift_pattern': 'A',
        'structure': 'B8X8B8X10B8U8',
        'barcode': (os.path.join(__srcdir, 'barcode', 'SO01V3', 'seekgene.txt'),),
        'match_type': (1,),
    },
    "DD-Q":{
        'shift': False,
        'structure': 'B17U12X7',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        'adapter1': [["GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "3"], ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ], 
        'adapter2': [["AAGCAGTGGTATCAACGCAGAGTACATGG", "5"], ["GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC", "5"], ], 
        'sc5p': None,
        'match_type': (1,),
    },
    "DD_AG":{
        'shift': False,
        'structure': 'B17U12',
        'barcode': (os.path.join(__srcdir, 'barcode', 'P3CBGB', 'P3CB.barcode.txt.gz'),),
        'adapter1': [["CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", "3"], ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "5"], ["TTTTTTTTTTTT", "5"] ], 
        'adapter2': [["GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG", "5"], ["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "3"], ["CTGTCTCTTATACACATCTACGAGCAACGACGGACG", "3"]], 
        'match_type': (1,),
    }
}
