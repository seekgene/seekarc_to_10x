import dnaio
import click
import os,sys
import json
from collections import defaultdict
@click.command()
@click.option('--fq1',help='fastq R1')
@click.option('--fq2',help='fastq R2')
@click.option('--outfq1',help='fastq R1 output')
@click.option('--outfq2',help='fastq R2 output')
@click.option('--out_json',help='output json')
@click.option('--out_html',help='output html')
@click.option('--params',help='other params')
@click.option('--fastp_path',help='fastp path')
def run_fastp(fq1,fq2,outfq1,outfq2,out_json,out_html,params,fastp_path):
    params = json.loads(params)
    fq3 = fq1
    fq4 = fq2
    count = 0
    length_required = params.get('length_required',60)
    f1 = dnaio.open(file1 = fq1, file2 = fq2)
    for line in f1:
        fq1 = line[0]
        fq2 = line[1]
        count += 1
        if count <= 10000:
            tmp = len(fq1.sequence)
            if tmp <= length_required: 
                length_required = tmp
            else:
                length_required = 60
        else:
            break
    params['length_required'] = length_required
    other_cmd = ' '.join([ f'--{_} {params[_]}' for _ in params.keys()])
    outdir = os.path.dirname(outfq1)
    os.makedirs(outdir,exist_ok = True)
    cmd = f"{fastp_path} -i {fq3} -I {fq4} -o {outfq1} -O {outfq2} -j {out_json} -h {out_html}  {other_cmd};"
    os.system(cmd)

if __name__ == "__main__":
    run_fastp()
