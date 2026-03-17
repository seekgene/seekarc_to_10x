import subprocess
import json
import time
import os
import argparse

# ------------------------------------------------------------------
# 命令行参数解析
# ------------------------------------------------------------------
def get_args():
    parser = argparse.ArgumentParser(description="SeekARC_to_10X Conversion Pipeline")
    
    # RNA Inputs
    parser.add_argument("--rna_r1", nargs='+', required=False, help="RNA R1 files (one or more)")
    parser.add_argument("--rna_r2", nargs='+', required=False, help="RNA R2 files (one or more)")
    
    # ATAC Inputs
    parser.add_argument("--atac_r1", nargs='+', required=False, help="ATAC R1 files (one or more)")
    parser.add_argument("--atac_r2", nargs='+', required=False, help="ATAC R2 files (one or more)")
    
    # Output and Config
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--core", type=int, default=12, help="Core number (default: 12)")
    parser.add_argument("--mode", choices=['multiome', 'atac'], default='multiome', 
                        help="Pipeline mode: 'multiome' (RNA+ATAC) or 'atac' (ATAC only). Default: multiome")
    
    args = parser.parse_args()

    # Manual validation based on mode
    if args.mode == 'multiome':
        if not (args.rna_r1 and args.rna_r2 and args.atac_r1 and args.atac_r2):
            parser.error("Mode 'multiome' requires --rna_r1, --rna_r2, --atac_r1, and --atac_r2.")
    elif args.mode == 'atac':
        if not (args.atac_r1 and args.atac_r2):
            parser.error("Mode 'atac' requires --atac_r1 and --atac_r2.")
            
    return args

# ------------------------------------------------------------------
# 辅助函数: 处理多个 Fastq 输入
# ------------------------------------------------------------------
def prepare_fastq_input(files, outdir, prefix):
    """
    检查输入文件列表。
    如果是单个文件，直接返回该文件路径。
    如果是多个文件，生成 cat 命令并返回合并后的文件路径。
    返回: (cat_cmd_string, final_file_path)
    """
    if not files:
        raise ValueError(f"No files provided for {prefix}")
        
    if len(files) == 1:
        return "", files[0]
    
    # 定义合并后的文件路径
    merged_dir = os.path.join(outdir, "merged")
    merged_path = os.path.join(merged_dir, f"{prefix}.fastq.gz")
    
    # 生成合并命令
    # 注意：确保目录存在
    files_quoted = [f'"{f}"' for f in files]
    cat_cmd = (
        f"mkdir -p {merged_dir} && "
        f"cat {' '.join(files_quoted)} > {merged_path}"
    )
    
    return cat_cmd, merged_path

# ------------------------------------------------------------------
# 主逻辑
# ------------------------------------------------------------------
def main():
    args = get_args()
    mode = args.mode
    
    # 校验输入文件数量一致性
    if mode == 'multiome':
        if not args.rna_r1 or not args.rna_r2 or not args.atac_r1 or not args.atac_r2:
             raise ValueError("Error: Mode 'multiome' requires all input types.")
        if len(args.rna_r1) != len(args.rna_r2):
            raise ValueError("Error: RNA R1 and R2 file counts do not match!")
        if len(args.atac_r1) != len(args.atac_r2):
            raise ValueError("Error: ATAC R1 and R2 file counts do not match!")
    elif mode == 'atac':
        if not args.atac_r1 or not args.atac_r2:
             raise ValueError("Error: Mode 'atac' requires ATAC input types.")
        if len(args.atac_r1) != len(args.atac_r2):
            raise ValueError("Error: ATAC R1 and R2 file counts do not match!")

    # 从 args 获取参数
    outdir = os.path.abspath(args.outdir)
    samplename = args.sample
    core = args.core
    
    # 工具与资源文件路径
    # 获取脚本所在目录
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 脚本文件
    wl_10x = os.path.join(current_dir, "lib/barcode/10x_ARC_barcode.tsv")
    wl_10x_atac = os.path.join(current_dir, "lib/barcode/10x_ATAC_barcode.tsv")
    wl_dd = os.path.join(current_dir, "lib/barcode/P3CB.barcode.txt.gz")
    fastp_py = os.path.join(current_dir, "lib/fastp/fastp.py")
    fastp_table_py = os.path.join(current_dir, "lib/fastp/fastp_table.py")
    conv_tool = os.path.join(current_dir, "bin/conv.0.1.2")
    fastp_tool = os.path.join(current_dir, "bin/fastp")
    seekarc_rna_script = os.path.join(current_dir, "lib/seekarctools_step1/seekarc_rna_step1_standalone.py")
    seekarc_atac_script = os.path.join(current_dir, "lib/seekarctools_step1/seekarc_atac_step1_standalone.py")
    
    # 验证关键工具和资源文件是否存在
    required_files = [wl_10x, wl_10x_atac, wl_dd, fastp_py, fastp_table_py, conv_tool, fastp_tool, seekarc_rna_script, seekarc_atac_script]
        
    for f in required_files:
        if not os.path.exists(f):
             print(f"Warning: Critical tool or resource file missing: {f}")
             # raise FileNotFoundError(f"Critical tool or resource file missing: {f}")
    
    # 辅助变量
    wl_rna = os.path.join(outdir, "data/rna_wl.tsv.gz")
    wl_atac = os.path.join(outdir, "data/atac_wl.tsv.gz")
    
    # ------------------------------------------------------------------
    # Step 0: 准备输入文件 (处理多 Lane 合并)
    # ------------------------------------------------------------------
    prep_cmds = []
    
    atac_fq1, atac_fq2 = "", ""
    rna_fq1, rna_fq2 = "", ""

    # ATAC (自动处理单文件/多文件：单文件返回原路径，多文件返回合并路径)
    cat_cmd_atac_r1, atac_fq1 = prepare_fastq_input(args.atac_r1, outdir, f"{samplename}_atac_S1_L001_R1_001")
    cat_cmd_atac_r2, atac_fq2 = prepare_fastq_input(args.atac_r2, outdir, f"{samplename}_atac_S1_L001_R2_001")
    if cat_cmd_atac_r1: prep_cmds.append(cat_cmd_atac_r1)
    if cat_cmd_atac_r2: prep_cmds.append(cat_cmd_atac_r2)
    
    if mode == 'multiome':
        # RNA (自动处理单文件/多文件：单文件返回原路径，多文件返回合并路径)
        cat_cmd_rna_r1, rna_fq1 = prepare_fastq_input(args.rna_r1, outdir, f"{samplename}_rna_S1_L001_R1_001")
        cat_cmd_rna_r2, rna_fq2 = prepare_fastq_input(args.rna_r2, outdir, f"{samplename}_rna_S1_L001_R2_001")
        if cat_cmd_rna_r1: prep_cmds.append(cat_cmd_rna_r1)
        if cat_cmd_rna_r2: prep_cmds.append(cat_cmd_rna_r2)
    
    # 如果有合并任务，先执行
    print("Step 0: Checking input files ...", flush=True)
    step0_cmd = ""
    if prep_cmds:
        step0_cmd = " && ".join(prep_cmds)
        print(f"Executing Step 0 Merge: {step0_cmd}", flush=True)
        result0 = subprocess.run(step0_cmd, shell=True, executable='/bin/bash')
        if result0.returncode != 0:
            print(f"Step 0 Failed, exit code: {result0.returncode}", flush=True)
            raise RuntimeError(f"Step 0 (Input Preparation) failed with exit code {result0.returncode}")
        print("Step 0 Merge Success.", flush=True)
    else:
        print("Step 0: Single input files detected, skipping merge.", flush=True)
    
    # ------------------------------------------------------------------
    # 第一步：Fastp 任务
    # ------------------------------------------------------------------
    
    # 构建 Fastp 命令列表
    fastp_cmds_list = [f"mkdir -p {outdir}/fastp"]
    
    # --- ATAC Fastp ---
    fastp_cmds_list.append(
        f"python {fastp_py} "
        f"--fq1 {atac_fq1} "
        f"--fq2 {atac_fq2} "
        f"--outfq1 {outdir}/fastp/{samplename}_atac_S1_L001_R1_001.fastq.gz "
        f"--outfq2 {outdir}/fastp/{samplename}_atac_S1_L001_R2_001.fastq.gz "
        f"--out_json {outdir}/fastp/{samplename}_atac_fastp.json "
        f"--out_html {outdir}/fastp/{samplename}_atac_fastp.html "
        f"--params '{{\"cut_tail_window_size\": 1, \"cut_tail_mean_quality\": 3, \"cut_tail\": \"\"}}' --fastp_path {fastp_tool} && "
        f"python {fastp_table_py} {outdir}/fastp/{samplename}_atac_fastp.json {outdir}/fastp/ {samplename}_atac"
    )

    if mode == 'multiome':
        # --- RNA Fastp ---
        fastp_cmds_list.append(
            f"python {fastp_py} "
            f"--fq1 {rna_fq1} "
            f"--fq2 {rna_fq2} "
            f"--outfq1 {outdir}/fastp/{samplename}_rna_S1_L001_R1_001.fastq.gz "
            f"--outfq2 {outdir}/fastp/{samplename}_rna_S1_L001_R2_001.fastq.gz "
            f"--out_json {outdir}/fastp/{samplename}_rna_fastp.json "
            f"--out_html {outdir}/fastp/{samplename}_rna_fastp.html "
            f"--params '{{\"cut_tail_window_size\": 1, \"cut_tail_mean_quality\": 3, \"cut_tail\": \"\"}}' --fastp_path {fastp_tool} && "
            f"python {fastp_table_py} {outdir}/fastp/{samplename}_rna_fastp.json {outdir}/fastp/ {samplename}_rna"
        )
    
    fastcmd = " && ".join(fastp_cmds_list)
    
    print("Step 1: Running Fastp ...", flush=True)
    result1 = subprocess.run(fastcmd, shell=True, executable='/bin/bash')
    
    if result1.returncode != 0:
        print(f"Step 1 Failed, exit code: {result1.returncode}", flush=True)
        raise RuntimeError(f"Step 1 (Fastp) failed with exit code {result1.returncode}")

    print("Step 1 Success. Start Step 2: Conversion ...", flush=True)

    # ------------------------------------------------------------------
    # 第二步：数据转换与预处理任务 (ATAC + RNA or ATAC only)
    # ------------------------------------------------------------------
    
    conversion_cmds = []

    # 1. 提取白名单
    cmd_wl = f"mkdir -p {outdir}/data/ && cd {outdir}/data/"
    if mode == 'multiome':
         cmd_wl += f" && awk '{{print $1}}' {wl_10x} | gzip > {wl_rna} && awk '{{print $2}}' {wl_10x} | gzip > {wl_atac}"
    else: # atac only
         cmd_wl += f" && cat {wl_10x_atac} | gzip > {wl_atac}"
    conversion_cmds.append(cmd_wl)
    
    if mode == 'multiome':
        # 2. 转化 RNA (生成 RNA map)
        cmd_rna_conv = (
            f"mkdir -p {outdir}/data/E/to10x && "
            f"{conv_tool} --fq1 {outdir}/fastp/{samplename}_rna_S1_L001_R1_001.fastq.gz --fq2 {outdir}/fastp/{samplename}_rna_S1_L001_R2_001.fastq.gz "
            f"--wl1 {wl_dd} --wl2 {wl_rna} --rs 17C+T -t {core} -o {outdir}/data/E/to10x"
        )
        conversion_cmds.append(cmd_rna_conv)
        
        # 3. 处理 Mapfile
        cmd_map = (
            f"cp {outdir}/data/E/to10x/map.tsv {outdir}/data/dd_rna_map.tsv && "
            f"""awk -F'\\t' 'NR==FNR{{a[$2]=$1;next}} {{if($1 in a)print a[$1]"\\t"$1"\\t"$2}}' {outdir}/data/dd_rna_map.tsv {wl_10x} > {outdir}/data/dd_rna_atac_map.tsv && """
            f"cut -f1,3 {outdir}/data/dd_rna_atac_map.tsv > {outdir}/data/dd_atac_map.tsv"
        )
        conversion_cmds.append(cmd_map)
        
        # 4. 转化 ATAC (使用 map)
        cmd_atac_conv = (
            f"mkdir -p {outdir}/data/A/to10x && "
            f"{conv_tool} --fq1 {outdir}/fastp/{samplename}_atac_S1_L001_R1_001.fastq.gz --fq2 {outdir}/fastp/{samplename}_atac_S1_L001_R2_001.fastq.gz "
            f"--map {outdir}/data/dd_atac_map.tsv --rs 17C+T -t {core} -o {outdir}/data/A/to10x"
        )
        conversion_cmds.append(cmd_atac_conv)
    else:
        # ATAC only 转化 (直接使用白名单)
        cmd_atac_conv = (
             f"mkdir -p {outdir}/data/A/to10x && "
             f"{conv_tool} --fq1 {outdir}/fastp/{samplename}_atac_S1_L001_R1_001.fastq.gz "
             f"--fq2 {outdir}/fastp/{samplename}_atac_S1_L001_R2_001.fastq.gz "
             f"--wl1 {wl_dd} --wl2 {wl_atac} --rs 17C+T -t {core} -o {outdir}/data/A/to10x && "
             f"mv {outdir}/data/A/to10x/map.tsv {outdir}/data/map.tsv"
        )
        conversion_cmds.append(cmd_atac_conv)
    
    # 5. Step1 去除 TSO (RNA + ATAC)
    if mode == 'multiome':
        cmd_step5_rna = (
            f"python {seekarc_rna_script} "
            f"--fq1 {outdir}/data/E/to10x/{samplename}_rna_S1_L001_R1_001.fastq.gz "
            f"--fq2 {outdir}/data/E/to10x/{samplename}_rna_S1_L001_R2_001.fastq.gz "
            f"--samplename {samplename} --chemistry custom --structure B16U12 --barcode {wl_rna} "
            f"--outdir {outdir}/data/E/ --core {core} --use_cr"
        )
        conversion_cmds.append(cmd_step5_rna)

    cmd_step5_atac = (
        f"python {seekarc_atac_script} "
        f"--fq1 {outdir}/data/A/to10x/{samplename}_atac_S1_L001_R1_001.fastq.gz "
        f"--fq2 {outdir}/data/A/to10x/{samplename}_atac_S1_L001_R2_001.fastq.gz "
        f"--samplename {samplename} --chemistry custom --structure B16U12 --barcode {wl_atac} "
        f"--outdir {outdir}/data/A/ --core {core} "
    )
    conversion_cmds.append(cmd_step5_atac)
    
    # 6. 改名 (RNA + ATAC)
    if mode == 'multiome':
        cmd_step6_rna = (
            f"mv {outdir}/data/E/step1/{samplename}_E_1.fq.gz {outdir}/data/E/step1/{samplename}_S1_L001_R1_001.fastq.gz && "
            f"mv {outdir}/data/E/step1/{samplename}_E_2.fq.gz {outdir}/data/E/step1/{samplename}_S1_L001_R2_001.fastq.gz "
        )
        conversion_cmds.append(cmd_step6_rna)

    cmd_step6_atac = (
        f"mv {outdir}/data/A/step1/{samplename}_A_cutR1.fastq.gz {outdir}/data/A/step1/{samplename}_S1_L001_R1_001.fastq.gz && "
        f"mv {outdir}/data/A/step1/{samplename}_A_cutR2.fastq.gz {outdir}/data/A/step1/{samplename}_S1_L001_R2_001.fastq.gz && "
        f"mv {outdir}/data/A/step1/{samplename}_A_cutR3.fastq.gz {outdir}/data/A/step1/{samplename}_S1_L001_R3_001.fastq.gz"
    )
    conversion_cmds.append(cmd_step6_atac)
    
    # 7. CellRanger
    if mode == 'multiome':
        # 生成 library.csv
        cmd_step7 = (
            f"echo 'fastqs,sample,library_type' > {outdir}/library.csv && "
            f"echo '{outdir}/data/E/step1,{samplename},Gene Expression' >> {outdir}/library.csv && "
            f"echo '{outdir}/data/A/step1,{samplename},Chromatin Accessibility' >> {outdir}/library.csv"
        )
        conversion_cmds.append(cmd_step7)
        
    else: # atac
        pass
    
    # 组合所有步骤命令
    wlconvertcmd = " && \\\n".join(conversion_cmds)
    
    # 执行流水线
    result2 = subprocess.run(wlconvertcmd, shell=True, executable='/bin/bash')
    
    if result2.returncode == 0:
        print("Pipeline Completed Successfully.", flush=True)
    else:
        print(f"Step 2 Failed, exit code: {result2.returncode}", flush=True)
        raise RuntimeError(f"Step 2 (Conversion) failed with exit code {result2.returncode}")

if __name__ == "__main__":
    main()
