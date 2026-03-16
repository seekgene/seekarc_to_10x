# SeekARC_to_10X: SeekARC 测序数据转 10x Genomics 转换流程

本流程专为解决 SeekARC 测序平台与 10x Genomics 下游分析生态系统之间的数据兼容性问题而设计。通过自动化的数据转换管道，流程能够高效处理 SeekARC 双组学（Multiome RNA+ATAC）及单组学（Single ATAC）原始数据。核心模块集成了多样本数据聚合、基于 Fastp 的高精度质量控制、以及严格的条形码映射与格式重构算法，确保输出数据在文件结构、序列格式及元数据层面完全符合 Cell Ranger ARC 和 Cell Ranger ATAC 的输入标准，从而无缝衔接主流的单细胞多组学分析工作流。

## 概述

SeekARC 平台产出的数据需要特定的预处理才能兼容广泛使用的 10x Genomics 分析生态系统。本流程通过以下步骤弥合了这一差异：

1.  **预处理 (Preprocessing)**：自动处理多文件输入（例如多 Lane 数据）并进行合并。
2.  **质量控制 (Quality Control)**：使用 `fastp` 去除接头并过滤低质量 Reads。
3.  **格式转换 (Format Conversion)**：
    *   将 SeekARC 条形码映射到 10x Genomics 白名单条形码。
    *   重构 Fastq 文件以符合 Cell Ranger 的输入规范。
    *   (Multiome 模式) 去除 RNA Reads 中的 TSO (Template Switch Oligo) 序列。

## 功能特性

*   **双模式支持 (Dual Mode Support)**：
    *   **Multiome 模式 (`--mode multiome`)**：处理配对的 RNA-seq 和 ATAC-seq 数据。
    *   **Single ATAC 模式 (`--mode atac`)**：仅处理 ATAC-seq 数据。
*   **强大的输入处理 (Robust Input Handling)**：支持复杂的输入场景，例如数据分布在多个 Lane 或文件中。只需提供所有文件路径，流程会自动处理合并。
*   **集成质量控制 (Integrated QC)**：为每个样本生成详细的 `fastp` 质控报告 (HTML/JSON)。

## 安装与环境配置

本流程依赖于特定的 Python 环境和若干外部工具。

### 1. Conda 环境

建议使用 Conda 管理依赖。项目根目录下提供了 `env.yaml` 文件用于创建环境。

```bash
# 创建环境
conda env create -f env.yaml

# 激活环境
conda activate seekarc_env
```

**核心依赖库：**
*   Python >= 3.10
*   cutadapt
*   click
*   dnaio
*   xopen
*   pigz

### 2. 外部工具

请确保以下工具已安装，并且在 `pipeline.py` 中正确配置了路径（或已添加至系统 `$PATH`）：

*   **Fastp**：已包含在项目的 `bin/` 目录下。
*   **conv.0.1.1**：已包含在项目的 `bin/` 目录下。

## 快速开始 (Quick Start)

项目提供了测试数据和一键运行脚本，方便您快速验证环境和流程。

*   **测试数据**：位于 `data/` 目录下。
*   **运行脚本**：`run.sh`。

```bash
# 确保已激活环境
conda activate seekarc_env

# 运行测试脚本
bash run.sh
```

该脚本会使用 `data/` 目录下的测试数据运行 Multiome 模式，并将结果输出到 `test_output/` 目录。

## 使用方法

使用 `python` 命令运行流程脚本。完整参数列表如下：

```bash
python pipeline.py \
    --mode <multiome|atac> \
    --outdir <输出目录> \
    --sample <样本名称> \
    --rna_r1 <RNA_R1_文件> [RNA_R1_文件 ...] \
    --rna_r2 <RNA_R2_文件> [RNA_R2_文件 ...] \
    --atac_r1 <ATAC_R1_文件> [ATAC_R1_文件 ...] \
    --atac_r2 <ATAC_R2_文件> [ATAC_R2_文件 ...] \
    --core <线程数>
```

### 参数说明

#### 必选参数

| 参数 | 说明 |
| :--- | :--- |
| `--outdir` | 结果输出的主目录，所有 QC、中间文件和 转换后的测序文件都存放在此。 |
| `--sample` | 样本的唯一标识符。该名称将用于命名输出文件。 |

#### 流程模式

| 参数 | 默认值 | 说明 |
| :--- | :--- | :--- |
| `--mode` | `multiome` | 指定分析模式：`multiome` (RNA+ATAC) 或 `atac` (仅 ATAC)。 |

#### 输入文件

对于每个输入参数，您可以提供多个文件（用空格分隔）。流程会自动合并它们。

**Multiome 模式 (`--mode multiome`) 需提供：**
*   `--rna_r1`: RNA Read 1 文件路径。
*   `--rna_r2`: RNA Read 2 文件路径。
*   `--atac_r1`: ATAC Read 1 文件路径。
*   `--atac_r2`: ATAC Read 2 文件路径。

**Single ATAC 模式 (`--mode atac`) 需提供：**
*   `--atac_r1`: ATAC Read 1 文件路径。
*   `--atac_r2`: ATAC Read 2 文件路径。
*   *(RNA 参数在此模式下将被忽略)*

#### 资源选项

| 参数 | 默认值 | 说明 |
| :--- | :--- | :--- |
| `--core` | `12` | 用于并行处理 (Fastp, 转换) 的 CPU 线程数。 |

## 运行示例

### 1. Multiome 分析 (多 Lane 输入)

```bash
python pipeline.py \
    --mode multiome \
    --outdir ./results/sample_multiome \
    --sample MySample \
    --rna_r1 ./raw/RNA_L1_R1.fq.gz ./raw/RNA_L2_R1.fq.gz \
    --rna_r2 ./raw/RNA_L1_R2.fq.gz ./raw/RNA_L2_R2.fq.gz \
    --atac_r1 ./raw/ATAC_L1_R1.fq.gz ./raw/ATAC_L2_R1.fq.gz \
    --atac_r2 ./raw/ATAC_L1_R2.fq.gz ./raw/ATAC_L2_R2.fq.gz \
    --core 24
```

### 2. Single ATAC Analysis

```bash
python pipeline.py \
    --mode atac \
    --outdir ./results/sample_atac \
    --sample MySample \
    --atac_r1 ./raw/ATAC_S1_R1.fq.gz \
    --atac_r2 ./raw/ATAC_S1_R2.fq.gz \
    --core 16
```

## 流程步骤详解

### 第一步：输入准备 (Input Preparation)
流程首先验证所有输入文件。如果为同一种 Read 类型提供了多个文件（例如 `L001_R1.fq.gz`, `L002_R1.fq.gz`），它们会被串联（concatenated）成一个文件，存放在 `{outdir}/merged/` 目录下。

### 第二步：质量控制 (Quality Control)
原始 Reads 使用 `fastp` 进行处理。
*   **修剪 (Trimming)**：去除接头序列和低质量碱基（尾部修剪）。
*   **报告 (Reporting)**：在 `{outdir}/fastp/` 生成 HTML 和 JSON 报告供检查。
*   **输出**：清洗后的 Clean Fastq 文件保存在 `{outdir}/fastp/`，用于下一步处理。

### 第三步：条形码翻译与格式转换 (Barcode Translation & Conversion)
这是核心转换步骤。
*   **白名单生成**：从提供的 10x Whitelist 文件中提取有效 Barcode。
*   **映射 (Mapping)**：
    *   **Multiome**：创建映射文件 (`dd_rna_atac_map.tsv`) 以关联 SeekARC Barcode 和 10x Barcode。由于 10x RNA 和 ATAC 文库的 Barcode 不同，该文件包含三列：SeekARC Barcode、10x RNA Barcode、10x ATAC Barcode。
    *   **ATAC**：使用 ATAC 白名单进行直接映射。在 `{outdir}/data/` 生成 `map.tsv` 文件。
*   **文件转换**：使用 `conv.0.1.1` 工具和自定义 Python 脚本进行以下处理，使其完全兼容 Cell Ranger：
    *   **条形码矫正 (Barcode Correction)**：基于 10x 白名单修正测序错误（允许 1 bp 错配）。
    *   **序列修剪 (Trimming)**：去除接头序列及 TSO (Template Switch Oligo) 序列（仅 Multiome RNA）。
    *   **格式重构 (Reformatting)**：转换 Fastq Header 和序列结构。

### 第四步：Cell Ranger 分析 (Analysis)
流程会自动生成所需的输入文件（Multiome 模式生成 `library.csv`，ATAC 模式生成转换后的 Fastq），但 **不会** 自动执行 Cell Ranger。请根据“运行 Cell Ranger”部分的说明手动运行。

## 输出目录结构

输出目录结构组织如下：

```text
{outdir}/
├── data/                         # 中间转换数据
│   ├── A/                        # ATAC 转换文件
│   ├── E/                        # RNA 转换文件 (仅 Multiome)
│   ├── map.tsv                   # Barcode 映射文件 (仅 ATAC 模式)
│   ├── dd_rna_atac_map.tsv     # Multiome Barcode 映射文件
│   └── ...
├── fastp/                        # 质量控制结果
│   ├── {sample}_atac_fastp.html  # ATAC 数据的 QC 报告
│   ├── {sample}_rna_fastp.html   # RNA 数据的 QC 报告
│   └── ...
├── merged/                       # 合并后的输入文件
│   └── ... (仅当提供了多个输入文件时存在)
└── library.csv                   # Cell Ranger Library CSV (仅 Multiome)
```

## 运行 Cell Ranger

流程成功完成后，您可以根据以下示例手动运行 Cell Ranger。

### 1. Multiome (Cell Ranger ARC)

Multiome 模式下，流程会在输出目录生成 `library.csv` 文件。

```bash
# 运行示例
/path/to/cellranger-arc count \
    --id=<样本名称> \
    --reference=/path/to/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
    --libraries=<输出目录>/library.csv \
    --localcores=12 \
    --localmem=64
```

### 2. Single ATAC (Cell Ranger ATAC)

Single ATAC 模式下，请使用位于 `<输出目录>/data/A/step1/` 的转换后 Fastq 文件。

```bash
# 运行示例
/path/to/cellranger-atac count \
    --id=<样本名称> \
    --reference=/path/to/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
    --fastqs=<输出目录>/data/A/step1/ \
    --chemistry=ARC-v1 \
    --sample=<样本名称> \
    --localcores=12 \
    --localmem=64
```
