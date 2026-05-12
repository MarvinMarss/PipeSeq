# PipeSeq

PipeSeq is a local Windows/WSL desktop workflow for RNA-seq processing. It provides a PyQt6 GUI
around SRA/FASTQ preparation, HISAT2 alignment, SAM/BAM processing, StringTie expression analysis,
featureCounts/DESeq2-based statistics and result visualization.

> Russian documentation: [README_RUS.md](README_RUS.md)

## What is included in git

This repository keeps the source code, launcher, documentation, examples and small static assets.
Large runtime files are deliberately not committed:

- raw FASTQ/FQ/SRA files;
- SAM/BAM outputs;
- genome FASTA files, GTF/GFF annotations and HISAT2 indexes;
- generated results and logs;
- local `Script/settings.json` / `Script/Mind.json`;
- downloaded FastQC/SRA Toolkit binaries.

The local folder layout is still represented by README placeholders in `Fastq/`, `Sra/`, `Genome/`,
`Genome/Index/`, `Output/`, `GTF/`, `Results/` and `Tool/`.

## Repository layout

```text
PipeSeq/
├── PipeSeq.bat                 # Windows launcher
├── Script/
│   ├── PipeSeq.py              # Main PyQt6 GUI
│   ├── align_hisat2.py         # FASTQ -> SAM with HISAT2
│   ├── process_sam_to_bam.py   # SAM -> sorted BAM
│   ├── stringtie_expression.py # StringTie expression step
│   ├── deseq2_analysis.py      # featureCounts + PyDESeq2 workflow
│   ├── extract_*.py            # Expression table extraction helpers
│   ├── pvalues_*.py            # Statistics helpers
│   ├── temp_card_p.py          # Heatmap/visualization helper
│   ├── settings.example.json   # Example pipeline configuration
│   └── Mind.example.json       # Example GUI memory file
├── Fastq/                      # Local FASTQ input, ignored by git
├── Sra/                        # Local SRA input/downloads, ignored by git
├── Genome/                     # Local genome/annotation files, ignored by git
├── Output/                     # Local SAM/BAM output, ignored by git
├── GTF/                        # Local StringTie output, ignored by git
├── Results/                    # Local final results, ignored by git
├── Tool/                       # Optional local third-party tools, ignored by git
└── docs/                       # Maintainer notes
```

## Requirements

- Windows with PowerShell.
- Python 3.10 or newer.
- WSL with Ubuntu for Linux command-line bioinformatics tools.
- Python packages from `requirements.txt`.
- WSL packages/tools: HISAT2, SAMtools, StringTie and Subread/featureCounts.
- Optional local Windows tools in `Tool/`: FastQC and SRA Toolkit.

## Installation

Open PowerShell as Administrator and install Python:

```powershell
winget install -e --id Python.Python.3.12
```

Create a virtual environment from the repository root:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Install Ubuntu in WSL:

```powershell
wsl --install -d Ubuntu
```

Inside Ubuntu, install the required command-line tools:

```bash
sudo apt update
sudo apt install -y git unzip wget curl hisat2 samtools stringtie subread
```

## Local data setup

1. Place FASTQ files in `Fastq/`, or place local SRA files in `Sra/`.
2. Place the reference genome FASTA and matching annotation GTF/GFF in `Genome/`.
3. Place a prebuilt HISAT2 index in `Genome/Index/`, or let the pipeline create it when supported.
4. Put optional local Windows tool bundles in `Tool/`, or install them elsewhere and configure paths
   in the GUI.

On first launch, PipeSeq creates `Script/settings.json` if it does not exist. You can also copy
`Script/settings.example.json` to `Script/settings.json` and edit the paths manually.

## Running

On Windows, double-click:

```text
PipeSeq.bat
```

Or run the GUI directly:

```powershell
python Script\PipeSeq.py
```

The launcher changes into `Script/` before starting the GUI because the pipeline scripts expect
their configuration files in that folder.

## Experiment IDs

Input experiment IDs and sample names in the GUI in this format:

```text
SRX8380271-HighLight1; SRX8380270-HighLight2; SRX8380269-HighLight3; SRX5120532-HighLightControl1; SRX5120531-HighLightControl2; SRX5120530-HighLightControl3
```

Replicate numbers should be placed at the end of the sample name. Local `.sra` files can also be
renamed using the same naming convention.

## Development checks

```powershell
python -m py_compile Script\*.py
```

Optional linting:

```powershell
pip install -r requirements-dev.txt
ruff check Script
```

## Notes

- HISAT2 index files usually use the basename `genome_index` and produce eight `.ht2` files.
- Windows paths are converted to WSL paths internally where required.
- See [docs/DATA_AND_TOOLS.md](docs/DATA_AND_TOOLS.md) for the repository data policy.

## License

MIT. See [LICENSE](LICENSE).
