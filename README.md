Environment Setup on a New PC
Before first use, all dependencies must be installed.

1. Open PowerShell as Administrator
Press Win + X → select Windows PowerShell (Admin) or Terminal (Admin).

2. Install Python and WSL (Ubuntu)
In PowerShell, run:

powershell

winget install -e --id Python.Python.3.12
Restart the terminal, then install Python packages:

bash

pip install PyQt6
pip install pandas numpy scipy seaborn matplotlib pyDESeq2
pip install pywin32
Install WSL with Ubuntu:

powershell

wsl --install -d Ubuntu
Reboot the PC if prompted.

3. After reboot — open PowerShell (Admin):

powershell

wsl --install -d Ubuntu
wsl
Ubuntu will launch. Follow the prompts to set a username and password (any memorable values).

4. Inside Ubuntu — set up the environment:

bash

# System update
sudo apt update && sudo apt upgrade -y
sudo add-apt-repository universe
sudo add-apt-repository multiverse
sudo apt update

# Install required packages
sudo apt install -y python3-venv build-essential zlib1g-dev \
  libbz2-dev liblzma-dev libncurses-dev \
  libcurl4-openssl-dev libssl-dev libsqlite3-dev wget curl \
  git unzip samtools hisat2 stringtie libgl1 libxkbcommon-x11-0

# Create Python virtual environment
python3 -m venv ~/pipeseq_env
source ~/pipeseq_env/bin/activate

# Upgrade pip and install Python packages
pip install --upgrade pip
pip install pandas numpy scipy seaborn matplotlib pyqt6 pyDESeq2
5. Running the Pipeline
Double-click PipeSeq.bat to launch the pipeline GUI.

6. On first launch, configure folder paths inside PipeSeq.
7. Place the appropriate .fa genome and corresponding .gtf annotation files in the "Genome" folder. If annotation errors occur, enable gtf.fix in PipeSeq.
8. To begin analysis, input experiment IDs and assigned names in this format:


SRX8380271-HighLight1; SRX8380270-HighLight2; SRX8380269-HighLight3; SRX5120532-HighLightControl1; SRX5120531-HighLightControl2; SRX5120530-HighLightControl3
Replicate numbers go at the end. You may also use local .sra files renamed accordingly.

Project Structure
bash

PipeSeq/
├── align_hisat2.py             
├── process_sam_to_bam.py       
├── stringtie_expression.py     
├── extract_fpkm.py             
├── pvalues_log2.py             
├── GTF_results_pvalues.py      
├── temp_card_p.py              
├── fix.gtf.py                  
├── run_gui.py                  
├── settings.json               
└── pipeline_log.txt            
settings.json Configuration
json

{
  "folders": {
    "fastq_folder": "path to FASTQ files",
    "bam_folder": "path for BAM output",
    "gtf_folder": "path to StringTie GTF files",
    "results_folder": "path to final tables",
    "genome_folder": "path to .fa and .gtf files",
    "genome_index": "path to HISAT2 index"
  },
  "options": {
    "delete_intermediate_files": true,
    "use_fdr_correction": true,
    "fix_genome": false
  },
  "stringtie": {
    "coverage_cutoff": 0.01
  },
  "gene_mapping": {
    "CHLRE_01g025050v5": "GATA-1",
    "CHLRE_10g435450v5": "GATA-2"
  },
  "visualization": {
    "show_p_values": true
  }
}
Running the Pipeline
bash

python run_gui.py
Or step-by-step:

bash

python align_hisat2.py
python process_sam_to_bam.py
python stringtie_expression.py
python extract_fpkm.py
python pvalues_log2.py
Features
Supports paired-end and single-end FASTQ.

Intelligent processing: sorting, skipping steps, temp file cleanup.

Genome index is auto-generated if missing.

StringTie sensitivity (-c) configurable via GUI.

Transparent logs and visualization.

Dependencies
Python 3.10+

PyQt6

HISAT2 (via WSL)

SAMtools (via WSL)

StringTie (via WSL)

Notes
HISAT2 requires 8 .ht2 index files with base name genome_index.

All WSL paths are auto-converted (e.g., /mnt/c/...).

If needed, use fix.gtf.py to correct annotation errors.

Uninstalling the Environment
1. Open PowerShell as Administrator

2. Remove WSL and Python:

powershell

winget uninstall Python.Python.3.12
pip uninstall -y PyQt6 pandas numpy scipy seaborn matplotlib pyDESeq2
wsl --uninstall
Wait for uninstallation to complete.

3. Confirm Ubuntu is fully removed:

powershell

wsl --unregister Ubuntu
4. If using a Python virtual environment in Ubuntu:

bash

deactivate
rm -rf ~/pipeseq_env
5. Remove all installed packages/tools in Ubuntu:

bash

sudo apt purge -y python3-venv build-essential zlib1g-dev \
  libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev \
  libcurl4-openssl-dev libssl-dev libsqlite3-dev wget curl \
  git unzip samtools hisat2 stringtie libgl1 libxkbcommon-x11-0
sudo apt autoremove -y
sudo apt clean
6. Remove PipeSeq files and configs:

bash

rm -rf ~/pipeSeq/settings.json
rm -rf ~/pipeSeq/results
rm -rf ~/pipeSeq/logs
rm -rf ~/pipeSeq
7. Uninstall SRA Toolkit (if installed):

powershell

winget uninstall SRA.SRA-Toolkit
8. Fully clean up WSL and Ubuntu (if no longer needed):

powershell

wsl --unregister Ubuntu
Support
For issues or questions, contact:
alexnerezenko@gmail.com
