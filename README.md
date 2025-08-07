Environment Setup on a New PC
Before first use, all dependencies must be installed.

1. Open PowerShell as Administrator
Press Win + X → select Windows PowerShell (Admin) or Terminal (Admin).

2. Install Python and WSL (Ubuntu)
In PowerShell, run:


winget install -e --id Python.Python.3.12
Restart the terminal, then install Python packages:



pip install PyQt6
pip install pandas numpy scipy seaborn matplotlib pyDESeq2
pip install pywin32
pip install statsmodels
Install WSL with Ubuntu:

wsl --install -d Ubuntu
Reboot the PC if prompted.

3. After reboot — open PowerShell (Admin):



wsl --install -d Ubuntu
wsl
Ubuntu will launch. Follow the prompts to set a username and password (any memorable values).

4. Inside Ubuntu — set up the environment:



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


python run_gui.py
Or step-by-step:



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

Removing the Environment from PC
To completely remove all installed programs and settings, follow these steps.

1. Open PowerShell as Administrator
Press Win + X → select Windows PowerShell (Admin) or Terminal (Admin).

2. Uninstall WSL and Ubuntu
Paste the following in PowerShell:

# Uninstall installed Python
winget uninstall Python.Python.3.12

# Uninstall all Python packages
pip uninstall -y PyQt6 pandas numpy scipy seaborn matplotlib pyDESeq2

pip uninstall statsmodels

# Uninstall WSL
wsl --uninstall
Wait for the uninstallation to finish.

3. Remove all Ubuntu settings and environments
After that, to delete all Ubuntu settings and environments, open PowerShell (Admin) again:

Make sure Ubuntu is completely removed.

4. Remove Python environment in Ubuntu
If you used a Python virtual environment, run the following commands:




# Deactivate the environment
deactivate

# Remove the isolated Python environment
rm -rf ~/pipeseq_env
5. Remove all installed packages and tools in Ubuntu



# Uninstall all installed packages
sudo apt purge -y python3-venv build-essential zlib1g-dev \
  libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev \
  libcurl4-openssl-dev libssl-dev libsqlite3-dev wget curl \
  git unzip samtools hisat2 stringtie libgl1 libxkbcommon-x11-0

# Clean up the system
sudo apt autoremove -y
sudo apt clean
6. Remove PipeSeq files and settings
To completely remove everything related to the pipeline, delete the following folders and files:

Delete the folder with genomes and annotations if it was downloaded.

Delete settings.json if you want to reset all settings:




rm -rf ~/pipeSeq/settings.json
Delete all intermediate files (e.g., result and log folders):




rm -rf ~/pipeSeq/results
rm -rf ~/pipeSeq/logs
7. Delete the PipeSeq program folder (if locally installed):



rm -rf ~/pipeSeq
8. Uninstall SRA Tools
If you installed the SRA Toolkit, you can remove it with:




winget uninstall SRA.SRA-Toolkit
9. Clean up the environment
If you're no longer going to use WSL and Ubuntu:

# Unregister Ubuntu from Windows
wsl --unregister Ubuntu
Support
alexnerezenko@gmail.com