# Contributing

Thank you for improving PipeSeq.

## Development setup

1. Install Python 3.10 or newer.
2. Create and activate a virtual environment.
3. Install dependencies:

```bash
pip install -r requirements-dev.txt
```

The pipeline also requires command-line bioinformatics tools in WSL or on `PATH`:
HISAT2, SAMtools, StringTie and Subread/featureCounts.

## Repository rules

- Keep source code and documentation in git.
- Do not commit FASTQ/SRA files, BAM/SAM files, genome FASTA/GTF files, HISAT2 indexes,
  generated result tables, logs, local `settings.json`, or third-party tool bundles.
- Use `Script/settings.example.json` as the template for local configuration.
- Keep pull requests focused on one behavior change or documentation update.

## Checks before a pull request

```bash
python -m py_compile Script/*.py
ruff check Script
```

If a full pipeline run is possible, mention the dataset, genome version and command/tool
versions used for verification.
