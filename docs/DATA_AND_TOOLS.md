# Data and Tool Policy

PipeSeq is source-controlled as code and documentation. Large local data and third-party binaries
stay outside git.

## Keep out of git

- FASTQ/FQ/SRA inputs
- SAM/BAM alignment outputs
- StringTie and featureCounts outputs
- Genome FASTA files and annotation GTF/GFF files
- HISAT2 indexes
- Local logs
- Local `Script/settings.json` and `Script/Mind.json`
- Downloaded third-party tool bundles

## Keep in git

- Python and batch source files
- Example configuration files
- README and usage documentation
- Small static assets such as `image.ico`
- GitHub templates and repository metadata

## Why this matters

GitHub has file-size limits, and sequencing/genome data can be several gigabytes. Keeping those
files local makes the repository cloneable and easier to review while preserving the same working
folder layout for users.
