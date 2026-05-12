# HISAT2 index folder

PipeSeq expects the HISAT2 index files here when a prebuilt index is used.

The expected basename is usually `genome_index`, producing files such as:

```text
genome_index.1.ht2
genome_index.2.ht2
...
genome_index.8.ht2
```

Index files are generated artifacts and are intentionally ignored by git.
