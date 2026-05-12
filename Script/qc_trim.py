
import os, sys, json, subprocess, glob, shutil, time, re, gzip, io, zipfile

SETTINGS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "settings.json")
LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "qc_trim_log.txt")

# ------------------------ Logging ------------------------

def log(msg: str):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(line + "\n")

def tool_versions():
    def ver(cmd):
        try:
            out = subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True)
            first = out.strip().splitlines()[0] if out.strip() else "ok"
            return f"{cmd[0]}: {first}"
        except Exception as e:
            return f"{cmd[0]}: {e}"
    log("Tool versions:")
    for cmd in (["cutadapt", "--version"], ["fastqc", "--version"], ["multiqc", "--version"]):
        log("  " + ver(cmd))

# ------------------------ Settings ------------------------

def load_settings():
    if not os.path.exists(SETTINGS_FILE):
        log(f"settings.json not found: {SETTINGS_FILE}")
        sys.exit(1)
    with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
        s = json.load(f)
    fastq_folder = s.get("folders", {}).get("fastq_folder")
    results_folder = s.get("folders", {}).get("results_folder")
    if not fastq_folder or not results_folder:
        log("folders.fastq_folder and folders.results_folder are required in settings.json")
        sys.exit(1)
    if not os.path.isdir(fastq_folder):
        log(f"FASTQ folder not found: {fastq_folder}")
        sys.exit(1)
    os.makedirs(results_folder, exist_ok=True)
    qc_opts = s.get("qc", {})
    # defaults for new options
    qc_opts.setdefault("auto", True)
    qc_opts.setdefault("enable_rrna_filter", False)
    qc_opts.setdefault("rrna_tool", "sortmerna")
    qc_opts.setdefault("rrna_db", "")
    return fastq_folder, results_folder, qc_opts

# ------------------------ FASTQ discovery ------------------------

def _is_pair1(name: str) -> bool:
    return name.endswith("_1.fastq") or name.endswith("_1.fastq.gz")

def _pair_base(name: str) -> str:
    if name.endswith("_1.fastq.gz"):
        return name[:-len("_1.fastq.gz")]
    if name.endswith("_1.fastq"):
        return name[:-len("_1.fastq")]
    return name

def _mate_name(base: str, gz: bool) -> str:
    return f"{base}_2.fastq.gz" if gz else f"{base}_2.fastq"

def find_fastqs(fastq_dir, sample_prefix=None):
    files = [f for f in os.listdir(fastq_dir) if (f.endswith(".fastq") or f.endswith(".fastq.gz"))]
    if sample_prefix:
        files = [f for f in files if f.startswith(sample_prefix)]
    paired, single = {}, []
    file_set = set(files)
    for f in files:
        # try to form pairs based on _1/_2 convention
        if _is_pair1(f):
            base = _pair_base(f)
            gz = f.endswith(".gz")
            mate = _mate_name(base, gz)
            if mate in file_set:
                paired[base] = (os.path.join(fastq_dir, f), os.path.join(fastq_dir, mate))
            else:
                single.append(os.path.join(fastq_dir, f))
        elif not (f.endswith("_2.fastq") or f.endswith("_2.fastq.gz")):
            single.append(os.path.join(fastq_dir, f))
    # deduplicate singles that might be part of pairs
    singles_clean = []
    paired_paths = {p for pair in paired.values() for p in pair}
    for p in set(os.path.join(fastq_dir, f) for f in files if not f.endswith("_2.fastq") and not f.endswith("_2.fastq.gz")):
        if p not in paired_paths and p not in singles_clean:
            singles_clean.append(p)
    return paired, singles_clean

# ------------------------ FastQC + MultiQC ------------------------

def run_fastqc(inputs, out_dir):
    if not inputs:
        return
    os.makedirs(out_dir, exist_ok=True)
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
    FASTQC_DIR = os.path.join(PROJECT_DIR, "Tool", "FastQC")

    candidates = [
        os.path.join(FASTQC_DIR, "fastqc.bat"),
        os.path.join(FASTQC_DIR, "run_fastqc.bat"),
        shutil.which("fastqc"),
    ]
    fastqc_path = next((p for p in candidates if p and os.path.exists(p)), None)
    if not fastqc_path:
        log("FastQC not found (Tool\\FastQC or PATH).")
        sys.exit(1)

    t0 = time.time(); window = 3600
    def _is_fresh(p):
        try:
            return (t0 - os.path.getmtime(p)) <= window
        except Exception:
            return False

    def collect(dst):
        moved = 0
        search_dirs = {FASTQC_DIR, os.path.dirname(os.path.abspath(inputs[0])), dst}
        patterns = ["*_fastqc.html", "*_fastqc.zip", "*_fastqc/fastqc_data.txt"]
        for src in list(search_dirs):
            for pat in patterns:
                for f in glob.glob(os.path.join(src, pat)):
                    if not _is_fresh(f):
                        continue
                    try:
                        # If it's already inside a folder like sample_fastqc/fastqc_data.txt, copy the whole folder
                        if f.endswith("fastqc_data.txt"):
                            folder = os.path.dirname(f)  # .../_fastqc
                            target = os.path.join(dst, os.path.basename(folder))
                            if os.path.abspath(folder) != os.path.abspath(target):
                                if os.path.exists(target):
                                    shutil.rmtree(target, ignore_errors=True)
                                shutil.copytree(folder, target)
                                moved += 1
                        else:
                            tgt = os.path.join(dst, os.path.basename(f))
                            if os.path.abspath(os.path.dirname(f)) != os.path.abspath(dst):
                                shutil.move(f, tgt)
                            moved += 1
                    except Exception as e:
                        log(f"collect warning: {e}")
        if moved == 0:
            for pat in ["*_fastqc.html", "*_fastqc.zip"]:
                for f in glob.glob(os.path.join(PROJECT_DIR, "**", pat), recursive=True):
                    if not _is_fresh(f):
                        continue
                    try:
                        tgt = os.path.join(dst, os.path.basename(f))
                        if os.path.abspath(os.path.dirname(f)) != os.path.abspath(dst):
                            shutil.move(f, tgt)
                        moved += 1
                    except Exception as e:
                        log(f"collect recursive warning: {e}")
        return moved

    cmd = (["cmd", "/c", fastqc_path] if fastqc_path.lower().endswith((".bat", ".cmd")) else [fastqc_path]) + inputs
    cp = subprocess.run(cmd, cwd=FASTQC_DIR)
    if cp.returncode not in (0, 1):
        log(f"FastQC exited with code {cp.returncode}")
    moved = collect(out_dir)
    have = any(glob.glob(os.path.join(out_dir, "*_fastqc.html")) + glob.glob(os.path.join(out_dir, "*_fastqc.zip")))
    if not (moved > 0 or have):
        log("FastQC produced no reports.");
        sys.exit(1)

def run_multiqc(qc_root):
    os.makedirs(qc_root, exist_ok=True)
    exe = shutil.which("multiqc")
    cmd = [exe, "-o", qc_root, qc_root] if exe else [sys.executable, "-m", "multiqc", "-o", qc_root, qc_root]
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        log(f"MultiQC unavailable ({e}); skipped.")

# ------------------------ FastQC parsing & auto-params ------------------------

def _read_fastqc_data_from_zip(zip_path):
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            for nm in zf.namelist():
                if nm.endswith("fastqc_data.txt"):
                    with zf.open(nm) as fh:
                        return fh.read().decode("utf-8", errors="ignore").splitlines()
    except Exception as e:
        log(f"zip parse warning: {e}")
    return []

def _fastqc_data_files(qc_dir):
    txts = []
    # unpacked folders
    for p in glob.glob(os.path.join(qc_dir, "*_fastqc", "fastqc_data.txt")):
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                txts.append(f.read().splitlines())
        except Exception as e:
            log(f"read fastqc_data warning: {e}")
    # zips
    for z in glob.glob(os.path.join(qc_dir, "*_fastqc.zip")):
        lines = _read_fastqc_data_from_zip(z)
        if lines:
            txts.append(lines)
    return txts

def parse_fastqc_data(qc_dir):
    datasets = _fastqc_data_files(qc_dir)
    adapters, two_color, q_profile = [], False, []
    for lines in datasets:
        # platform heuristic
        if any(("NextSeq" in ln) or ("NovaSeq" in ln) for ln in lines):
            two_color = True
        # Per base sequence quality
        try:
            # Find module start
            idx_header = None
            for i, l in enumerate(lines):
                if l.strip().startswith("#Base"):
                    idx_header = i
            if idx_header is not None:
                i = idx_header + 1
                while i < len(lines) and not lines[i].startswith(">>END_MODULE"):
                    parts = lines[i].split()
                    # expected columns: Base  Mean  Median ...
                    if len(parts) >= 6:
                        base_bin = parts[0]  # e.g. "1-10"
                        try:
                            medianQ = float(parts[5])  # 6th column usually Median
                        except Exception:
                            medianQ = None
                        if medianQ is not None:
                            q_profile.append((base_bin, medianQ))
                    i += 1
        except Exception:
            pass
        # Overrepresented sequences / Adapter Content
        in_over = False
        for l in lines:
            if l.startswith(">>Overrepresented sequences"):
                in_over = True;
                continue
            if in_over and l.startswith(">>END_MODULE"):
                in_over = False
            if in_over and "\t" in l:
                # columns: Sequence  Count  Percentage  Possible Source
                try:
                    seq = l.split("\t")[0].strip()
                    src = l.split("\t")[-1].lower()
                    if "adapter" in src and 18 <= len(seq) <= 36:
                        adapters.append(seq)
                except Exception:
                    pass
    # unique keep order
    adapters = list(dict.fromkeys(adapters))
    return adapters, {"two_color": two_color}, q_profile

def _top_end(b: str) -> int:
    # base bin like "1-10" or "101-110" -> return upper bound
    try:
        return int(str(b).split("-")[-1])
    except Exception:
        try:
            return int(b)
        except Exception:
            return 0

def recommend_q_minlen(q_profile, default_q=20, default_minlen=50):
    if not q_profile:
        return default_q, default_minlen
    # target median Q ~ 25 as compromise for RNA-seq
    target = 25
    try:
        L = max((_top_end(b) for (b, med) in q_profile if med >= target), default_minlen)
        # clamp to reasonable range
        L = max(35, min(L, 150))
        return 25, L
    except Exception:
        return default_q, default_minlen

def auto_params(pre_qc_dir):
    adapters, flags, qprof = parse_fastqc_data(pre_qc_dir)
    q, minlen = recommend_q_minlen(qprof, default_q=20, default_minlen=50)
    nextseq_trim = 20 if flags.get("two_color", False) else None
    if not adapters:
        adapters = [
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",  # R1 (TruSeq)
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"   # R2
        ]
    a1 = adapters[0]
    a2 = adapters[1] if len(adapters) > 1 else adapters[0]
    return {
        "q": q,
        "minlen_pe": minlen,
        "minlen_se": max(30, minlen - 15),
        "adapter_r1": a1,
        "adapter_r2": a2,
        "nextseq_trim": nextseq_trim
    }

# ------------------------ Cutadapt wrappers ------------------------

def cutadapt_pe(r1_in, r2_in, r1_out, r2_out, a1, a2, q=20, minlen=50, nextseq=None):
    if not a1: a1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    if not a2: a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    cmd = ["cutadapt", "-j", "4", "-a", a1, "-A", a2, "-q", f"{q},{q}", "--minimum-length", str(minlen)]
    if nextseq:
        cmd += ["--nextseq-trim", str(nextseq)]
    cmd += ["-o", r1_out, "-p", r2_out, r1_in, r2_in]
    log("Running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)

def cutadapt_se(r_in, r_out, a, q=20, minlen=35, nextseq=None):
    if not a: a = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    cmd = ["cutadapt", "-j", "4", "-a", a, "-q", str(q), "--minimum-length", str(minlen)]
    if nextseq:
        cmd += ["--nextseq-trim", str(nextseq)]
    cmd += ["-o", r_out, r_in]
    log("Running: " + " ".join(cmd))
    subprocess.run(cmd, check=True)

# ------------------------ Optional rRNA/organellar filtering (stub/hook) ------------------------

def rrna_filter_if_enabled(paths, qc_opts):
    if not qc_opts.get("enable_rrna_filter", False):
        return
    tool = qc_opts.get("rrna_tool", "sortmerna")
    db = qc_opts.get("rrna_db", "")
    if not db:
        log("rRNA filter requested but 'rrna_db' not specified; skipping.")
        return
    log(f"rRNA filtering enabled ({tool}) — not implemented here; provide your own hook/command.")
    # Placeholder: users can inject SortMeRNA/BBduk here if needed.

# ------------------------ Main ------------------------

def main():
    sample_prefix = None
    TRIM = False
    KEEP_RAW = True

    args = sys.argv[1:]
    if args and not args[0].startswith("--"):
        sample_prefix = args[0]; args = args[1:]
    i = 0
    while i < len(args):
        if args[i] == "--trim" and i+1 < len(args):
            TRIM = args[i+1].lower() == "on"; i += 2
        elif args[i] == "--keep-raw" and i+1 < len(args):
            KEEP_RAW = args[i+1].lower() == "yes"; i += 2
        else:
            i += 1

    # housekeeping
    if os.path.exists(LOG_FILE):
        try: os.remove(LOG_FILE)
        except Exception: pass
    tool_versions()

    fastq_dir, results_dir, qc_opts = load_settings()
    qc_root = os.path.join(results_dir, "qc")
    qc_pre  = os.path.join(qc_root, "raw")
    qc_post = os.path.join(qc_root, "trimmed")
    os.makedirs(qc_pre, exist_ok=True); os.makedirs(qc_post, exist_ok=True)

    paired, single = find_fastqs(fastq_dir, sample_prefix)
    pre_inputs = [p for pair in paired.values() for p in pair] + list(single)

    if not pre_inputs:
        log("No FASTQ files found.");
        sys.exit(1)

    log(f"Discovered {len(paired)} paired and {len(single)} single-end files.")
    run_fastqc(pre_inputs, qc_pre)

    # Optional auto-parameterization
    auto = None
    if qc_opts.get("auto", True):
        try:
            auto = auto_params(qc_pre)
            log(f"Auto-parameters: q={auto['q']}, minlen_pe={auto['minlen_pe']}, minlen_se={auto['minlen_se']}, nextseq_trim={auto['nextseq_trim']}, adapters=({auto['adapter_r1'][:10]}..., {auto['adapter_r2'][:10]}...)")
        except Exception as e:
            log(f"Auto-parameterization failed: {e} (falling back to defaults)")

    def pick(key, default):
        val = qc_opts.get(key, None)
        if val is None and auto is not None:
            return auto.get(key, default)
        return val if val is not None else default

    trim_q   = pick("trim_q", 20)
    minlenPE = pick("trim_minlen_pe", 50)
    minlenSE = pick("trim_minlen_se", 35)
    a1 = pick("adapter_r1", None)
    a2 = pick("adapter_r2", None)
    nextseq = pick("nextseq_trim", None)

    if TRIM:
        trimmed_dir = os.path.join(fastq_dir, "fastq_trimmed"); os.makedirs(trimmed_dir, exist_ok=True)
        backup_dir  = os.path.join(fastq_dir, "fastq_backup")
        if KEEP_RAW: os.makedirs(backup_dir, exist_ok=True)

        # Paired-end
        for base, (r1, r2) in paired.items():
            r1_tmp = os.path.join(trimmed_dir, os.path.basename(r1))
            r2_tmp = os.path.join(trimmed_dir, os.path.basename(r2))
            cutadapt_pe(r1, r2, r1_tmp, r2_tmp, a1=a1, a2=a2, q=trim_q, minlen=minlenPE, nextseq=nextseq)

            if KEEP_RAW:
                shutil.move(r1, os.path.join(backup_dir, os.path.basename(r1)))
                shutil.move(r2, os.path.join(backup_dir, os.path.basename(r2)))
            else:
                os.remove(r1); os.remove(r2)
            shutil.move(r1_tmp, r1); shutil.move(r2_tmp, r2)

        # Single-end
        for f in single:
            r_tmp = os.path.join(trimmed_dir, os.path.basename(f))
            cutadapt_se(f, r_tmp, a=a1, q=trim_q, minlen=minlenSE, nextseq=nextseq)
            if KEEP_RAW:
                shutil.move(f, os.path.join(backup_dir, os.path.basename(f)))
            else:
                os.remove(f)
            shutil.move(r_tmp, f)

        # Optional rRNA/organellar filtering (hook)
        rrna_filter_if_enabled(pre_inputs, qc_opts)

        post_inputs = [p for pair in paired.values() for p in pair] + list(single)
        run_fastqc(post_inputs, qc_post)

    run_multiqc(qc_root)
    log("QC/trim finished.")

if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        log(f"External command failed: {e}"); sys.exit(1)
    except Exception as e:
        log(f"Unhandled error: {e}"); sys.exit(1)
