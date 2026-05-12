import argparse
import csv
import json
import re
import sys
from datetime import datetime
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SETTINGS_FILE = SCRIPT_DIR / "settings.json"


MOVE_REPLICATE_AFTER_CONTROL_PREFIXES = (
    "Acidicrep",
    "CO2rep",
    "Chlamylinnemannicu",
    "Chlamylinnemanniwt",
)

IMPLICIT_SINGLE_REPLICATE_PREFIXES = (
    "05DarkAnoxic",
    "6DarkAnoxic",
)


def load_counts_dir() -> Path:
    with SETTINGS_FILE.open("r", encoding="utf-8") as fh:
        settings = json.load(fh)
    results_dir = Path(settings["folders"]["results_folder"])
    return results_dir / "Counts"


def split_count_filename(name: str):
    prefix = "gene_counts_"
    if not name.startswith(prefix):
        return None
    if name.endswith(".txt.summary"):
        sample = name[len(prefix):-len(".txt.summary")]
        ext = ".txt.summary"
        return sample, ext
    if name.endswith(".txt"):
        sample = name[len(prefix):-len(".txt")]
        ext = ".txt"
        return sample, ext
    return None


def rewrite_sample_name(sample: str) -> str:
    for prefix in IMPLICIT_SINGLE_REPLICATE_PREFIXES:
        treated_old = f"{prefix}_single_sorted"
        control_old = f"{prefix}Control_single_sorted"
        if sample == treated_old:
            return f"{prefix}1_single_sorted"
        if sample == control_old:
            return f"{prefix}Control1_single_sorted"

    m = re.match(r"^(CC125Diflegilation60min(?:Control)?\d+)_1(_single_sorted)$", sample)
    if m:
        return f"{m.group(1)}{m.group(2)}"

    for prefix in MOVE_REPLICATE_AFTER_CONTROL_PREFIXES:
        m = re.match(
            rf"^({re.escape(prefix)})(\d+)Control(_(?:paired|single)_sorted)$",
            sample,
        )
        if m:
            return f"{m.group(1)}Control{m.group(2)}{m.group(3)}"

    return sample


def build_plan(counts_dir: Path):
    plan = []
    for path in sorted(counts_dir.iterdir(), key=lambda p: p.name.lower()):
        if not path.is_file():
            continue
        parsed = split_count_filename(path.name)
        if parsed is None:
            continue
        sample, ext = parsed
        new_sample = rewrite_sample_name(sample)
        if new_sample == sample:
            continue
        new_name = f"gene_counts_{new_sample}{ext}"
        target = path.with_name(new_name)
        plan.append((path, target))
    return plan


def validate_plan(plan):
    problems = []
    seen_targets = {}
    for src, dst in plan:
        if dst.exists() and dst != src:
            problems.append(f"Target already exists: {dst}")
        key = str(dst).lower()
        previous = seen_targets.get(key)
        if previous is not None and previous != src:
            problems.append(f"Two sources map to one target: {previous} and {src} -> {dst}")
        seen_targets[key] = src
    return problems


def write_log(counts_dir: Path, plan):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = counts_dir / f"rename_counts_log_{stamp}.tsv"
    with log_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["old_name", "new_name"])
        for src, dst in plan:
            writer.writerow([src.name, dst.name])
    return log_path


def main():
    parser = argparse.ArgumentParser(
        description="Rename count files so sample names match PipeSeq README expectations."
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually rename files. Without this flag, only print the plan.",
    )
    args = parser.parse_args()

    counts_dir = load_counts_dir()
    if not counts_dir.exists():
        print(f"Counts directory not found: {counts_dir}", file=sys.stderr)
        sys.exit(1)

    plan = build_plan(counts_dir)
    if not plan:
        print("No files need renaming.")
        return

    problems = validate_plan(plan)
    if problems:
        print("Cannot continue because of conflicts:", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
        sys.exit(2)

    print("Planned renames:")
    for src, dst in plan:
        print(f"  {src.name} -> {dst.name}")

    if not args.apply:
        print(f"Dry-run complete. Total files to rename: {len(plan)}")
        return

    log_path = write_log(counts_dir, plan)
    for src, dst in plan:
        src.rename(dst)
    print(f"Applied {len(plan)} renames.")
    print(f"Log written to: {log_path}")


if __name__ == "__main__":
    main()
