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


def load_gtf_dir() -> Path:
    with SETTINGS_FILE.open("r", encoding="utf-8") as fh:
        settings = json.load(fh)
    return Path(settings["folders"]["gtf_folder"])


def split_gtf_filename(name: str):
    if name.endswith("_coverage.tsv"):
        sample = name[:-len("_coverage.tsv")]
        ext = "_coverage.tsv"
        return sample, ext
    if name.endswith(".gtf"):
        sample = name[:-len(".gtf")]
        ext = ".gtf"
        return sample, ext
    return None


def rewrite_sample_name(sample: str) -> str:
    for prefix in IMPLICIT_SINGLE_REPLICATE_PREFIXES:
        treated_old = f"{prefix}_sorted"
        control_old = f"{prefix}Control_sorted"
        if sample == treated_old:
            return f"{prefix}1_sorted"
        if sample == control_old:
            return f"{prefix}Control1_sorted"

    m = re.match(r"^(CC125Diflegilation60min(?:Control)?\d+)_1(_sorted)$", sample)
    if m:
        return f"{m.group(1)}{m.group(2)}"

    for prefix in MOVE_REPLICATE_AFTER_CONTROL_PREFIXES:
        m = re.match(rf"^({re.escape(prefix)})(\d+)Control(_sorted)$", sample)
        if m:
            return f"{m.group(1)}Control{m.group(2)}{m.group(3)}"

    return sample


def build_plan(gtf_dir: Path):
    plan = []
    for path in sorted(gtf_dir.iterdir(), key=lambda p: p.name.lower()):
        if not path.is_file():
            continue
        parsed = split_gtf_filename(path.name)
        if parsed is None:
            continue
        sample, ext = parsed
        new_sample = rewrite_sample_name(sample)
        if new_sample == sample:
            continue
        target = path.with_name(f"{new_sample}{ext}")
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


def write_log(gtf_dir: Path, plan):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = gtf_dir / f"rename_gtf_log_{stamp}.tsv"
    with log_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["old_name", "new_name"])
        for src, dst in plan:
            writer.writerow([src.name, dst.name])
    return log_path


def main():
    parser = argparse.ArgumentParser(
        description="Rename GTF files so sample names match PipeSeq README expectations."
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually rename files. Without this flag, only print the plan.",
    )
    args = parser.parse_args()

    gtf_dir = load_gtf_dir()
    if not gtf_dir.exists():
        print(f"GTF directory not found: {gtf_dir}", file=sys.stderr)
        sys.exit(1)

    plan = build_plan(gtf_dir)
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

    log_path = write_log(gtf_dir, plan)
    for src, dst in plan:
        src.rename(dst)
    print(f"Applied {len(plan)} renames.")
    print(f"Log written to: {log_path}")


if __name__ == "__main__":
    main()
