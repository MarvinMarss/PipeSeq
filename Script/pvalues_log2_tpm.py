import os, json
import pandas as pd
import re

SETTINGS_FILE = "settings.json"

def log(m): print(m)

def load_results_folder():
    with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
        st = json.load(f)
    rf = st.get("folders", {}).get("results_folder")
    if not rf: raise ValueError("results_folder not set")
    return rf

def get_unique_filename(base, ext, folder):
    k=1; name=f"{base}{ext}"
    while os.path.exists(os.path.join(folder, name)):
        name=f"{base}_{k}{ext}"; k+=1
    return os.path.join(folder, name)

def find_latest_matching_file(folder, pattern):
    matches = []
    for name in os.listdir(folder):
        if re.fullmatch(pattern, name):
            path = os.path.join(folder, name)
            matches.append((os.path.getmtime(path), name, path))
    if not matches:
        raise FileNotFoundError(f"Нет файлов, подходящих под шаблон: {pattern}")
    matches.sort()
    return matches[-1][2]

def extract_order(name: str) -> float:
    base = name.split("_")[0]
    m = re.search(r"(\d+(?:\.\d+)?)", base)
    return float(m.group(1)) if m else float('inf')

def main():
    rf = load_results_folder()

    pfile = find_latest_matching_file(rf, r"GTF_results_pvalues_tpm(?:_\d+)?\.txt")
    l2file = find_latest_matching_file(rf, r"GTF_results_log2_tpm(?:_\d+)?\.txt")
    log(f"Используем p-values (TPM): {os.path.basename(pfile)}")
    log(f"Используем log2 (TPM): {os.path.basename(l2file)}")

    p_df = pd.read_csv(pfile, sep="\t")
    l2_df = pd.read_csv(l2file, sep="\t")

    for c in ("Base Name","Gene ID","GATA Name"):
        p_df[c] = p_df[c].astype(str).str.strip()
        l2_df[c]= l2_df[c].astype(str).str.strip()

    merged = pd.merge(p_df, l2_df, on=["Base Name","Gene ID","GATA Name"], how="inner")
    if merged.empty:
        log("После объединения нет совпадений."); raise SystemExit(1)

    merged["GATA Order"] = merged["GATA Name"].apply(extract_order)
    merged = merged.sort_values(by=["Base Name","GATA Order"]).drop(columns=["GATA Order"]).drop_duplicates()

    outfile = get_unique_filename("Stringtie_TPM", ".txt", rf)
    merged.to_csv(outfile, sep="\t", index=False)
    log(f"Объединённый файл сохранён: {outfile}")

if __name__ == "__main__":
    main()
