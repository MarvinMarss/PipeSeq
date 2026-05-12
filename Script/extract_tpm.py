import os, json, re
import pandas as pd
import numpy as np
from collections import defaultdict

SETTINGS_FILE = "settings.json"
LOG_FILE = "extract_tpm_log.txt"

def log(msg):
    print(msg)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(msg + "\n")

def load_settings():
    if not os.path.exists(SETTINGS_FILE):
        raise FileNotFoundError("settings.json not found")
    with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
        st = json.load(f)
    results_folder = st["folders"].get("results_folder")
    gene_mapping = st.get("gene_mapping", {})
    if not results_folder:
        raise ValueError("results_folder is not set in settings.json")
    return results_folder, gene_mapping

def get_unique_filename(base_name, extension, folder):
    counter = 1
    file_name = f"{base_name}{extension}"
    while os.path.exists(os.path.join(folder, file_name)):
        file_name = f"{base_name}_{counter}{extension}"
        counter += 1
    return os.path.join(folder, file_name)

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

def normalize_base_name(file_name):
    # идентично extract_fpkm.py
    m = re.match(r"^(.*?)(Control)?\d+_sorted\.gtf$", file_name)
    if m:
        prefix, control = m.groups()
        return f"{prefix}{'_Control' if control else ''}"
    return file_name

def main():
    results_folder, gene_mapping = load_settings()
    gata_order = list(gene_mapping.values())

    src = find_latest_matching_file(results_folder, r"GTF_results_fpkm_all(?:_\d+)?\.txt")
    log(f"Используем входной файл TPM/FPKM: {src}")

    df_all = pd.read_csv(src, sep="\t")
    req = {"File","Gene ID","GATA Name","TPM"}
    missing = req - set(df_all.columns)
    if missing:
        raise ValueError(f"Во входном файле отсутствуют столбцы: {missing}")

    # 1) Полный файл TPM (зеркало *_fpkm_all, но только TPM)
    df_all_tpm = df_all[["File","Gene ID","GATA Name","TPM"]].copy()
    df_all_tpm["GATA Order"] = df_all_tpm["GATA Name"].apply(
        lambda x: gata_order.index(x.split("_")[0]) if x.split("_")[0] in gata_order else -1
    )
    df_all_tpm = df_all_tpm.sort_values(by=["GATA Order","File"]).drop(columns=["GATA Order"])
    out_all = get_unique_filename("GTF_results_tpm_all", ".txt", results_folder)
    df_all_tpm.to_csv(out_all, sep="\t", index=False)

    # 2) Средние по Base Name/IsControl
    df_all_tpm["Base Name"] = df_all_tpm["File"].apply(normalize_base_name)
    df_all_tpm["IsControl"] = df_all_tpm["Base Name"].str.endswith("_Control")
    avg_df = (df_all_tpm
              .groupby(["Base Name","IsControl","Gene ID","GATA Name"])["TPM"]
              .mean().reset_index())
    avg_df["GATA Order"] = avg_df["GATA Name"].apply(
        lambda x: gata_order.index(x.split("_")[0]) if x.split("_")[0] in gata_order else -1
    )
    avg_df = avg_df.sort_values(by=["GATA Order","Base Name","IsControl"]).drop(columns=["GATA Order"])
    out_avg = get_unique_filename("GTF_results_tpm_avg", ".txt", results_folder)
    avg_df.to_csv(out_avg, sep="\t", index=False)

    # 3) log2(Exp/Control) по TPM
    log2_rows = []
    exp_rows = avg_df[~avg_df["IsControl"]]
    for _, r in exp_rows.iterrows():
        base_clean = r["Base Name"].replace("_Control","")
        gene_id = r["Gene ID"]; name = r["GATA Name"]; tpm_exp = r["TPM"]
        ctrl = avg_df[(avg_df["Base Name"]==f"{base_clean}_Control") &
                      (avg_df["IsControl"]) &
                      (avg_df["Gene ID"]==gene_id) &
                      (avg_df["GATA Name"]==name)]
        if not ctrl.empty:
            tpm_ctrl = ctrl.iloc[0]["TPM"]
            log2val = np.log2(tpm_exp/tpm_ctrl) if (tpm_exp>0 and tpm_ctrl>0) else 0.0
        else:
            log2val = 0.0
            log(f"Контроль не найден для: {base_clean}, {gene_id}, {name}")
        log2_rows.append({"Base Name": base_clean,
                          "Gene ID": gene_id,
                          "GATA Name": name,
                          "log2(Exp/Control)": log2val})
    log2_df = pd.DataFrame(log2_rows)
    log2_df["GATA Order"] = log2_df["GATA Name"].apply(
        lambda x: gata_order.index(x.split("_")[0]) if x.split("_")[0] in gata_order else -1
    )
    log2_df = log2_df.sort_values(by=["GATA Order","Base Name"]).drop(columns=["GATA Order"])
    out_log2 = get_unique_filename("GTF_results_log2_tpm", ".txt", results_folder)
    log2_df.to_csv(out_log2, sep="\t", index=False)

    log(f"GTF_results_tpm_all сохранён: {out_all}")
    log(f"GTF_results_tpm_avg сохранён: {out_avg}")
    log(f"GTF_results_log2_tpm сохранён: {out_log2}")

if __name__ == "__main__":
    main()
