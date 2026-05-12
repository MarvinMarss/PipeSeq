import os
import json
import pandas as pd
import numpy as np
from collections import defaultdict
import re

SETTINGS_FILE = "settings.json"
LOG_FILE = "extract_fpkm_log.txt"

def log(message):
    print(message)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(message + "\n")

def load_settings():
    if not os.path.exists(SETTINGS_FILE):
        log(f"Файл настроек {SETTINGS_FILE} не найден!")
        raise FileNotFoundError("settings.json not found")

    with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
        settings = json.load(f)

    gtf_folder = settings["folders"].get("gtf_folder")
    output_folder = settings["folders"].get("results_folder")
    gene_mapping = settings.get("gene_mapping", {})

    if not gtf_folder or not output_folder:
        log("В settings.json не указаны gtf_folder или results_folder!")
        raise ValueError("Missing folders in settings")

    return gtf_folder, output_folder, gene_mapping

def get_unique_filename(base_name, extension, folder):
    counter = 1
    file_name = f"{base_name}{extension}"
    while os.path.exists(os.path.join(folder, file_name)):
        file_name = f"{base_name}_{counter}{extension}"
        counter += 1
    return os.path.join(folder, file_name)

def normalize_base_name(file_name):

    match = re.match(r"^(.*?)(Control)?\d+_sorted\.gtf$", file_name)
    if match:
        prefix, control = match.groups()
        return f"{prefix}{'_Control' if control else ''}"
    return file_name

def extract_fpkm(gtf_folder, output_folder, gene_mapping):
    gata_order = list(gene_mapping.values())
    results = []
    gata4_count = defaultdict(int)

    for file_name in os.listdir(gtf_folder):
        if file_name.endswith(".gtf"):
            file_path = os.path.join(gtf_folder, file_name)
            with open(file_path, "r") as gtf_file:
                for line in gtf_file:
                    for gene_id, gene_name in gene_mapping.items():
                        if f'gene_id "{gene_id}"' in line and "FPKM" in line:
                            try:
                                fpkm_value = float([field for field in line.strip().split(";") if "FPKM" in field][0].split('"')[1])
                                tpm_value = float([field for field in line.strip().split(";") if "TPM" in field][0].split('"')[1])
                            except (ValueError, IndexError):
                                continue
                            name = gene_name
                            if name == "GATA-4":
                                gata4_count[file_name] += 1
                                name = f"GATA-4_t{gata4_count[file_name]}"
                            results.append({
                                "File": file_name,
                                "Gene ID": gene_id,
                                "GATA Name": name,
                                "FPKM": fpkm_value,
                                "TPM": tpm_value
                            })

    df = pd.DataFrame(results)
    df["GATA Order"] = df["GATA Name"].apply(lambda x: gata_order.index(x.split("_")[0]) if x.split("_")[0] in gata_order else -1)
    df = df.sort_values(by=["GATA Order", "File"]).drop(columns=["GATA Order"])

    output_fpkm_all = get_unique_filename("GTF_results_fpkm_all", ".txt", output_folder)
    df.to_csv(output_fpkm_all, sep="\t", index=False)

    df["Base Name"] = df["File"].apply(normalize_base_name)
    df["IsControl"] = df["Base Name"].str.endswith("_Control")

    avg_df = df.groupby(["Base Name", "IsControl", "Gene ID", "GATA Name"])["FPKM"].mean().reset_index()
    avg_df["GATA Order"] = avg_df["GATA Name"].apply(lambda x: gata_order.index(x.split("_")[0]) if x.split("_")[0] in gata_order else -1)
    avg_df = avg_df.sort_values(by=["GATA Order", "Base Name", "IsControl"]).drop(columns=["GATA Order"])

    output_fpkm_avg = get_unique_filename("GTF_results_fpkm_avg", ".txt", output_folder)
    avg_df.to_csv(output_fpkm_avg, sep="\t", index=False)

    log2_results = []


    exp_rows_all = avg_df[~avg_df["IsControl"]]

    for _, row in exp_rows_all.iterrows():
        base_clean = row["Base Name"].replace("_Control", "")
        gene_id = row["Gene ID"]
        name = row["GATA Name"]
        fpkm_exp = row["FPKM"]

        ctrl_rows = avg_df[
            (avg_df["Base Name"] == f"{base_clean}_Control") &
            (avg_df["IsControl"]) &
            (avg_df["Gene ID"] == gene_id) &
            (avg_df["GATA Name"] == name)
        ]

        if not ctrl_rows.empty:
            fpkm_ctrl = ctrl_rows.iloc[0]["FPKM"]
            if fpkm_exp > 0 and fpkm_ctrl > 0:
                log2_val = np.log2(fpkm_exp / fpkm_ctrl)
            else:
                log2_val = 0.0
        else:
            log2_val = 0.0
            log(f"Не найден контроль для: {base_clean}, {gene_id}, {name}")

        log2_results.append({
            "Base Name": base_clean,
            "Gene ID": gene_id,
            "GATA Name": name,
            "log2(Exp/Control)": log2_val
        })

    log2_df = pd.DataFrame(log2_results)
    log2_df["GATA Order"] = log2_df["GATA Name"].apply(lambda x: gata_order.index(x.split("_")[0]) if x.split("_")[0] in gata_order else -1)
    log2_df = log2_df.sort_values(by=["GATA Order", "Base Name"]).drop(columns=["GATA Order"])

    output_log2 = get_unique_filename("GTF_results_log2", ".txt", output_folder)
    log2_df.to_csv(output_log2, sep="\t", index=False)

    log(f"GTF_results_fpkm_all сохранён: {output_fpkm_all}")
    log(f"GTF_results_fpkm_avg сохранён: {output_fpkm_avg}")
    log(f"GTF_results_log2 сохранён: {output_log2}")

if __name__ == "__main__":
    gtf_folder, output_folder, gene_mapping = load_settings()
    extract_fpkm(gtf_folder, output_folder, gene_mapping)