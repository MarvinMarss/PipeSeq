import os
import json
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import sys

SETTINGS_FILE = "settings.json"
DEBUG_LOG_FILE = "GTF_results_pvalues_log.txt"

def log(message, logfile=DEBUG_LOG_FILE):
    print(message)
    with open(logfile, "a", encoding="utf-8") as f:
        f.write(message + "\n")

def load_settings():
    if not os.path.exists(SETTINGS_FILE):
        log(f"Файл настроек {SETTINGS_FILE} не найден!")
        sys.exit(1)
    with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
        settings = json.load(f)
    folders = settings.get("folders", {})
    results_folder = folders.get("results_folder")
    if not results_folder:
        log("Не задан путь к папке GTF (results_folder) в settings.json!")
        sys.exit(1)
    return results_folder

def get_unique_filename(base_name, extension, folder):
    counter = 1
    file_name = f"{base_name}{extension}"
    while os.path.exists(os.path.join(folder, file_name)):
        file_name = f"{base_name}_{counter}{extension}"
        counter += 1
    return os.path.join(folder, file_name)

def get_base_name(file_name):
    name = file_name.replace("_merged.gtf", "").replace("_midel_merged.gtf", "").replace("_sorted.gtf", "")
    return name.rstrip("0123456789")

def main():
    if os.path.exists(DEBUG_LOG_FILE):
        os.remove(DEBUG_LOG_FILE)
    log("Начинаем обработку данных...")

    results_folder = load_settings()
    input_file = os.path.join(results_folder, "GTF_results_fpkm_all.txt")
    if not os.path.exists(input_file):
        log(f"Входной файл {input_file} не найден!")
        sys.exit(1)
    log(f"Загружаем данные из {input_file}")
    df = pd.read_csv(input_file, sep="\t")
    log(f"Загружено {len(df)} строк.")

    df["Base Name"] = df["File"].apply(get_base_name)

    exp_df = df[~df["Base Name"].str.contains("Control")]
    ctrl_df = df[df["Base Name"].str.contains("Control")]
    log(f"Экспериментальных записей: {len(exp_df)}; Контрольных: {len(ctrl_df)}")

    grouped_exp = exp_df.groupby(["Base Name", "Gene ID", "GATA Name"])["FPKM"].apply(list).reset_index()
    grouped_ctrl = ctrl_df.groupby(["Base Name", "Gene ID", "GATA Name"])["FPKM"].apply(list).reset_index()
    grouped_ctrl["Base Name"] = grouped_ctrl["Base Name"].str.replace("Control", "", regex=False)

    merged = pd.merge(grouped_exp, grouped_ctrl, on=["Base Name", "Gene ID", "GATA Name"],
                      suffixes=("_exp", "_ctrl"), how="inner")
    log(f"После объединения осталось {len(merged)} групп для расчёта p-value.")

    if merged.empty:
        log("Недостаточно совпадений для расчёта p-value. Создаётся заглушка на основе всех уникальных комбинаций.")
        unique_combos = df[["Base Name", "Gene ID", "GATA Name"]].drop_duplicates()
        unique_combos["p-value"] = 1.0
        final_df = unique_combos
    else:
        merged["FPKM_exp"] = merged["FPKM_exp"].apply(lambda x: list(map(float, x)))
        merged["FPKM_ctrl"] = merged["FPKM_ctrl"].apply(lambda x: list(map(float, x)))

        def calc_pvalue(row):
            exp_values = row["FPKM_exp"]
            ctrl_values = row["FPKM_ctrl"]
            if len(exp_values) > 1 and len(ctrl_values) > 1:
                return ttest_ind(exp_values, ctrl_values, equal_var=False)[1]
            else:
                return 1.0

        merged["p-value"] = merged.apply(calc_pvalue, axis=1)
        merged["p-value"] = merged["p-value"].fillna(1.0)
        final_df = merged[["Base Name", "Gene ID", "GATA Name", "p-value"]].drop_duplicates()

    final_df = final_df.sort_values(by=["Base Name", "Gene ID", "GATA Name"])
    output_file = get_unique_filename("GTF_results_pvalues", ".txt", results_folder)
    final_df.to_csv(output_file, sep="\t", index=False)
    log(f"Результаты p-value сохранены в {output_file}")

if __name__ == "__main__":
    main()
