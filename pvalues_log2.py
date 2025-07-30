import os
import json
import pandas as pd
import numpy as np
import re

SETTINGS_FILE = "settings.json"

def log(message):
    print(message)

def load_settings():
    if not os.path.exists(SETTINGS_FILE):
        log(f"Файл настроек {SETTINGS_FILE} не найден!")
        exit(1)
    with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
        settings = json.load(f)
    results_folder = settings.get("folders", {}).get("results_folder")
    if not results_folder:
        log("В settings.json не задан путь к results_folder.")
        exit(1)
    return results_folder

def get_unique_filename(base_name, extension, folder):
    counter = 1
    file_name = f"{base_name}{extension}"
    while os.path.exists(os.path.join(folder, file_name)):
        file_name = f"{base_name}_{counter}{extension}"
        counter += 1
    return os.path.join(folder, file_name)

def extract_order(name):
    base = name.split("_")[0]
    match = re.search(r"(\d+(?:\.\d+)?)", base)
    return float(match.group(1)) if match else float('inf')

def main():
    input_folder = load_settings()

    pvalue_files = sorted([f for f in os.listdir(input_folder) if f.startswith("GTF_results_pvalues") and f.endswith(".txt")])
    log2_files = sorted([f for f in os.listdir(input_folder) if f.startswith("GTF_results_log2") and f.endswith(".txt")])

    log(f"Найдено файлов p-values: {len(pvalue_files)} -> {pvalue_files}")
    log(f"Найдено файлов log2: {len(log2_files)} -> {log2_files}")

    if not pvalue_files or not log2_files:
        log("Нет нужных файлов для объединения.")
        exit(1)

    pvalue_dfs = []
    for f in pvalue_files:
        try:
            df = pd.read_csv(os.path.join(input_folder, f), sep="\t")
            log(f"Загружен p-value файл: {f}, строк: {len(df)}")
            pvalue_dfs.append(df)
        except Exception as e:
            log(f"Ошибка при загрузке {f}: {e}")

    log2_dfs = []
    for f in log2_files:
        try:
            df = pd.read_csv(os.path.join(input_folder, f), sep="\t")
            log(f"Загружен log2 файл: {f}, строк: {len(df)}")
            log2_dfs.append(df)
        except Exception as e:
            log(f"Ошибка при загрузке {f}: {e}")

    if not pvalue_dfs or not log2_dfs:
        log("Не удалось загрузить данные.")
        exit(1)

    pvalues_df = pd.concat(pvalue_dfs, ignore_index=True)
    log2_df = pd.concat(log2_dfs, ignore_index=True)

    pvalues_df["Base Name"] = pvalues_df["Base Name"].str.strip()
    log2_df["Base Name"] = log2_df["Base Name"].str.strip()

    merged_df = pd.merge(pvalues_df, log2_df, on=["Base Name", "Gene ID", "GATA Name"], how="inner")

    if merged_df.empty:
        log("После объединения нет совпадений по Base Name, Gene ID и GATA Name.")
        exit(1)

    merged_df["GATA Order"] = merged_df["GATA Name"].apply(extract_order)
    merged_df = merged_df.sort_values(by=["Base Name", "GATA Order"]).drop(columns=["GATA Order"])
    merged_df = merged_df.drop_duplicates()

    output_file = get_unique_filename("Stringtie", ".txt", input_folder)
    merged_df.to_csv(output_file, sep="\t", index=False)

    log(f"Объединённый файл сохранён: {output_file}")

if __name__ == "__main__":
    main()
