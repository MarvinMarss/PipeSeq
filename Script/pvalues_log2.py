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

def extract_order(name):
    base = name.split("_")[0]
    match = re.search(r"(\d+(?:\.\d+)?)", base)
    return float(match.group(1)) if match else float('inf')

def main():
    input_folder = load_settings()

    pvalue_file = find_latest_matching_file(input_folder, r"GTF_results_pvalues(?:_\d+)?\.txt")
    log2_file = find_latest_matching_file(input_folder, r"GTF_results_log2(?:_\d+)?\.txt")

    log(f"Используем p-value файл: {os.path.basename(pvalue_file)}")
    log(f"Используем log2 файл: {os.path.basename(log2_file)}")

    pvalues_df = pd.read_csv(pvalue_file, sep="\t")
    log2_df = pd.read_csv(log2_file, sep="\t")

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
