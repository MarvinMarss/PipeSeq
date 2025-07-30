import os
import json
import pandas as pd
import re

SETTINGS_FILE = "settings.json"
LOG_FILE = "extract_deseq2_log.txt"

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
    
    output_folder = settings["folders"].get("results_folder")
    gene_mapping = settings.get("gene_mapping", {})
    
    if not output_folder:
        log("В settings.json не указана папка results_folder!")
        raise ValueError("Missing results_folder")
    
    return output_folder, gene_mapping

def extract_gata_number(gata_name):
    """
    Извлекает числовую часть из строки GATA Name.
    Например, для "GATA-4_t1" вернется 4.
    Если число не найдено, возвращает 999.
    """
    m = re.search(r'GATA[-\s]*(\d+)', str(gata_name))
    if m:
        return int(m.group(1))
    else:
        return 999

def extract_genes(output_folder, gene_mapping):

    all_files = os.listdir(output_folder)
    result_files = [f for f in all_files if f.startswith("results_Deseq2") and f.endswith(".tsv")]
    
    if not result_files:
        log("В папке нет файлов, начинающихся на 'results_Deseq2' с расширением .tsv")
        raise FileNotFoundError("No results_Deseq2*.tsv files found")
    
    dfs = []
    for file in result_files:
        file_path = os.path.join(output_folder, file)
        try:
            df = pd.read_csv(file_path, sep='\t')
        except Exception as e:
            log(f"Не удалось прочитать файл {file}: {e}")
            continue


        if "Gene ID" not in df.columns:
            if df.index.name in ["gene", "Gene ID"]:
                df = df.reset_index()
            else:
                log(f"Файл {file} не содержит колонки 'Gene ID'. Пропускаю его.")
                continue


        filtered_df = df[df["Gene ID"].isin(gene_mapping.keys())].copy()
        filtered_df["GATA Name"] = filtered_df["Gene ID"].map(gene_mapping)
        

        base_name = ""
        if file.startswith("results_Deseq2_"):
            base_name = file[len("results_Deseq2_"):]
        elif file.startswith("results_Deseq2"):
            base_name = file[len("results_Deseq2"):]
        else:
            base_name = file
        base_name = base_name.replace(".tsv", "").strip("_ ")
        

        filtered_df["Base Name"] = base_name
        

        required_columns = ["Base Name", "Gene ID", "GATA Name", "p-value", "log2(Exp/Control)"]
        missing_columns = [col for col in required_columns if col not in filtered_df.columns]
        if missing_columns:
            log(f"В файле {file} отсутствуют колонки: {', '.join(missing_columns)}. Пропускаю этот файл.")
            continue
        
        filtered_df = filtered_df[required_columns]
        dfs.append(filtered_df)
    
    if not dfs:
        log("Нет данных для объединения после фильтрации по интересующим генам.")
        raise ValueError("No data extracted from results files")
    

    combined_df = pd.concat(dfs, ignore_index=True)


    combined_df["GATA_num"] = combined_df["GATA Name"].apply(extract_gata_number)
    combined_df = combined_df.sort_values(by=["Base Name", "GATA_num"])
    combined_df = combined_df.drop(columns=["GATA_num"])
    
    output_file = os.path.join(output_folder, "Deseq2.txt")
    combined_df.to_csv(output_file, sep="\t", index=False)
    log(f"Итоговый файл сохранён: {output_file}")

if __name__ == "__main__":
    output_folder, gene_mapping = load_settings()
    extract_genes(output_folder, gene_mapping)
