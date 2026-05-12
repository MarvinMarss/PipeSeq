import os
import pandas as pd
from statsmodels.stats.multitest import multipletests


input_file     = "/mnt/c/Sra_tool/output/GTF_results_pvalues.txt"
output_folder  = "/mnt/c/Sra_tool/output"
debug_log_file = os.path.join(output_folder, "debug_log.txt")


def get_unique_filename(base_name, extension, folder):
    counter    = 1
    file_name  = f"{base_name}{extension}"
    while os.path.exists(os.path.join(folder, file_name)):
        file_name = f"{base_name}_{counter}{extension}"
        counter  += 1
    return os.path.join(folder, file_name)

output_file = get_unique_filename("GTF_results_fdr_pvalues", ".txt", output_folder)


with open(debug_log_file, "w", encoding="utf-8") as log:
    log.write("=== Начинаем корректировку FDR ===\n")

    df = pd.read_csv(input_file, sep="\t")
    log.write(f"Загружено {len(df)} строк из {input_file}\n")


    df["p-value"] = df["p-value"].astype(float)


    _, qvals, _, _ = multipletests(df["p-value"].values,
                                   alpha=0.05,
                                   method="fdr_bh")
    df["FDR"] = qvals
    log.write("Выполнена FDR-коррекция методом Бенджамини–Хохберга\n")


    gata_order = {
        "GATA-1": 1,  "GATA-2": 2,  "GATA-3": 3,  "GATA-4": 4,
        "GATA-5": 5,  "GATA-6": 6,  "GATA-7": 7,  "GATA-8": 8,
        "GATA-9": 9,  "GATA-10": 10,"GATA-11": 11,"GATA-12": 12
    }
    df["GATA Order"] = df["GATA Name"].map(gata_order)
    df = df.sort_values(by=["GATA Order"]).drop(columns=["GATA Order"])
    log.write("Отсортировано по порядку GATA\n")


    df[["Base Name", "Gene ID", "GATA Name", "FDR"]].to_csv(
        output_file, sep="\t", index=False
    )
    log.write(f"FDR-показатели сохранены в {output_file}\n")
    log.write("=== Завершено ===\n")

print(f"Готово. Логи в {debug_log_file}")
