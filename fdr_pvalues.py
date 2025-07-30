import os
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind


input_file = "/mnt/c/Sra_tool/output/GTF_results_fpkm_all.txt"
output_folder = "/mnt/c/Sra_tool/output"
debug_log_file = os.path.join(output_folder, "debug_log.txt")


def get_unique_filename(base_name, extension, folder):
    counter = 1
    file_name = f"{base_name}{extension}"
    while os.path.exists(os.path.join(folder, file_name)):
        file_name = f"{base_name}_{counter}{extension}"
        counter += 1
    return os.path.join(folder, file_name)


output_file = get_unique_filename("GTF_results_fdr_pvalues", ".txt", output_folder)


with open(debug_log_file, "w") as log:
    log.write("Начинаем обработку данных...\n")


    df = pd.read_csv(input_file, sep="\t")
    log.write(f"Загружено {len(df)} строк из {input_file}\n")


    def get_base_name(file_name):
        return file_name.replace("_merged.gtf", "").rstrip("0123456789")


    df["Base Name"] = df["File"].apply(get_base_name)


    experimental_groups = df[~df["Base Name"].str.contains("Control")]
    control_groups = df[df["Base Name"].str.contains("Control")]

    log.write(f"Экспериментов: {len(experimental_groups)}\n")
    log.write(f"Контролей: {len(control_groups)}\n")


    grouped_exp = experimental_groups.groupby(["Base Name", "Gene ID", "GATA Name"])["FPKM"].apply(list).reset_index()
    grouped_ctrl = control_groups.groupby(["Base Name", "Gene ID", "GATA Name"])["FPKM"].apply(list).reset_index()


    grouped_ctrl["Base Name"] = grouped_ctrl["Base Name"].str.replace("Control", "", regex=False)


    missing_controls = set(grouped_exp["Base Name"]) - set(grouped_ctrl["Base Name"])
    log.write(f"Пропущенные контрольные группы: {missing_controls}\n")


    merged_df = pd.merge(grouped_exp, grouped_ctrl, on=["Base Name", "Gene ID", "GATA Name"], suffixes=("_exp", "_ctrl"), how="inner")

    log.write(f"После объединения {len(merged_df)} записей\n")


    merged_df["FPKM_exp"] = merged_df["FPKM_exp"].apply(lambda x: list(map(float, x)) if isinstance(x, list) else [])
    merged_df["FPKM_ctrl"] = merged_df["FPKM_ctrl"].apply(lambda x: list(map(float, x)) if isinstance(x, list) else [])


    p_values = []
    for _, row in merged_df.iterrows():
        exp_values = row["FPKM_exp"]
        ctrl_values = row["FPKM_ctrl"]

        if len(exp_values) > 1 and len(ctrl_values) > 1:  
            p_val = ttest_ind(exp_values, ctrl_values, equal_var=False)[1]
        else:  
            p_val = 1.0

        p_values.append(p_val)

    merged_df["p-value"] = p_values


    m = (merged_df["p-value"] < 1.0).sum()  

    log.write(f"Всего независимых тестов (m) = {m}\n")


    p_values_sorted_idx = np.argsort(merged_df["p-value"].values)
    q_values = np.zeros(len(p_values))

    for rank, i in enumerate(p_values_sorted_idx):
        q_values[i] = (merged_df["p-value"].iloc[i] * m) / (rank + 1)


    for i in range(len(q_values) - 2, -1, -1):
        q_values[i] = min(q_values[i], q_values[i + 1])


    merged_df["FDR"] = q_values


    gata_order = {
        "GATA-1": 1, "GATA-2": 2, "GATA-3": 3, "GATA-4": 4,
        "GATA-5": 5, "GATA-6": 6, "GATA-7": 7, "GATA-8": 8,
        "GATA-9": 9, "GATA-10": 10, "GATA-11": 11, "GATA-12": 12
    }
    merged_df["GATA Order"] = merged_df["GATA Name"].map(gata_order)
    merged_df = merged_df.sort_values(by=["GATA Order"]).drop(columns=["GATA Order"])


    merged_df["FDR"] = merged_df["FDR"].fillna(1.0)


    final_df = merged_df[["Base Name", "Gene ID", "GATA Name", "FDR"]]


    final_df.to_csv(output_file, sep="\t", index=False)
    log.write(f"FDR-корректированные p-values сохранены в {output_file}\n")

print(f"Код выполнен. Проверьте лог в {debug_log_file}")
