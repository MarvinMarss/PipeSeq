import os, json, sys
import pandas as pd
from scipy.stats import shapiro, ttest_ind, mannwhitneyu

SETTINGS_FILE = "settings.json"
DEBUG_LOG_FILE = "GTF_results_pvalues_tpm_log.txt"

def log(m):
    print(m)
    with open(DEBUG_LOG_FILE, "a", encoding="utf-8") as f:
        f.write(m + "\n")

def load_results():
    if not os.path.exists(SETTINGS_FILE):
        log("settings.json not found"); sys.exit(1)
    with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
        st = json.load(f)
    rf = st.get("folders", {}).get("results_folder")
    if not rf: log("results_folder not set"); sys.exit(1)
    return rf

def get_unique_filename(base, ext, folder):
    k=1; name=f"{base}{ext}"
    while os.path.exists(os.path.join(folder, name)):
        name=f"{base}_{k}{ext}"; k+=1
    return os.path.join(folder, name)

def find_latest_matching_file(folder, pattern):
    matches = []
    for name in os.listdir(folder):
        if __import__("re").fullmatch(pattern, name):
            path = os.path.join(folder, name)
            matches.append((os.path.getmtime(path), name, path))
    if not matches:
        raise FileNotFoundError(f"Нет файлов, подходящих под шаблон: {pattern}")
    matches.sort()
    return matches[-1][2]

def get_base_name(fn):
    name = fn.replace("_merged.gtf","").replace("_midel_merged.gtf","").replace("_sorted.gtf","")
    return name.rstrip("0123456789")

def main():
    if os.path.exists(DEBUG_LOG_FILE): os.remove(DEBUG_LOG_FILE)
    rf = load_results()
    infile = find_latest_matching_file(rf, r"GTF_results_tpm_all(?:_\d+)?\.txt")
    log(f"Используем входной TPM-файл: {infile}")

    df = pd.read_csv(infile, sep="\t")
    if not set(["File","Gene ID","GATA Name","TPM"]).issubset(df.columns):
        log("Отсутствуют необходимые столбцы (File, Gene ID, GATA Name, TPM)"); sys.exit(1)

    df["Base Name"] = df["File"].apply(get_base_name)
    exp_df  = df[~df["Base Name"].str.contains("Control")]
    ctrl_df = df[ df["Base Name"].str.contains("Control")]

    gexp  = exp_df.groupby(["Base Name","Gene ID","GATA Name"])["TPM"].apply(list).reset_index()
    gctrl = ctrl_df.groupby(["Base Name","Gene ID","GATA Name"])["TPM"].apply(list).reset_index()
    gctrl["Base Name"] = gctrl["Base Name"].str.replace("Control","", regex=False)

    merged = pd.merge(gexp, gctrl, on=["Base Name","Gene ID","GATA Name"],
                      suffixes=("_exp","_ctrl"), how="inner")

    if merged.empty:
        log("Нет пар для расчёта p-value, создаю заглушки.")
        out = pd.DataFrame(df[["Base Name","Gene ID","GATA Name"]].drop_duplicates())
        out["p-value"] = 1.0
    else:
        def pv(row):
            xe = list(map(float, row["TPM_exp"]))
            xc = list(map(float, row["TPM_ctrl"]))
            ne, nc = len(xe), len(xc)
            if ne < 2 or nc < 2: return 1.0
            if ne >= 3 and nc >= 3:
                pe, pc = shapiro(xe)[1], shapiro(xc)[1]
                if pe > 0.05 and pc > 0.05:
                    return ttest_ind(xe, xc, equal_var=False)[1]
                return mannwhitneyu(xe, xc, alternative="two-sided")[1]
            return ttest_ind(xe, xc, equal_var=False)[1]
        merged["p-value"] = merged.apply(pv, axis=1)
        out = merged[["Base Name","Gene ID","GATA Name","p-value"]].drop_duplicates()

    out = out.sort_values(["Base Name","Gene ID","GATA Name"])
    outfile = get_unique_filename("GTF_results_pvalues_tpm", ".txt", rf)
    out.to_csv(outfile, sep="\t", index=False)
    log(f"Готово: {outfile}")

if __name__ == "__main__":
    main()
