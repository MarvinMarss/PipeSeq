import os
import re
import subprocess
import json
import sys
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


import tkinter as tk
from tkinter import messagebox

SETTINGS_FILE = "settings.json"
LOG_FILE = "deseq2_analysis_log.txt"


def log(message):
    print(message)
    with open(LOG_FILE, "a", encoding="utf-8") as log_file:
        log_file.write(message + "\n")


def load_settings():
    if not os.path.exists(SETTINGS_FILE):
        log(f"Файл настроек {SETTINGS_FILE} не найден!")
        sys.exit(1)

    with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
        settings = json.load(f)

    folders = settings.get("folders", {})
    bam_folder = folders.get("bam_folder")
    genome_folder = folders.get("genome_folder")
    results_folder = folders.get("results_folder")

    if not bam_folder or not genome_folder or not results_folder:
        log("Не заданы bam_folder, genome_folder или results_folder в settings.json!")
        sys.exit(1)

    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    return bam_folder, genome_folder, results_folder, settings


def convert_path_to_wsl(win_path):
    abs_path = os.path.abspath(win_path)
    path = abs_path.replace("\\", "/")
    if ":" in path:
        drive, rest = path.split(":", 1)
        return f"/mnt/{drive.lower()}{rest}"
    return path


def run_command(command):
    log(f"Запуск в WSL:\n{command}")
    try:
        subprocess.run(f'wsl bash -c "{command}"', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        log(f"Ошибка выполнения команды: {e}")
        sys.exit(1)



def parse_sample_name(filename):

    full_sample = filename

    base = filename
    for suffix in ["_paired_sorted", "_single_sorted"]:
        if base.endswith(suffix):
            base = base[:-len(suffix)]

    m = re.match(r"^(.*?)[Cc]ontrol(\d*)$", base)
    if m:
        experiment = m.group(1)
        replicate = int(m.group(2)) if m.group(2) != "" else None
        condition = "control"
        sample = base  
        return {"sample": sample, "full_sample": full_sample, "experiment": experiment, "replicate": replicate, "condition": condition}
    else:

        m = re.match(r"^(.*?)(\d+)$", base)
        if m:
            experiment = m.group(1)
            replicate = int(m.group(2))
            condition = "treated"
            sample = base  
            return {"sample": sample, "full_sample": full_sample, "experiment": experiment, "replicate": replicate, "condition": condition}
        else:
            experiment = base
            condition = "treated"
            return {"sample": base, "full_sample": full_sample, "experiment": experiment, "replicate": None, "condition": condition}



def parse_sample_info(bam_files):
    samples = []
    for bam in bam_files:
        info = parse_sample_name(os.path.splitext(bam)[0])
        info["bam_file"] = bam
        samples.append(info)
    return pd.DataFrame(samples)



def run_featurecounts_individual(bam_file, bam_folder, results_folder, annotation_gtf_wsl, extra_options=""):
    full_sample = os.path.splitext(bam_file)[0]
    output_counts = os.path.join(results_folder, f"gene_counts_{full_sample}.txt")
    output_counts_wsl = convert_path_to_wsl(output_counts)
    bam_path = os.path.join(bam_folder, bam_file)
    bam_path_wsl = convert_path_to_wsl(bam_path)
    cmd = (
        f"featureCounts -a {annotation_gtf_wsl} -o {output_counts_wsl} "
        f"{extra_options} -T 4 -g gene_id -t exon -s 0 {bam_path_wsl}"
    )
    run_command(cmd)
    return output_counts


def prepare_count_matrix(count_file):
    df = pd.read_csv(count_file, sep='\t', comment='#', skiprows=1)
    df = df.drop(columns=["Chr", "Start", "End", "Strand", "Length"])
    df = df.rename(columns={"Geneid": "gene"})
    return df.set_index("gene")



def run_global_deseq2(count_matrix, sample_table):
    dds = DeseqDataSet(
        counts=count_matrix.T,
        metadata=sample_table,
        design_factors="group",
        refit_cooks=True,
        n_cpus=4
    )
    dds.deseq2()
    return dds


def format_and_save_results(results, mapping, results_folder, base_name, output_filename):
    results = results.reset_index()
    results = results.rename(columns={
        "log2FoldChange": "log2(Exp/Control)",
        "padj": "FDR",
        "pvalue": "p-value",
        "gene": "Gene ID"
    })
    results["GATA Name"] = results["Gene ID"].map(mapping)
    results["Base Name"] = base_name
    cols = ["Base Name", "Gene ID", "GATA Name", "p-value", "log2(Exp/Control)"]
    output_path = os.path.join(results_folder, output_filename)
    results[cols].to_csv(output_path, sep="\t", index=False)
    log(f"Результаты сохранены в файл: {output_path}")


def move_counts_files(results_folder):
    counts_folder = os.path.join(results_folder, "Counts")
    if not os.path.exists(counts_folder):
        os.makedirs(counts_folder)
        log(f"Создана папка для counts: {counts_folder}")

    for file in os.listdir(results_folder):
        if file.startswith("gene_counts_") and (file.endswith(".txt") or file.endswith(".txt.summary")):
            src = os.path.join(results_folder, file)
            dst = os.path.join(counts_folder, file)
            try:
                os.rename(src, dst)
                log(f"Перемещён файл {file} в папку Counts")
            except Exception as e:
                log(f"Не удалось переместить файл {file}: {e}")


def main():
    bam_folder, genome_folder, results_folder, settings = load_settings()
    counts_folder = os.path.join(results_folder, "Counts")
    

    skip_counting = False
    if os.path.exists(counts_folder):
        count_files = [f for f in os.listdir(counts_folder) if f.startswith("gene_counts_") and (f.endswith(".txt") or f.endswith(".txt.summary"))]
        if count_files:
            root = tk.Tk()
            root.withdraw()  
            answer = messagebox.askyesno("Подсчёт файлов", "Найдены файлы подсчёта в папке Counts.\nПропустить этап подсчёта и использовать существующие файлы?")
            if answer:
                skip_counting = True
            root.destroy()


    all_bam_files = [f for f in os.listdir(bam_folder) if f.endswith("_sorted.bam")]
    if not all_bam_files:
        log("Нет _sorted.bam файлов в папке!")
        sys.exit(1)

    log(f"Найдено {len(all_bam_files)} BAM-файлов.")


    annotation_gtf = None
    for file in os.listdir(genome_folder):
        if file.endswith(".gtf"):
            annotation_gtf = os.path.join(genome_folder, file)
            break
    if not annotation_gtf:
        log("GTF-файл аннотации не найден!")
        sys.exit(1)
    annotation_gtf_wsl = convert_path_to_wsl(annotation_gtf)


    sample_df = parse_sample_info(all_bam_files)
    sample_df["group"] = sample_df["experiment"] + "_" + sample_df["condition"]


    count_file_dict = {}
    if not skip_counting:
        for idx, row in sample_df.iterrows():
            bam_file = row["bam_file"]
            extra_options = "-p" if "_paired" in bam_file.lower() else ""
            output_file = run_featurecounts_individual(bam_file, bam_folder, results_folder, annotation_gtf_wsl, extra_options)
            count_file_dict[row["full_sample"]] = output_file
    else:
        log("Пропускаем этап подсчёта. Используем файлы из папки Counts.")
        for idx, row in sample_df.iterrows():
            candidate = os.path.join(counts_folder, f"gene_counts_{row['full_sample']}.txt")
            if os.path.exists(candidate):
                count_file_dict[row["full_sample"]] = candidate
            else:
                log(f"Файл подсчёта для образца {row['full_sample']} не найден в папке Counts.")
                sys.exit(1)


    count_dfs = []
    for sample in sample_df["full_sample"]:
        count_file = count_file_dict[sample]
        df = prepare_count_matrix(count_file)
        original_col = df.columns[0]
        df = df.rename(columns={original_col: sample})
        count_dfs.append(df)
    global_counts = pd.concat(count_dfs, axis=1)
    global_counts_filename = os.path.join(results_folder, "global_merged_counts.tsv")
    global_counts.to_csv(global_counts_filename, sep="\t")
    log(f"Глобальные count данные сохранены в файл: {global_counts_filename}")


    sample_table = sample_df.set_index("full_sample")

    log("Запуск глобальной нормализации DESeq2...")
    dds = run_global_deseq2(global_counts, sample_table)

    gene_mapping = settings.get("gene_mapping", {})


    for experiment in sample_df["experiment"].unique():
        log(f"Извлечение результата для эксперимента: {experiment}")
        group_treated = f"{experiment}_treated"
        group_control = f"{experiment}_control"

        exp_samples = sample_table[sample_table["experiment"] == experiment]
        if group_treated not in exp_samples["group"].values or group_control not in exp_samples["group"].values:
            log(f"Для эксперимента {experiment} отсутствует один из классов (treated или control). Пропускаю.")
            continue

        cond_counts = exp_samples["condition"].value_counts()
        insufficient_reps = cond_counts.min() < 3

        stat_res = DeseqStats(dds, contrast=("group", group_treated, group_control))
        stat_res.summary()
        res = stat_res.results_df

        if insufficient_reps:
            log(f"В эксперименте {experiment} менее 3 биологических повторов. Устанавливаю p-value и FDR = 1.")
            res["pvalue"] = 1
            res["padj"] = 1

        base_name = experiment
        output_filename = f"results_Deseq2_{experiment}.tsv"
        format_and_save_results(res, gene_mapping, results_folder, base_name, output_filename)

    log("Глобальный анализ DESeq2 завершён для всех экспериментов.")

    if not skip_counting:
        move_counts_files(results_folder)


if __name__ == "__main__":
    main()
