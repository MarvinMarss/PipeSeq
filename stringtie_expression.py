import os
import subprocess
import json
import sys

SETTINGS_FILE = "settings.json"
LOG_FILE = "stringtie_expression_log.txt"

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
    gtf_folder = folders.get("gtf_folder")
    genome_folder = folders.get("genome_folder")

    if not bam_folder or not gtf_folder or not genome_folder:
        log("Не заданы пути (bam_folder, gtf_folder, genome_folder) в settings.json!")
        sys.exit(1)


    reference_gtf = None
    for file in os.listdir(genome_folder):
        if file.endswith(".gtf"):
            reference_gtf = os.path.join(genome_folder, file)
            break

    if not reference_gtf or not os.path.exists(reference_gtf):
        log(f"GTF-файл генома не найден в папке {genome_folder}")
        sys.exit(1)

    return bam_folder, gtf_folder, reference_gtf, settings

def convert_path_to_wsl(win_path):
    abs_path = os.path.abspath(win_path)
    path = abs_path.replace("\\", "/")
    if ":" in path:
        drive, rest = path.split(":", 1)
        wsl_path = f"/mnt/{drive.lower()}{rest}"
        log(f"Путь {win_path} -> WSL {wsl_path}")
        return wsl_path
    return path

def run_command(command):
    log(f"Запуск в WSL:\n{command}")
    try:
        subprocess.run(f'wsl bash -c "{command}"', shell=True, check=True)
    except subprocess.CalledProcessError as e:
        log(f"Ошибка выполнения команды: {e}")
        sys.exit(1)

def calculate_expression_with_stringtie():
    bam_folder, gtf_target_folder, reference_gtf, settings = load_settings()

    if not os.path.exists(bam_folder):
        log(f"Папка BAM/Output не найдена: {bam_folder}")
        sys.exit(1)

    if not os.path.exists(gtf_target_folder):
        os.makedirs(gtf_target_folder)
        log(f"Создана папка для результатов GTF: {gtf_target_folder}")

    bam_folder_wsl = convert_path_to_wsl(bam_folder)
    gtf_target_folder_wsl = convert_path_to_wsl(gtf_target_folder)
    reference_gtf_wsl = convert_path_to_wsl(reference_gtf)

    sorted_bam_files = [f for f in os.listdir(bam_folder) if f.endswith("_sorted.bam")]

    if not sorted_bam_files:
        log("Нет файлов _sorted.bam для анализа.")
        sys.exit(1)

    log(f"Найдено {len(sorted_bam_files)} файлов _sorted.bam для обработки.")

    options = settings.get("options", {})
    use_strict_annotation = options.get("strict_annotation", False)
    stringtie_c = options.get("stringtie_sensitivity", None)

    stringtie_flags = ""
    if use_strict_annotation:
        stringtie_flags += " -e"
    if stringtie_c is not None:
        stringtie_flags += f" -c {stringtie_c}"

    for bam_file in sorted_bam_files:

        base_name = os.path.splitext(bam_file)[0]
        condition_name = base_name.replace("_paired", "").replace("_single", "")

        bam_path_wsl = f"{bam_folder_wsl}/{bam_file}"
        gtf_output_wsl = f"{gtf_target_folder_wsl}/{condition_name}.gtf"
        coverage_output_wsl = f"{gtf_target_folder_wsl}/{condition_name}_coverage.tsv"

        log(f"Обработка {bam_file} с StringTie...")

        command = (
            f"stringtie {bam_path_wsl} "
            f"-G {reference_gtf_wsl} "
            f"-o {gtf_output_wsl} "
            f"{stringtie_flags} --rf -A {coverage_output_wsl}"
        )
        run_command(command)

        log(f"Завершена обработка: {bam_file}")

    log("Все файлы обработаны StringTie!")

if __name__ == "__main__":
    calculate_expression_with_stringtie()
    log("Скрипт успешно завершён!")
