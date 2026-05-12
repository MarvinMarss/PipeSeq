import os
import subprocess
import json
import sys

SETTINGS_FILE = "settings.json"
LOG_FILE = "process_sam_to_bam_log.txt"

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
    delete_intermediate = settings.get("options", {}).get("delete_intermediate_files", False)

    if not bam_folder:
        log("Не задан путь к папке BAM/Output в settings.json!")
        sys.exit(1)

    return bam_folder, delete_intermediate

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

def process_files():
    bam_folder, delete_intermediate = load_settings()

    if not os.path.exists(bam_folder):
        log(f"Папка BAM/Output не найдена: {bam_folder}")
        sys.exit(1)


    sample_filter = None
    if len(sys.argv) > 1:
        sample_filter = sys.argv[1].strip()
        log(f"Фильтр по образцу: {sample_filter}")

    bam_folder_wsl = convert_path_to_wsl(bam_folder)
    files_in_folder = os.listdir(bam_folder)

    if sample_filter:
        sam_files = sorted(f for f in files_in_folder if f.endswith(".sam") and f.startswith(sample_filter))
    else:
        sam_files = sorted(f for f in files_in_folder if f.endswith(".sam"))
    bam_files = []


    if sam_files:
        log(f"Найдено {len(sam_files)} SAM файлов. Конвертация в BAM...")
        for sam_file in sam_files:
            sample_name = os.path.splitext(sam_file)[0]
            sam_path = os.path.join(bam_folder, sam_file)
            sam_path_wsl = f"{bam_folder_wsl}/{sam_file}"
            bam_path_wsl = f"{bam_folder_wsl}/{sample_name}.bam"

            log(f"\nКонвертация {sam_file} -> {sample_name}.bam")
            command = f"samtools view -@ 4 -S -b {sam_path_wsl} -o {bam_path_wsl}"
            run_command(command)

            if delete_intermediate:
                try:
                    os.remove(sam_path)
                    log(f"Удалён исходный SAM: {sam_path}")
                except Exception as e:
                    log(f"Не удалось удалить SAM: {sam_path} - {e}")

            bam_files.append(f"{sample_name}.bam")

    else:
        log("SAM файлов не найдено.")
        if sample_filter:

            bam_files = [f for f in files_in_folder if f.endswith(".bam") and f.startswith(sample_filter) and not f.endswith("_sorted.bam")]
            if not bam_files:
                log("Нет BAM файлов для сортировки по заданному фильтру. Завершаем.")
                sys.exit(0)
            log(f"Найдено {len(bam_files)} BAM файлов для сортировки...")
        else:
            bam_files = [f for f in files_in_folder if f.endswith(".bam") and not f.endswith("_sorted.bam")]
            if not bam_files:
                log("Нет BAM файлов для сортировки. Завершаем.")
                sys.exit(0)
            log(f"Найдено {len(bam_files)} BAM файлов для сортировки...")


    for bam_file in bam_files:
        sample_name = os.path.splitext(bam_file)[0]
        bam_path = os.path.join(bam_folder, bam_file)
        bam_path_wsl = f"{bam_folder_wsl}/{bam_file}"
        sorted_bam_path_wsl = f"{bam_folder_wsl}/{sample_name}_sorted.bam"

        log(f"\nСортировка {bam_file} -> {sample_name}_sorted.bam")
        command = f"samtools sort -@ 4 -o {sorted_bam_path_wsl} {bam_path_wsl}"
        run_command(command)

        if delete_intermediate:
            try:
                os.remove(bam_path)
                log(f"Удалён исходный BAM: {bam_path}")
            except Exception as e:
                log(f"Не удалось удалить BAM: {bam_path} - {e}")

    log("\nОбработка файлов завершена!")

if __name__ == "__main__":
    process_files()
