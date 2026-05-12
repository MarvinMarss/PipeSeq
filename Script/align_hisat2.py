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

def find_fasta_file(folder):
    for file in os.listdir(folder):
        if file.endswith(".fa") or file.endswith(".fasta"):
            fasta_path = os.path.join(folder, file)
            log(f"Найден FASTA файл: {fasta_path}")
            return fasta_path
    return None

def build_hisat2_index(fasta_file, index_base_path):
    fasta_wsl = convert_path_to_wsl(fasta_file)
    index_wsl = convert_path_to_wsl(index_base_path)
    command = f"hisat2-build {fasta_wsl} {index_wsl}"
    log(f"Строим индекс HISAT2 в папке индекса:\n{command}")
    try:
        subprocess.run(f'wsl bash -c "{command}"', shell=True, check=True)
        log(f"Индекс HISAT2 успешно создан: {index_base_path}")
    except subprocess.CalledProcessError as e:
        log(f"Ошибка при создании индекса: {e}")
        sys.exit(1)

def check_or_create_hisat2_index(genome_folder, genome_index_folder):
    required_extensions = [f".{i}.ht2" for i in range(1, 9)]
    base_name = "genome_index"
    index_base_path = os.path.join(genome_index_folder, base_name)
    missing_files = []
    for ext in required_extensions:
        index_file_path = os.path.join(genome_index_folder, f"{base_name}{ext}")
        if not os.path.exists(index_file_path):
            missing_files.append(f"{base_name}{ext}")
    if not missing_files:
        log(f"Индекс HISAT2 уже существует: {index_base_path}")
        return index_base_path
    log(f"Индекс не найден или не полный. Будем строить новый.")
    fasta_file = find_fasta_file(genome_folder)
    if not fasta_file:
        log(f"FASTA файл (.fa/.fasta) не найден в папке {genome_folder}!")
        sys.exit(1)
    build_hisat2_index(fasta_file, index_base_path)
    return index_base_path

def load_settings():
    if not os.path.exists(SETTINGS_FILE):
        log(f"Файл настроек {SETTINGS_FILE} не найден!")
        sys.exit(1)
    with open(SETTINGS_FILE, "r", encoding="utf-8") as f:
        settings = json.load(f)
    folders = settings.get("folders", {})
    fastq_folder = folders.get("fastq_folder")
    bam_folder = folders.get("bam_folder")
    genome_folder = folders.get("genome_folder")
    genome_index_folder = folders.get("genome_index")
    if not fastq_folder or not bam_folder or not genome_folder or not genome_index_folder:
        log("В settings.json должны быть пути: fastq_folder, bam_folder, genome_folder и genome_index!")
        sys.exit(1)
    index_base = check_or_create_hisat2_index(genome_folder, genome_index_folder)
    return fastq_folder, bam_folder, index_base

def align_with_hisat2():

    sample_filter = None
    if len(sys.argv) > 1:
        sample_filter = sys.argv[1].strip()
        log(f"Фильтр по образцу: {sample_filter}")

    fastq_folder, bam_folder, genome_index_base = load_settings()
    if not os.path.exists(fastq_folder):
        log(f"Папка FASTQ не найдена: {fastq_folder}")
        sys.exit(1)
    if not os.path.exists(bam_folder):
        log(f"Папка BAM не найдена: {bam_folder}")
        sys.exit(1)
    fastq_files = [f for f in os.listdir(fastq_folder) if f.endswith(".fastq")]
    if sample_filter:

        fastq_files = [f for f in fastq_files if f.startswith(sample_filter)]
    if not fastq_files:
        log("Нет файлов .fastq для обработки по заданному фильтру.")
        sys.exit(1)
    fastq_folder_wsl = convert_path_to_wsl(fastq_folder)
    bam_folder_wsl = convert_path_to_wsl(bam_folder)
    genome_index_wsl = convert_path_to_wsl(genome_index_base)
    paired_files = {}
    single_files = []
    for f in fastq_files:
        if "_1.fastq" in f:
            base = f.replace("_1.fastq", "")
            pair = f.replace("_1.fastq", "_2.fastq")
            if pair in fastq_files:
                paired_files[base] = (f, pair)
            else:
                single_files.append(f)
        elif "_2.fastq" not in f:
            single_files.append(f)
    log(f"Найдено {len(paired_files)} paired-end и {len(single_files)} одиночных fastq файлов для обработки.")

    for base, (r1, r2) in paired_files.items():
        r1_path = f"{fastq_folder_wsl}/{r1}"
        r2_path = f"{fastq_folder_wsl}/{r2}"
        output_sam = f"{bam_folder_wsl}/{base}_paired.sam"
        log(f"\nВыравнивание парных файлов: {r1} + {r2}")
        command = (
            f"hisat2 -p 4 -x {genome_index_wsl} -1 {r1_path} -2 {r2_path} -S {output_sam}"
        )
        run_command(command)

    for f in single_files:
        sample_name = os.path.splitext(f)[0]
        fastq_path = f"{fastq_folder_wsl}/{f}"
        output_sam = f"{bam_folder_wsl}/{sample_name}_single.sam"
        log(f"\nВыравнивание одиночного файла: {f}")
        command = (
            f"hisat2 -p 4 -x {genome_index_wsl} -U {fastq_path} -S {output_sam}"
        )
        run_command(command)
    log("\nВыравнивание всех файлов завершено!")

if __name__ == "__main__":
    align_with_hisat2()
