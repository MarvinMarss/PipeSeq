Установка окружения на новом ПК
Перед первым запуском необходимо установить все зависимости.

1. Открой PowerShell от имени администратора
Нажми Win + X → выбери Windows PowerShell (Admin) или Терминал (Admin).

2. Установи WSL и Ubuntu
Вставь в PowerShell:

winget install -e --id Python.Python.3.12

Перезагрузи терминал 

pip install PyQt6

pip install pandas numpy scipy seaborn matplotlib pyDESeq2

pip install pywin32


wsl --install -d Ubuntu
Дождись завершения установки. Перезагрузи ПК, если попросят.

3. После перезагрузки — снова PowerShell (Admin):
wsl --install -d Ubuntu
wsl
Откроется Ubuntu. Дождись приглашения ввести имя пользователя и пароль (введи любые, которые запомнишь).

4. Внутри Ubuntu — настрой рабочее окружение:

# Обновление системы
sudo apt update && sudo apt upgrade -y

sudo add-apt-repository universe

sudo add-apt-repository multiverse

sudo apt update

# Установка необходимых пакетов и инструментов
sudo apt install -y python3-venv build-essential zlib1g-dev \
  libbz2-dev liblzma-dev libncurses-dev \
  libcurl4-openssl-dev libssl-dev libsqlite3-dev wget curl \
  git unzip samtools hisat2 stringtie libgl1 libxkbcommon-x11-0

# Создание изолированного Python-окружения - перезагрузите терминал 

python3 -m venv ~/pipeseq_env
source ~/pipeseq_env/bin/activate

# Обновление pip и установка нужных Python-библиотек
pip install --upgrade pip
pip install pandas numpy scipy seaborn matplotlib pyqt6 pyDESeq2

5. Запуск пайплайна
После завершения всех шагов:

Дважды щёлкни PipeSeq.bat — это запустит оболочку пайплайна.

6. После запуска на новом устройстве, всегда сначала задай расположение папок внутри PipeSeq.
7. Для корректной работы программ скачай и помести в папку Genome файл .fa одной версии генома объекта интереса и его .gtf аннотации. В случае ошибок с анотированием включить внутри PipeSeq gtf.fix.
8. Можно приступать к работе, записав в окно 1 ID эксперимента-ncbi-имя, которое ты ему присвоишь, в следующем формате:
Копировать
SRX8380271-HighLight1; SRX8380270-HighLight2; SRX8380269-HighLight3; SRX5120532-HighLightControl1; SRX5120531-HighLightControl2; SRX5120530-HighLightControl3
Номер повторности пишется в конце, как в примере.

Или можешь использовать локальные файлы sra, присвоив им перед этим имена в той же форме.

## Структура проекта

```bash
PipeSeq/
├── align_hisat2.py             # Выравнивание FASTQ -> SAM (HISAT2)
├── process_sam_to_bam.py       # SAM -> BAM + сортировка
├── stringtie_expression.py     # Анализ экспрессии (StringTie)
├── extract_fpkm.py             # Извлечение FPKM + логарифм
├── pvalues_log2.py             # Объединение результатов
├── GTF_results_pvalues.py      # Альтернатива для расчёта p-values
├── temp_card_p.py              # Визуализация (heatmap)
├── fix.gtf.py                  # (опционально) корректировка GTF
├── run_gui.py                  # GUI-интерфейс (PyQt6)
├── settings.json               # Конфигурация путей и параметров
└── pipeline_log.txt            # Лог выполнения
```

---

## Настройки `settings.json`

```json
{
  "folders": {
    "fastq_folder": "путь к FASTQ-файлам",
    "bam_folder": "путь к BAM-выходу",
    "gtf_folder": "путь к GTF от StringTie",
    "results_folder": "путь к итоговым таблицам",
    "genome_folder": "путь к .fa и .gtf",
    "genome_index": "путь к индексу HISAT2"
  },
  "options": {
    "delete_intermediate_files": true,
    "use_fdr_correction": true,
    "fix_genome": false
  },
  "stringtie": {
    "coverage_cutoff": 0.01
  },
  "gene_mapping": {
    "CHLRE_01g025050v5": "GATA-1",
    "CHLRE_10g435450v5": "GATA-2"
  },
  "visualization": {
    "show_p_values": true
  }
}
```

---

## Запуск пайплайна

```bash
python run_gui.py
```

Или запуск по шагам:

```bash
python align_hisat2.py
python process_sam_to_bam.py
python stringtie_expression.py
python extract_fpkm.py
python pvalues_log2.py
```

---

##  Особенности

- Поддержка paired-end и одиночных FASTQ-файлов.
- Умная обработка: сортировка, пропуск лишних шагов, удаление временных файлов.
- Геномный индекс строится автоматически при необходимости.
- Поддержка настройки чувствительности StringTie (`-c`) из GUI.
- Прозрачная визуализация и сохранение логов.

---

## Зависимости

- Python 3.10+
- PyQt6
- HISAT2 (в WSL)
- SAMtools (в WSL)
- StringTie (в WSL)

---

##  Примечания

- HISAT2 использует базовое имя `genome_index`, и ожидает 8 файлов: `.1.ht2`, ..., `.8.ht2`
- Все пути внутри WSL автоматически трансформируются (`/mnt/c/...`)
- Файл GTF можно поправить скриптом `fix.gtf.py`, если есть ошибки аннотации.


Удаление окружения на ПК
Для полного удаления всех установленных программ и настроек следуй этим шагам.

1. Открой PowerShell от имени администратора
Нажми Win + X → выбери Windows PowerShell (Admin) или Терминал (Admin).

2. Удаление WSL и Ubuntu
Вставь в PowerShell:

powershell

# Удаление установленного Python
winget uninstall Python.Python.3.12

# Удаление всех Python-пакетов
pip uninstall -y PyQt6 pandas numpy scipy seaborn matplotlib pyDESeq2

# Удаление WSL
wsl --uninstall
Дождись завершения удаления.

3. После этого, чтобы удалить все настройки и окружения Ubuntu, снова открой PowerShell (Admin):
powershell

powershell
Убедись, что Ubuntu полностью удалена.

4. Удаление Python-окружения в Ubuntu:
Если ты использовал виртуальное окружение Python, выполни следующие команды:

bash

# Деактивация окружения
deactivate

# Удаление изолированного окружения Python
rm -rf ~/pipeseq_env
5. Удаление всех установленных пакетов и инструментов в Ubuntu:
bash

# Удаление всех установленных пакетов
sudo apt purge -y python3-venv build-essential zlib1g-dev \
  libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev \
  libcurl4-openssl-dev libssl-dev libsqlite3-dev wget curl \
  git unzip samtools hisat2 stringtie libgl1 libxkbcommon-x11-0

# Очистка системы от ненужных файлов
sudo apt autoremove -y
sudo apt clean
6. Удаление файлов и настроек PipeSeq
Чтобы полностью удалить все, что связано с пайплайном, удалите следующие папки и файлы:

Удалите папку с геномами и аннотациями, если она была скачана.

Удалите файл settings.json, если хотите сбросить все настройки:

bash

rm -rf ~/pipeSeq/settings.json
Удалите все промежуточные файлы (например, файлы в папках с результатами и логами):

bash

rm -rf ~/pipeSeq/results
rm -rf ~/pipeSeq/logs
7. Удаление папки с программой PipeSeq (если она была установлена локально):
bash

rm -rf ~/pipeSeq
8. Удаление установленных SRA инструментов:
Если установил SRA Toolkit, можно выполнить команду для его удаления:

bash

winget uninstall SRA.SRA-Toolkit
9. Очистка среды
Если ты не собираешься использовать WSL и Ubuntu больше:

powershell

# Удаление Ubuntu из Windows
wsl --unregister Ubuntu

## Поддержка

alexnerezenko@gmail.com
