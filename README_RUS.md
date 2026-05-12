# PipeSeq

PipeSeq - локальный Windows/WSL-пайплайн с графическим интерфейсом для обработки RNA-seq данных.
Он объединяет подготовку SRA/FASTQ, выравнивание HISAT2, обработку SAM/BAM, анализ StringTie,
featureCounts/DESeq2-статистику и визуализацию результатов.

Английская документация: [README.md](README.md)

## Что хранится в git

В репозитории хранятся код, launcher, документация, примеры конфигурации и небольшие статические
файлы. Большие рабочие данные специально не добавляются в git:

- сырые FASTQ/FQ/SRA файлы;
- SAM/BAM результаты выравнивания;
- FASTA геномы, GTF/GFF аннотации и HISAT2 индексы;
- сгенерированные результаты и логи;
- локальные `Script/settings.json` и `Script/Mind.json`;
- скачанные бинарные сборки FastQC и SRA Toolkit.

При этом локальная структура папок сохранена: `Fastq/`, `Sra/`, `Genome/`, `Genome/Index/`,
`Output/`, `GTF/`, `Results/` и `Tool/` представлены README-файлами.

## Структура репозитория

```text
PipeSeq/
├── PipeSeq.bat                 # Windows-запуск GUI
├── Script/
│   ├── PipeSeq.py              # Главный PyQt6 GUI
│   ├── align_hisat2.py         # FASTQ -> SAM через HISAT2
│   ├── process_sam_to_bam.py   # SAM -> отсортированный BAM
│   ├── stringtie_expression.py # Шаг StringTie
│   ├── deseq2_analysis.py      # featureCounts + PyDESeq2
│   ├── extract_*.py            # Извлечение таблиц экспрессии
│   ├── pvalues_*.py            # Статистические расчеты
│   ├── temp_card_p.py          # Heatmap/визуализация
│   ├── settings.example.json   # Пример конфигурации
│   └── Mind.example.json       # Пример файла памяти GUI
├── Fastq/                      # Локальные FASTQ, игнорируются git
├── Sra/                        # Локальные SRA, игнорируются git
├── Genome/                     # Геном/аннотация, игнорируются git
├── Output/                     # SAM/BAM output, игнорируется git
├── GTF/                        # StringTie output, игнорируется git
├── Results/                    # Итоговые результаты, игнорируются git
├── Tool/                       # Локальные сторонние инструменты, игнорируются git
└── docs/                       # Заметки для сопровождения
```

## Требования

- Windows с PowerShell.
- Python 3.10 или новее.
- WSL с Ubuntu для Linux-инструментов.
- Python-пакеты из `requirements.txt`.
- Инструменты в WSL: HISAT2, SAMtools, StringTie и Subread/featureCounts.
- Опционально локально в `Tool/`: FastQC и SRA Toolkit.

## Установка

Открой PowerShell от имени администратора и установи Python:

```powershell
winget install -e --id Python.Python.3.12
```

Создай виртуальное окружение из корня репозитория:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Установи Ubuntu в WSL:

```powershell
wsl --install -d Ubuntu
```

Внутри Ubuntu установи нужные bioinformatics-инструменты:

```bash
sudo apt update
sudo apt install -y git unzip wget curl hisat2 samtools stringtie subread
```

## Подготовка локальных данных

1. Помести FASTQ-файлы в `Fastq/` или локальные SRA-файлы в `Sra/`.
2. Помести FASTA геном и соответствующую GTF/GFF аннотацию в `Genome/`.
3. Помести готовый HISAT2 index в `Genome/Index/` или дай пайплайну создать его, если это доступно.
4. Локальные Windows-сборки FastQC/SRA Toolkit можно положить в `Tool/` или указать путь к ним в GUI.

При первом запуске PipeSeq создаст `Script/settings.json`, если файла нет. Также можно скопировать
`Script/settings.example.json` в `Script/settings.json` и вручную заменить пути.

## Запуск

В Windows дважды щелкни:

```text
PipeSeq.bat
```

Или запусти GUI напрямую:

```powershell
python Script\PipeSeq.py
```

`PipeSeq.bat` перед запуском переходит в папку `Script/`, потому что скрипты пайплайна ожидают
конфигурационные файлы именно там.

## Формат ID экспериментов

Вводи ID экспериментов и имена образцов в GUI в таком формате:

```text
SRX8380271-HighLight1; SRX8380270-HighLight2; SRX8380269-HighLight3; SRX5120532-HighLightControl1; SRX5120531-HighLightControl2; SRX5120530-HighLightControl3
```

Номер повторности лучше ставить в конце имени образца. Локальные `.sra` файлы можно переименовывать
по той же схеме.

## Проверки для разработки

```powershell
python -m py_compile Script\*.py
```

Опциональная проверка стиля:

```powershell
pip install -r requirements-dev.txt
ruff check Script
```

## Примечания

- HISAT2 index обычно имеет базовое имя `genome_index` и состоит из восьми `.ht2` файлов.
- Windows-пути внутри пайплайна преобразуются в WSL-пути там, где это требуется.
- Политика хранения данных описана в [docs/DATA_AND_TOOLS.md](docs/DATA_AND_TOOLS.md).

## Лицензия

MIT. См. [LICENSE](LICENSE).
