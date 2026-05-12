import sys
import os
import subprocess
import time
import json
import glob
import shutil
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QPushButton,
    QTextEdit, QMessageBox, QFileDialog, QProgressBar, QCheckBox,
    QLineEdit, QGroupBox
)
from PyQt6.QtGui import QFont, QPixmap
from PyQt6.QtCore import Qt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MEMORY_FILE = os.path.join(SCRIPT_DIR, "Mind.json")
LOG_FILE = os.path.join(SCRIPT_DIR, "process_log.txt")
PIPELINE_SETTINGS_FILE = os.path.join(SCRIPT_DIR, "settings.json")

def log(message: str) -> None:
    with open(LOG_FILE, 'a', encoding='utf-8') as f:
        f.write(message + '\n')
    print(message)

def normalize_win_path(path: str) -> str:
    return os.path.normpath(path)

class SRAConverterApp(QWidget):
    def __init__(self):
        super().__init__()

        self.sratoolkit_path = ""
        self.sra_download_folder = ""
        self.fastq_output_folder = ""
        self.delete_sra_after_conversion = False
        self.memory = {}
        self.load_memory()

        self.pipeline_settings = {}
        self.load_pipeline_settings()
        self.init_ui()

    # ---------------- I/O настроек ----------------

    def load_memory(self):
        self.memory = {}
        if os.path.exists(MEMORY_FILE):
            try:
                with open(MEMORY_FILE, 'r', encoding='utf-8') as f:
                    data = f.read().strip()
                    if data:
                        self.memory = json.loads(data)
                        self.sratoolkit_path = self.memory.get("sratoolkit_path", "")
                        self.sra_download_folder = self.memory.get("sra_download_folder", "")
                        self.fastq_output_folder = self.memory.get("fastq_output_folder", "")
                        self.delete_sra_after_conversion = self.memory.get("delete_sra_after_conversion", False)
                    else:
                        log("Mind.json is empty. Paths will be chosen manually.")
            except json.JSONDecodeError:
                log("Mind.json is corrupted. Ignoring and continuing with empty paths.")
        else:
            log("Mind.json file not found. Starting with empty paths.")

    def save_memory(self):
        self.memory["sratoolkit_path"] = self.sratoolkit_path
        self.memory["sra_download_folder"] = self.sra_download_folder
        self.memory["fastq_output_folder"] = self.fastq_output_folder
        self.memory["delete_sra_after_conversion"] = self.delete_sra_after_conversion
        with open(MEMORY_FILE, 'w', encoding='utf-8') as f:
            json.dump(self.memory, f, indent=4)

    def load_pipeline_settings(self):
        if os.path.exists(PIPELINE_SETTINGS_FILE):
            try:
                with open(PIPELINE_SETTINGS_FILE, 'r', encoding='utf-8') as f:
                    self.pipeline_settings = json.load(f)
            except json.JSONDecodeError:
                log("settings.json is corrupted. Loading default values.")
                self.pipeline_settings = self.default_pipeline_settings()
                self.save_pipeline_settings()
        else:
            self.pipeline_settings = self.default_pipeline_settings()
            self.save_pipeline_settings()

    def default_pipeline_settings(self):
        return {
            "folders": {
                "fastq_folder": "",
                "bam_folder": "",
                "gtf_folder": "",
                "results_folder": "",
                "genome_folder": "",
                "genome_index": ""
            },
            "options": {
                "delete_intermediate_files": False,
                "fix_genome": False,
                "use_stringtie": True,
                "use_deseq2": False,
                "strict_annotation": False,
                "stringtie_sensitivity": ""
            },
            "gene_mapping": {},
            "visualization": {"show_p_values": True},
            "qc": {
                "enable_fastqc": True,
                "enable_trimming": False,
                "keep_raw_fastq": True,
                "adapter_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                "adapter_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
                "trim_q": 20,
                "trim_minlen_pe": 50,
                "trim_minlen_se": 35
            }
        }

    def save_pipeline_settings(self):
        with open(PIPELINE_SETTINGS_FILE, "w", encoding='utf-8') as f:
            json.dump(self.pipeline_settings, f, indent=4)

    # ---------------- UI ----------------

    def select_folder_dialog(self, title="Select folder", start_path=None):
        folder = QFileDialog.getExistingDirectory(
            self, title, start_path or "C:/", QFileDialog.Option.ShowDirsOnly
        )
        if folder and os.path.exists(folder):
            return normalize_win_path(folder)
        QMessageBox.warning(self, "Selection Error", "Folder not chosen or inaccessible.")
        return ""

    def select_combined_fastq_folder(self):
        default_path = self.fastq_output_folder or self.pipeline_settings["folders"].get("fastq_folder", "")
        folder = self.select_folder_dialog("Select FASTQ folder", default_path)
        if folder:
            self.fastq_output_folder = folder
            self.save_memory()
            self.pipeline_settings["folders"]["fastq_folder"] = folder
            self.save_pipeline_settings()
            if hasattr(self, "btn_fastq"):
                self.btn_fastq.setText(f"FASTQ ({os.path.basename(folder)})")

    def select_pipeline_folder(self, key, label_text):
        folder = self.select_folder_dialog(f"Select folder for {label_text}", self.pipeline_settings["folders"].get(key, ""))
        if folder:
            self.pipeline_settings["folders"][key] = folder
            self.save_pipeline_settings()
            log(f"Folder for {label_text} updated: {folder}")
            button = getattr(self, f"btn_{key}", None)
            if button:
                button.setText(f"{label_text} ({os.path.basename(folder)})")

    def toggle_pipeline_option(self, key, state):
        self.pipeline_settings["options"][key] = bool(state)
        self.save_pipeline_settings()
        log(f"Option {key} set to {bool(state)}")

    def update_pipeline_sensitivity(self):
        value = self.pipeline_sensitivity_input.text().strip()
        try:
            if value:
                float_val = float(value)
                self.pipeline_settings["options"]["stringtie_sensitivity"] = float_val
                self.save_pipeline_settings()
                log(f"📎 StringTie sensitivity (-c) set to {float_val}")
            else:
                if "stringtie_sensitivity" in self.pipeline_settings["options"]:
                    del self.pipeline_settings["options"]["stringtie_sensitivity"]
                    self.save_pipeline_settings()
                    log("📎 Removed -c parameter (sensitivity)")
        except ValueError:
            QMessageBox.warning(self, "Error", "Sensitivity value (-c) must be a number.")

    def edit_pipeline_gene_mapping(self):
        dialog = QWidget()
        dialog.setWindowTitle("Gene List Editor")
        dialog.resize(600, 500)

        main_layout = QVBoxLayout()
        list_widget = QTextEdit()
        list_widget.setReadOnly(True)
        gene_mapping = self.pipeline_settings.get("gene_mapping", {})
        mapping_text = "\n".join([f"{gid} ➜ {name}" for gid, name in gene_mapping.items()])
        list_widget.setText(mapping_text)

        main_layout.addWidget(QLabel("Current gene list:"))
        main_layout.addWidget(list_widget)

        from PyQt6.QtWidgets import QFormLayout
        form_layout = QFormLayout()
        new_gene_id = QLineEdit()
        new_gene_name = QLineEdit()
        form_layout.addRow("New Gene ID:", new_gene_id)
        form_layout.addRow("New Gene Name:", new_gene_name)
        main_layout.addLayout(form_layout)

        buttons_layout = QHBoxLayout()
        add_btn = QPushButton("Add gene")
        def add_gene():
            gid = new_gene_id.text().strip()
            gname = new_gene_name.text().strip()
            if not gid or not gname:
                QMessageBox.warning(dialog, "Error", "Specify both Gene ID and Gene Name.")
                return
            self.pipeline_settings["gene_mapping"][gid] = gname
            self.save_pipeline_settings()
            current_text = list_widget.toPlainText()
            list_widget.setText(current_text + f"\n{gid} ➜ {gname}")
            new_gene_id.clear()
            new_gene_name.clear()
            log(f"Gene added: {gid} ➜ {gname}")
        add_btn.clicked.connect(add_gene)
        remove_btn = QPushButton("Delete gene")
        remove_btn.clicked.connect(lambda: QMessageBox.information(dialog, "Instruction", "To delete a gene, edit the settings.json file manually."))
        buttons_layout.addWidget(add_btn)
        buttons_layout.addWidget(remove_btn)
        main_layout.addLayout(buttons_layout)
        dialog.setLayout(main_layout)
        dialog.show()

    def init_ui(self):
        self.setStyleSheet("""
            QWidget { background-color: #2b2b40; color: #f5f5fa; font-family: Arial; font-size: 14px; }
            QPushButton { background-color: #3b3b5c; border: 1px solid #444; padding: 4px; border-radius: 6px; }
            QPushButton:hover { background-color: #50507a; }
            QTextEdit { background-color: #3b3b5c; color: #f5f5fa; border-radius: 5px; padding: 5px; }
            QLineEdit { background-color: #3b3b5c; color: #f5f5fa; padding: 5px; border-radius: 5px; }
            QLabel { color: #f5f5fa; }
        """)

        main_layout = QVBoxLayout()
        main_layout.setSpacing(8)

        title = QLabel("PipeSeq")
        title.setFont(QFont("Arial", 22, QFont.Weight.Bold))
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        main_layout.addWidget(title)
        cat_image = QLabel()
        pixmap = QPixmap(os.path.join(SCRIPT_DIR, "cat.png"))
        if not pixmap.isNull():
            cat_image.setPixmap(pixmap.scaled(80, 80, Qt.AspectRatioMode.KeepAspectRatio))
            cat_image.setAlignment(Qt.AlignmentFlag.AlignCenter)
            main_layout.addWidget(cat_image)

        group_folders = QGroupBox("Folder Selection")
        grid_folders = QGridLayout()
        grid_folders.setHorizontalSpacing(10)
        grid_folders.setVerticalSpacing(6)
        grid_folders.setContentsMargins(10, 10, 10, 10)

        self.btn_sratoolkit = QPushButton("SRA Toolkit"); self.btn_sratoolkit.clicked.connect(self.select_sratoolkit)
        self.btn_sra = QPushButton("SRA Files"); self.btn_sra.clicked.connect(self.select_sra_folder)
        self.btn_fastq = QPushButton("FASTQ"); self.btn_fastq.clicked.connect(self.select_combined_fastq_folder)
        self.btn_bam = QPushButton("BAM/Output"); self.btn_bam.clicked.connect(lambda: self.select_pipeline_folder("bam_folder", "BAM/Output"))
        self.btn_gtf = QPushButton("GTF"); self.btn_gtf.clicked.connect(lambda: self.select_pipeline_folder("gtf_folder", "GTF"))
        self.btn_results = QPushButton("Results"); self.btn_results.clicked.connect(lambda: self.select_pipeline_folder("results_folder", "Results"))
        self.btn_genome = QPushButton("GENOME"); self.btn_genome.clicked.connect(lambda: self.select_pipeline_folder("genome_folder", "GENOME"))
        self.btn_genome_index = QPushButton("HISAT2 index"); self.btn_genome_index.clicked.connect(lambda: self.select_pipeline_folder("genome_index", "HISAT2 index"))

        grid_folders.addWidget(self.btn_sratoolkit, 0, 0)
        grid_folders.addWidget(self.btn_sra, 0, 1)
        grid_folders.addWidget(self.btn_fastq, 1, 0)
        grid_folders.addWidget(self.btn_bam, 1, 1)
        grid_folders.addWidget(self.btn_gtf, 2, 0)
        grid_folders.addWidget(self.btn_results, 2, 1)
        grid_folders.addWidget(self.btn_genome, 3, 0)
        grid_folders.addWidget(self.btn_genome_index, 3, 1)
        group_folders.setLayout(grid_folders)
        main_layout.addWidget(group_folders)

        group_pipeline = QGroupBox("Pipeline Settings")
        grid_pipeline = QGridLayout()
        grid_pipeline.setHorizontalSpacing(10); grid_pipeline.setVerticalSpacing(6)
        grid_pipeline.setContentsMargins(10, 10, 10, 10)

        chk_delete_intermediate = QCheckBox("Delete intermediate files")
        chk_delete_intermediate.setChecked(self.pipeline_settings["options"].get("delete_intermediate_files", False))
        chk_delete_intermediate.stateChanged.connect(lambda state: self.toggle_pipeline_option("delete_intermediate_files", state))
        chk_fix_genome = QCheckBox("Fix genome (fix.gtf.py)")
        chk_fix_genome.setChecked(self.pipeline_settings["options"].get("fix_genome", False))
        chk_fix_genome.stateChanged.connect(lambda state: self.toggle_pipeline_option("fix_genome", state))
        chk_strict_annotation = QCheckBox("Strict annotation (-e)")
        chk_strict_annotation.setChecked(self.pipeline_settings["options"].get("strict_annotation", False))
        chk_strict_annotation.stateChanged.connect(lambda state: self.toggle_pipeline_option("strict_annotation", state))
        chk_stringtie = QCheckBox("StringTie")
        chk_stringtie.setChecked(self.pipeline_settings["options"].get("use_stringtie", True))
        chk_stringtie.stateChanged.connect(lambda state: self.toggle_pipeline_option("use_stringtie", state))
        chk_deseq2 = QCheckBox("DESeq2")
        chk_deseq2.setChecked(self.pipeline_settings["options"].get("use_deseq2", False))
        chk_deseq2.stateChanged.connect(lambda state: self.toggle_pipeline_option("use_deseq2", state))

        grid_pipeline.addWidget(chk_delete_intermediate, 0, 0)
        grid_pipeline.addWidget(chk_fix_genome, 0, 1)
        grid_pipeline.addWidget(chk_strict_annotation, 1, 0)
        grid_pipeline.addWidget(chk_stringtie, 1, 1)
        grid_pipeline.addWidget(chk_deseq2, 2, 0, 1, 2)

        hbox_sensitivity = QHBoxLayout()
        lbl_sensitivity = QLabel("StringTie (-c):")
        self.pipeline_sensitivity_input = QLineEdit(); self.pipeline_sensitivity_input.setPlaceholderText("For example: 0.1")
        current_s = self.pipeline_settings["options"].get("stringtie_sensitivity", "")
        if current_s: self.pipeline_sensitivity_input.setText(str(current_s))
        self.pipeline_sensitivity_input.editingFinished.connect(self.update_pipeline_sensitivity)
        hbox_sensitivity.addWidget(lbl_sensitivity); hbox_sensitivity.addWidget(self.pipeline_sensitivity_input); hbox_sensitivity.setSpacing(6)

        btn_gene_mapping = QPushButton("Edit gene list"); btn_gene_mapping.clicked.connect(self.edit_pipeline_gene_mapping)

        vbox_pipeline = QVBoxLayout()
        vbox_pipeline.addLayout(grid_pipeline)
        vbox_pipeline.addLayout(hbox_sensitivity)
        vbox_pipeline.addWidget(btn_gene_mapping, alignment=Qt.AlignmentFlag.AlignCenter)
        group_pipeline.setLayout(vbox_pipeline)
        main_layout.addWidget(group_pipeline)

        group_sra = QGroupBox("SRA Conversion")
        vbox_sra = QVBoxLayout(); vbox_sra.setSpacing(6); vbox_sra.setContentsMargins(10, 10, 10, 10)

        lbl_samples = QLabel("Enter pairs SRXID-SampleName\n(format: SRX8380271-HighLight1; SRX5120532-HighLightControl1):")
        vbox_sra.addWidget(lbl_samples)
        self.input_samples = QTextEdit(); self.input_samples.setFixedHeight(60)
        vbox_sra.addWidget(self.input_samples)

        hbox_sra_options = QHBoxLayout()
        chk_delete_sra = QCheckBox("Delete SRA after conversion")
        chk_delete_sra.setChecked(self.delete_sra_after_conversion)
        chk_delete_sra.stateChanged.connect(lambda state: self.toggle_delete_sra(state))
        hbox_sra_options.addWidget(chk_delete_sra)
        self.progress_label = QLabel("Awaiting actions..."); hbox_sra_options.addWidget(self.progress_label, stretch=1)
        vbox_sra.addLayout(hbox_sra_options)

        self.progress_bar = QProgressBar(); self.progress_bar.setValue(0)
        vbox_sra.addWidget(self.progress_bar)

        hbox_sra_buttons = QHBoxLayout()
        btn_download = QPushButton("Download and convert"); btn_download.clicked.connect(self.start_process)
        btn_convert_local = QPushButton("Convert local SRA"); btn_convert_local.clicked.connect(self.convert_existing_sra_files)
        hbox_sra_buttons.addWidget(btn_download); hbox_sra_buttons.addWidget(btn_convert_local)
        vbox_sra.addLayout(hbox_sra_buttons)

        group_sra.setLayout(vbox_sra)
        main_layout.addWidget(group_sra)

        btn_run_pipeline = QPushButton("Run PipeSeq-2"); btn_run_pipeline.clicked.connect(self.run_pipeline)
        main_layout.addWidget(btn_run_pipeline, alignment=Qt.AlignmentFlag.AlignCenter)

        self.setLayout(main_layout)
        self.setWindowTitle("PipeSeq-1")
        self.resize(600, 650)

    # ---------------- события UI ----------------

    def toggle_delete_sra(self, state):
        self.delete_sra_after_conversion = (state == Qt.CheckState.Checked.value)
        self.save_memory()
        log(f"Delete SRA after conversion: {self.delete_sra_after_conversion}")

    def select_sratoolkit(self):
        folder = self.select_folder_dialog("Select SRA Toolkit bin folder", self.sratoolkit_path)
        if folder:
            self.sratoolkit_path = folder
            self.save_memory()

    def select_sra_folder(self):
        folder = self.select_folder_dialog("Select SRA download folder", self.sra_download_folder)
        if folder:
            self.sra_download_folder = folder
            self.save_memory()

    # ---------------- основной поток ----------------

    def start_process(self):
        if not all([self.sratoolkit_path, self.sra_download_folder, self.fastq_output_folder]):
            QMessageBox.warning(self, "Error", "Fill in all paths before starting.")
            return

        samples_raw = self.input_samples.toPlainText().strip()
        if not samples_raw:
            QMessageBox.warning(self, "Error", "Enter pairs SRXID-SampleName")
            return

        rename_mapping = [entry.strip() for entry in samples_raw.split(";") if entry.strip()]
        if not rename_mapping:
            QMessageBox.warning(self, "Error", "No valid entries for processing")
            return

        self.progress_bar.setMaximum(len(rename_mapping))
        self.progress_bar.setValue(0)

        for i, entry in enumerate(rename_mapping, 1):
            try:
                runid, samplename = entry.split("-")
            except ValueError:
                log(f"Invalid format for string: {entry}")
                continue

            self.progress_label.setText(f"Downloading {runid}... [{i}/{len(rename_mapping)}]")
            before_prefetch = set(os.listdir(self.sra_download_folder))
            prefetch_path = os.path.join(self.sratoolkit_path, "prefetch.exe")
            if not self.execute_command_with_error_handling(
                [prefetch_path, runid, "--output-directory", self.sra_download_folder],
                f"Prefetch {runid}"
            ):
                log(f"Prefetch step for {runid} skipped.")
                continue

            after_prefetch = set(os.listdir(self.sra_download_folder))
            new_dirs = list(after_prefetch - before_prefetch)
            if not new_dirs:
                log(f"No new folder found after prefetch for {runid}")
                continue

            downloaded_folder = new_dirs[0]
            downloaded_folder_path = os.path.join(self.sra_download_folder, downloaded_folder)
            sra_files = [f for f in os.listdir(downloaded_folder_path) if f.endswith('.sra')]
            if not sra_files:
                log(f".sra file not found in {downloaded_folder_path}")
                continue

            sra_file_name = sra_files[0]
            src_file_path = os.path.join(downloaded_folder_path, sra_file_name)
            dst_file_path = os.path.join(self.sra_download_folder, f"{samplename}.sra")

            skip_this_sample = False
            while True:
                try:
                    shutil.move(src_file_path, dst_file_path)
                    os.rmdir(downloaded_folder_path)
                    break
                except Exception as e:
                    choice = self.handle_error(f"Error moving {src_file_path}: {e}")
                    if choice == "retry":
                        continue
                    elif choice == "skip":
                        skip_this_sample = True
                        break
                    elif choice == "abort":
                        sys.exit("Processes terminated by the user.")
            if skip_this_sample:
                log("File moving step skipped, moving to next sample.")
                continue

            self.convert_sra_file(dst_file_path, samplename)
            self.progress_bar.setValue(i)

        self.progress_label.setText("Done.")
        QMessageBox.information(self, "Done", "Process completed!")
        self.run_pipeline()

    def convert_existing_sra_files(self):
        sra_files = [f for f in os.listdir(self.sra_download_folder) if f.endswith('.sra')]
        if not sra_files:
            QMessageBox.warning(self, "No files", "No SRA files for conversion.")
            return

        self.progress_bar.setMaximum(len(sra_files))
        self.progress_bar.setValue(0)

        for i, sra_file in enumerate(sra_files, 1):
            sra_path = os.path.join(self.sra_download_folder, sra_file)
            sample_name = sra_file.replace('.sra', '')

            self.progress_label.setText(f"Converting {sample_name}... [{i}/{len(sra_files)}]")
            self.convert_sra_file(sra_path, sample_name)
            self.progress_bar.setValue(i)

        self.progress_label.setText("Done.")
        QMessageBox.information(self, "Done", "Conversion completed!")

        pipeline_script = os.path.join(SCRIPT_DIR, "runpaiplain.py")
        subprocess.run([sys.executable, pipeline_script], check=True)
        self.close()

    # ---------------- QC / Trim / MultiQC ----------------

    def run_qc_and_trim(self, sample_name: str):
        """Оркестрация: FastQC (pre) → [Cutadapt] → FastQC (post) → MultiQC."""
        qc_opts = self.pipeline_settings.get("qc", {})
        folders = self.pipeline_settings.get("folders", {})
        fastq_dir = self.fastq_output_folder or folders.get("fastq_folder", "")
        results_dir = folders.get("results_folder", "")

        if not fastq_dir or not os.path.isdir(fastq_dir):
            log("QC: FASTQ folder not set or inaccessible.")
            return
        if not results_dir:
            log("QC: results_folder not set; using FASTQ folder for reports.")
            results_dir = fastq_dir

        r1 = os.path.join(fastq_dir, f"{sample_name}_1.fastq")
        r2 = os.path.join(fastq_dir, f"{sample_name}_2.fastq")
        se = os.path.join(fastq_dir, f"{sample_name}.fastq")
        is_pe = os.path.exists(r1) and os.path.exists(r2)
        is_se = os.path.exists(se)
        if not (is_pe or is_se):
            log(f"QC: FASTQ for sample {sample_name} not found.")
            return

        qc_root = os.path.join(results_dir, "qc")
        qc_pre = os.path.join(qc_root, "raw")
        qc_post = os.path.join(qc_root, "trimmed")
        os.makedirs(qc_pre, exist_ok=True)
        os.makedirs(qc_post, exist_ok=True)

        # FastQC (pre)
        if qc_opts.get("enable_fastqc", True):
            files_pre = [r1, r2] if is_pe else [se]
            self._fastqc(files_pre, qc_pre)

        # Тримминг (по желанию)
        if qc_opts.get("enable_trimming", False):
            keep_raw = qc_opts.get("keep_raw_fastq", True)
            backup_dir = os.path.join(fastq_dir, "fastq_backup")
            tmp_dir = os.path.join(fastq_dir, "fastq_trimmed")
            os.makedirs(tmp_dir, exist_ok=True)
            if keep_raw:
                os.makedirs(backup_dir, exist_ok=True)

            if is_pe:
                r1_out_tmp = os.path.join(tmp_dir, f"{sample_name}_1.fastq")
                r2_out_tmp = os.path.join(tmp_dir, f"{sample_name}_2.fastq")
                self._cutadapt_pe(
                    r1, r2, r1_out_tmp, r2_out_tmp,
                    a1=qc_opts.get("adapter_r1"),
                    a2=qc_opts.get("adapter_r2"),
                    q=qc_opts.get("trim_q", 20),
                    minlen=qc_opts.get("trim_minlen_pe", 50)
                )
                if keep_raw:
                    shutil.move(r1, os.path.join(backup_dir, os.path.basename(r1)))
                    shutil.move(r2, os.path.join(backup_dir, os.path.basename(r2)))
                else:
                    os.remove(r1); os.remove(r2)
                shutil.move(r1_out_tmp, r1)
                shutil.move(r2_out_tmp, r2)
            else:
                se_out_tmp = os.path.join(tmp_dir, f"{sample_name}.fastq")
                self._cutadapt_se(
                    se, se_out_tmp,
                    a=qc_opts.get("adapter_r1"),
                    q=qc_opts.get("trim_q", 20),
                    minlen=qc_opts.get("trim_minlen_se", 35)
                )
                if keep_raw:
                    shutil.move(se, os.path.join(backup_dir, os.path.basename(se)))
                else:
                    os.remove(se)
                shutil.move(se_out_tmp, se)

            # FastQC (post)
            if qc_opts.get("enable_fastqc", True):
                files_post = [r1, r2] if is_pe else [se]
                self._fastqc(files_post, qc_post)

        # MultiQC (агрегация)
        if qc_opts.get("enable_fastqc", True):
            self._multiqc(qc_root)

                # Очистка артефактов QC при включённой настройке
        if self.pipeline_settings["options"].get("delete_intermediate_files", False):
            report = os.path.join(qc_root, "multiqc_report.html")
            if os.path.exists(report):
                self._cleanup_qc_artifacts(qc_root, sample_name)
            else:
                log("QC cleanup: multiqc_report.html не найден -> пропуск удаления QC-артефактов.")

    def _fastqc(self, files, out_dir, threads=4):
        """Устойчивый запуск FastQC на Windows-сборке; сбор отчётов в out_dir независимо от фактического места вывода."""
        import time

        os.makedirs(out_dir, exist_ok=True)
        PROJECT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
        FASTQC_DIR = os.path.join(PROJECT_DIR, "Tool", "FastQC")

        candidates = [
            os.path.join(FASTQC_DIR, "fastqc.bat"),
            os.path.join(FASTQC_DIR, "run_fastqc.bat"),
            shutil.which("fastqc"),
        ]
        fastqc_path = next((p for p in candidates if p and os.path.exists(p)), None)
        if not fastqc_path:
            raise FileNotFoundError("FastQC not found. Install into C:\\PipeSeq-3\\Tool\\FastQC or add to PATH.")

        t0 = time.time()
        window = 3600  # 1 час

        def _is_fresh(p):
            try:
                return (t0 - os.path.getmtime(p)) <= window
            except Exception:
                return False

        def collect(dst_dir) -> int:
            moved = 0
            search_dirs = {FASTQC_DIR, os.path.dirname(os.path.abspath(files[0])), dst_dir}
            pats = ["*_fastqc.html", "*_fastqc.zip"]

            # прямой поиск
            for src in list(search_dirs):
                for pat in pats:
                    for f in glob.glob(os.path.join(src, pat)):
                        if not _is_fresh(f):
                            continue
                        try:
                            tgt = os.path.join(dst_dir, os.path.basename(f))
                            if os.path.abspath(os.path.dirname(f)) != os.path.abspath(dst_dir):
                                shutil.move(f, tgt)
                            moved += 1
                        except Exception:
                            pass
            # при необходимости — рекурсивный обход по проекту
            if moved == 0:
                for pat in pats:
                    for f in glob.glob(os.path.join(PROJECT_DIR, "**", pat), recursive=True):
                        if not _is_fresh(f):
                            continue
                        try:
                            tgt = os.path.join(dst_dir, os.path.basename(f))
                            if os.path.abspath(os.path.dirname(f)) != os.path.abspath(dst_dir):
                                shutil.move(f, tgt)
                            moved += 1
                        except Exception:
                            pass
            return moved

        # запуск без флагов (данная сборка может игнорировать CLI-опции)
        cmd = (["cmd", "/c", fastqc_path] if fastqc_path.lower().endswith((".bat", ".cmd")) else [fastqc_path]) + files
        self.execute_command_with_error_handling(cmd, f"FastQC > {os.path.basename(out_dir)}", cwd=FASTQC_DIR, ok_returncodes=(0, 1))

        moved = collect(out_dir)
        have_reports = any(glob.glob(os.path.join(out_dir, "*_fastqc.html")) + glob.glob(os.path.join(out_dir, "*_fastqc.zip")))
        if not (moved > 0 or have_reports):
            raise RuntimeError("FastQC produced no reports (searched Tool\\FastQC, FASTQ folder, project tree).")

    def _multiqc(self, qc_root: str, sample_name: str | None = None) -> None:
        """
        Агрегация FastQC-отчётов. Создаёт временный file-list, принудительно
        перезаписывает итоговый отчёт 'multiqc_report.html' и не прерывает конвейер.
        """
        import sys, os, glob, shutil, subprocess, tempfile, time

        os.makedirs(qc_root, exist_ok=True)

        # --- формируем подмножество входов (только текущий образец; иначе свежие за 6 ч) ---
        candidates: list[str] = []
        search_dirs = [os.path.join(qc_root, "raw"), os.path.join(qc_root, "trimmed")]
        search_dirs = [d for d in search_dirs if os.path.isdir(d)]

        def add_glob(pat: str) -> None:
            for d in search_dirs:
                candidates.extend(glob.glob(os.path.join(d, pat)))

        if sample_name:
            add_glob(f"{sample_name}_*_fastqc.zip")
            add_glob(f"{sample_name}_fastqc.zip")
            for d in search_dirs:
                candidates.extend(glob.glob(os.path.join(d, f"{sample_name}_*_fastqc", "fastqc_data.txt")))
                candidates.extend(glob.glob(os.path.join(d, f"{sample_name}_fastqc", "fastqc_data.txt")))
        else:
            t0, window = time.time(), 6 * 3600
            add_glob("*_fastqc.zip")
            candidates = [p for p in candidates if (t0 - os.path.getmtime(p)) <= window]
            for d in search_dirs:
                for p in glob.glob(os.path.join(d, "*_fastqc", "fastqc_data.txt")):
                    try:
                        if (t0 - os.path.getmtime(p)) <= window:
                            candidates.append(p)
                    except Exception:
                        pass

        if not candidates:
            try: log("MultiQC: нет релевантных отчётов для агрегации — шаг пропущен.")
            except Exception: pass
            return

        # --- временный file-list внутри qc_root и гарантированное удаление ---
        fd, filelist = tempfile.mkstemp(prefix="multiqc_filelist_", suffix=".txt", dir=qc_root, text=True)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                for p in sorted(set(map(os.path.abspath, candidates))):
                    fh.write(p + "\n")

            exe = shutil.which("multiqc")
            base = ([exe] if exe else [sys.executable, "-m", "multiqc"])
            cmd = base + [
                "-q",            # тихий режим
                "-f",            # force overwrite
                "-n", "multiqc_report",  # единое имя отчёта и data-папки
                "-o", qc_root,
                "--ignore", "*_fastqc.html",
                "--ignore", "multiqc_*",
                "--file-list", filelist
            ]
            log(f"Executing command: {' '.join(cmd)} (stage: MultiQC)")
            subprocess.run(cmd, cwd=qc_root, check=False)
        except KeyboardInterrupt:
            log("MultiQC: прервано пользователем; стадия пропущена, конвейер продолжается.")
        except Exception as e:
            log(f"MultiQC: ошибка '{e}'; стадия пропущена.")
        finally:
            try:
                os.remove(filelist)  # не оставляем multiqc_filelist_*.txt
                log(f"MultiQC: удалён временный список {filelist}")
            except Exception:
                pass


    def _cleanup_qc_artifacts(self, qc_root: str, sample_name: str | None = None) -> None:
        """
        Удаляет артефакты FastQC/MultiQC. Сохраняет только один итоговый отчёт
        'multiqc_report.html'. FASTQ не затрагиваются.
        """
        import os, glob, shutil

        def rm_file(p):
            try: os.remove(p); log(f"QC cleanup: deleted file {p}")
            except FileNotFoundError: pass
            except Exception as e: log(f"QC cleanup warning (file): {p} -> {e}")

        def rm_dir(p):
            try: shutil.rmtree(p, ignore_errors=True); log(f"QC cleanup: deleted dir {p}")
            except Exception as e: log(f"QC cleanup warning (dir): {p} -> {e}")

        # --- FastQC: html/zip + каталоги образца или все, если sample_name=None
        pat_files = ["*_fastqc.html", "*_fastqc.zip"]
        pat_dirs  = ["*_fastqc"]
        if sample_name:
            pat_files = [
                f"{sample_name}_fastqc.html",  f"{sample_name}_1_fastqc.html",  f"{sample_name}_2_fastqc.html",
                f"{sample_name}_fastqc.zip",   f"{sample_name}_1_fastqc.zip",   f"{sample_name}_2_fastqc.zip"
            ]
            pat_dirs  = [f"{sample_name}_fastqc", f"{sample_name}_1_fastqc", f"{sample_name}_2_fastqc"]

        for sub in ("raw", "trimmed"):
            d = os.path.join(qc_root, sub)
            if not os.path.isdir(d): continue
            for pat in pat_files:
                for p in glob.glob(os.path.join(d, pat)): rm_file(p)
            for pat in pat_dirs:
                for p in glob.glob(os.path.join(d, pat)): rm_dir(p)

        # --- MultiQC: сохраняем только один финальный отчёт (multiqc_report.html) ---
        # Удаляем вспомогательные каталоги и filelist-файлы
        for pat in ("multiqc_data*", "multiqc_plots*", "multiqc_fastqc_*", "multiqc.log", "multiqc_filelist_*.txt"):
            for p in glob.glob(os.path.join(qc_root, pat)):
                rm_dir(p) if os.path.isdir(p) else rm_file(p)

        # Удаляем любые дубликаты отчётов, оставляя фиксированный 'multiqc_report.html'
        for p in glob.glob(os.path.join(qc_root, "multiqc_report*.html")):
            if os.path.basename(p) != "multiqc_report.html":
                rm_file(p)

        # Временные каталоги тримминга (если пустые)
        fastq_dir = self.fastq_output_folder or self.pipeline_settings.get("folders", {}).get("fastq_folder", "")
        for name in ("fastq_trimmed", "fastq_backup"):
            p = os.path.join(fastq_dir, name)
            try:
                if os.path.isdir(p) and not os.listdir(p):
                    rm_dir(p)
            except Exception:
                pass

    def _cutadapt_pe(self, r1_in, r2_in, r1_out, r2_out, a1, a2, q=20, minlen=50):
        if not a1: a1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
        if not a2: a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        cmd = [
            "cutadapt", "-j", "4",
            "-a", a1, "-A", a2,
            "-q", f"{q},{q}", "--minimum-length", str(minlen),
            "-o", r1_out, "-p", r2_out, r1_in, r2_in
        ]
        self.execute_command_with_error_handling(cmd, "Cutadapt (PE)")

    def _cutadapt_se(self, r_in, r_out, a, q=20, minlen=35):
        if not a: a = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
        cmd = [
            "cutadapt", "-j", "4",
            "-a", a,
            "-q", str(q), "--minimum-length", str(minlen),
            "-o", r_out, r_in
        ]
        self.execute_command_with_error_handling(cmd, "Cutadapt (SE)")

    def convert_sra_file(self, sra_file_path, sample_name):
        fasterq_path = os.path.join(self.sratoolkit_path, "fasterq-dump.exe")
        log(f"Starting conversion {sample_name} → FASTQ")
        command = [fasterq_path, sra_file_path, "--split-files", "--outdir", self.fastq_output_folder]
        if not self.execute_command_with_error_handling(command, f"Conversion {sample_name}"):
            log(f"Conversion skipped for {sample_name}")
            return

        if self.delete_sra_after_conversion:
            while True:
                try:
                    os.remove(sra_file_path)
                    log(f"SRA file deleted: {sra_file_path}")
                    break
                except Exception as e:
                    choice = self.handle_error(f"Error deleting {sra_file_path}: {e}")
                    if choice == "retry":
                        continue
                    elif choice == "skip":
                        break
                    elif choice == "abort":
                        sys.exit("Processes terminated by the user.")

        # ВАЖНО: QC запускается всегда
        self.run_qc_and_trim(sample_name)

        align_script = os.path.join(SCRIPT_DIR, "align_hisat2.py")
        log(f"Starting alignment for {sample_name}")
        if not self.execute_command_with_error_handling([sys.executable, align_script, sample_name], f"Alignment {sample_name}"):
            log(f"Alignment skipped for {sample_name}")
            return

        process_script = os.path.join(SCRIPT_DIR, "process_sam_to_bam.py")
        log(f"Starting SAM → BAM conversion for {sample_name}")
        if not self.execute_command_with_error_handling([sys.executable, process_script, sample_name], f"SAM → BAM conversion {sample_name}"):
            log(f"SAM → BAM conversion skipped for {sample_name}")

        if self.pipeline_settings["options"].get("delete_intermediate_files", False):
            fastq_files = [f for f in os.listdir(self.fastq_output_folder) if f.startswith(sample_name) and f.endswith('.fastq')]
            for fastq in fastq_files:
                fastq_path = os.path.join(self.fastq_output_folder, fastq)
                try:
                    os.remove(fastq_path)
                    log(f"Intermediate FASTQ file deleted: {fastq_path}")
                except Exception as e:
                    log(f"Error deleting file {fastq_path}: {e}")

    def run_pipeline(self):
        pipeline_script = os.path.join(SCRIPT_DIR, "run_pipeline_remaining.py")
        while True:
            try:
                subprocess.run([sys.executable, pipeline_script], check=True)
                self.close()
                break
            except subprocess.CalledProcessError as e:
                choice = self.handle_error(f"Pipeline startup error: {e}")
                if choice == "retry":
                    continue
                elif choice == "skip":
                    log("Pipeline startup step skipped.")
                    self.close()
                    break
                elif choice == "abort":
                    sys.exit("Processes terminated by the user.")

    def handle_error(self, error_details):
        msg_box = QMessageBox(self)
        msg_box.setIcon(QMessageBox.Icon.Critical)
        msg_box.setWindowTitle("Step Error")
        msg_box.setText("An error occurred:")
        msg_box.setInformativeText(str(error_details))
        retry_button = msg_box.addButton("Retry", QMessageBox.ButtonRole.AcceptRole)
        skip_button = msg_box.addButton("Skip step", QMessageBox.ButtonRole.DestructiveRole)
        abort_button = msg_box.addButton("Abort processes", QMessageBox.ButtonRole.RejectRole)
        msg_box.setDefaultButton(retry_button)
        msg_box.exec()
        clicked = msg_box.clickedButton()
        if clicked == retry_button:
            return "retry"
        elif clicked == skip_button:
            return "skip"
        elif clicked == abort_button:
            return "abort"
        return "abort"

    def execute_command_with_error_handling(self, command, stage, cwd=None, ok_returncodes=(0,)):
        while True:
            try:
                log(f"Executing command: {' '.join(command)} (stage: {stage})")
                cp = subprocess.run(command, cwd=cwd)
                if cp.returncode in ok_returncodes:
                    return True
                raise subprocess.CalledProcessError(cp.returncode, command)
            except subprocess.CalledProcessError as e:
                choice = self.handle_error(f"Command error: {' '.join(command)}\nDescription: {e}")
                if choice == "retry":
                    continue
                elif choice == "skip":
                    return False
                elif choice == "abort":
                    sys.exit("Processes terminated by the user.")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SRAConverterApp()
    window.show()
    sys.exit(app.exec())
