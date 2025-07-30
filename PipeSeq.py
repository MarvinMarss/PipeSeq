import sys
import os
import subprocess
import time
import json
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

def log(message):
    with open(LOG_FILE, 'a', encoding='utf-8') as f:
        f.write(message + '\n')
    print(message)

def run_command(command, cwd=None):
    """Execute command directly, without WSL"""
    log(f"Executing command: {' '.join(command)}")
    try:
        subprocess.run(command, check=True, cwd=cwd)
    except subprocess.CalledProcessError as e:
        log(f"Error executing command: {e}")
        raise

def normalize_win_path(path):
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
            "visualization": {
                "show_p_values": True
            }
        }

    def save_pipeline_settings(self):
        with open(PIPELINE_SETTINGS_FILE, "w", encoding='utf-8') as f:
            json.dump(self.pipeline_settings, f, indent=4)

    def select_folder_dialog(self, title="Select folder", start_path=None):
        folder = QFileDialog.getExistingDirectory(
            self,
            title,
            start_path or "C:/",
            QFileDialog.Option.ShowDirsOnly
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

        self.btn_sratoolkit = QPushButton("SRA Toolkit")
        self.btn_sratoolkit.clicked.connect(self.select_sratoolkit)
        self.btn_sra = QPushButton("SRA Files")
        self.btn_sra.clicked.connect(self.select_sra_folder)
        self.btn_fastq = QPushButton("FASTQ")
        self.btn_fastq.clicked.connect(self.select_combined_fastq_folder)
        self.btn_bam = QPushButton("BAM/Output")
        self.btn_bam.clicked.connect(lambda: self.select_pipeline_folder("bam_folder", "BAM/Output"))
        self.btn_gtf = QPushButton("GTF")
        self.btn_gtf.clicked.connect(lambda: self.select_pipeline_folder("gtf_folder", "GTF"))
        self.btn_results = QPushButton("Results")
        self.btn_results.clicked.connect(lambda: self.select_pipeline_folder("results_folder", "Results"))
        self.btn_genome = QPushButton("GENOME")
        self.btn_genome.clicked.connect(lambda: self.select_pipeline_folder("genome_folder", "GENOME"))
        self.btn_genome_index = QPushButton("HISAT2 index")
        self.btn_genome_index.clicked.connect(lambda: self.select_pipeline_folder("genome_index", "HISAT2 index"))

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
        grid_pipeline.setHorizontalSpacing(10)
        grid_pipeline.setVerticalSpacing(6)
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
        self.pipeline_sensitivity_input = QLineEdit()
        self.pipeline_sensitivity_input.setPlaceholderText("For example: 0.1")
        current_s = self.pipeline_settings["options"].get("stringtie_sensitivity", "")
        if current_s:
            self.pipeline_sensitivity_input.setText(str(current_s))
        self.pipeline_sensitivity_input.editingFinished.connect(self.update_pipeline_sensitivity)
        hbox_sensitivity.addWidget(lbl_sensitivity)
        hbox_sensitivity.addWidget(self.pipeline_sensitivity_input)
        hbox_sensitivity.setSpacing(6)

        btn_gene_mapping = QPushButton("Edit gene list")
        btn_gene_mapping.clicked.connect(self.edit_pipeline_gene_mapping)

        vbox_pipeline = QVBoxLayout()
        vbox_pipeline.addLayout(grid_pipeline)
        vbox_pipeline.addLayout(hbox_sensitivity)
        vbox_pipeline.addWidget(btn_gene_mapping, alignment=Qt.AlignmentFlag.AlignCenter)
        group_pipeline.setLayout(vbox_pipeline)
        main_layout.addWidget(group_pipeline)

        group_sra = QGroupBox("SRA Conversion")
        vbox_sra = QVBoxLayout()
        vbox_sra.setSpacing(6)
        vbox_sra.setContentsMargins(10, 10, 10, 10)

        lbl_samples = QLabel("Enter pairs SRXID-SampleName\n(format: SRX8380271-HighLight1; SRX5120532-HighLightControl1):")
        vbox_sra.addWidget(lbl_samples)
        self.input_samples = QTextEdit()
        self.input_samples.setFixedHeight(60)
        vbox_sra.addWidget(self.input_samples)
        
        hbox_sra_options = QHBoxLayout()
        chk_delete_sra = QCheckBox("Delete SRA after conversion")
        chk_delete_sra.setChecked(self.delete_sra_after_conversion)
        chk_delete_sra.stateChanged.connect(lambda state: self.toggle_delete_sra(state))
        hbox_sra_options.addWidget(chk_delete_sra)
        
        self.progress_label = QLabel("Awaiting actions...")
        hbox_sra_options.addWidget(self.progress_label, stretch=1)
        vbox_sra.addLayout(hbox_sra_options)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        vbox_sra.addWidget(self.progress_bar)
        
        hbox_sra_buttons = QHBoxLayout()
        btn_download = QPushButton("Download and convert")
        btn_download.clicked.connect(self.start_process)
        btn_convert_local = QPushButton("Convert local SRA")
        btn_convert_local.clicked.connect(self.convert_existing_sra_files)
        hbox_sra_buttons.addWidget(btn_download)
        hbox_sra_buttons.addWidget(btn_convert_local)
        vbox_sra.addLayout(hbox_sra_buttons)
        
        group_sra.setLayout(vbox_sra)
        main_layout.addWidget(group_sra)

        btn_run_pipeline = QPushButton("Run PipeSeq-2")
        btn_run_pipeline.clicked.connect(self.run_pipeline)
        main_layout.addWidget(btn_run_pipeline, alignment=Qt.AlignmentFlag.AlignCenter)

        self.setLayout(main_layout)
        self.setWindowTitle("PipeSeq-1")
        self.resize(600, 650)

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
            if not self.execute_command_with_error_handling([prefetch_path, runid, "--output-directory", self.sra_download_folder],
                                                              f"Prefetch {runid}"):
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

        align_script = os.path.join(SCRIPT_DIR, "align_hisat2.py")
        log(f"Starting alignment for {sample_name}")
        if not self.execute_command_with_error_handling([sys.executable, align_script, sample_name],
                                                        f"Alignment {sample_name}"):
            log(f"Alignment skipped for {sample_name}")
            return

        process_script = os.path.join(SCRIPT_DIR, "process_sam_to_bam.py")
        log(f"Starting SAM → BAM conversion for {sample_name}")
        if not self.execute_command_with_error_handling([sys.executable, process_script, sample_name],
                                                        f"SAM → BAM conversion {sample_name}"):
            log(f"SAM → BAM conversion skipped for {sample_name}")

        if self.pipeline_settings["options"].get("delete_intermediate_files", False):
            fastq_files = [f for f in os.listdir(self.fastq_output_folder)
                           if f.startswith(sample_name) and f.endswith('.fastq')]
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

    def execute_command_with_error_handling(self, command, stage_desc, cwd=None):
        while True:
            try:
                log(f"Executing command: {' '.join(command)} (stage: {stage_desc})")
                subprocess.run(command, check=True, cwd=cwd)
                return True
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
