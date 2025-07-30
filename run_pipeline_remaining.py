import sys
import os
import json
import subprocess
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QLabel, QPushButton, QTextEdit, QMessageBox,
    QProgressBar
)
from PyQt6.QtGui import QFont
from PyQt6.QtCore import Qt

script_dir = os.path.dirname(os.path.abspath(__file__))
SETTINGS_FILE = os.path.join(script_dir, "settings.json")
PIPELINE_LOG = os.path.join(script_dir, "run_pipeline_log.txt")


STEPS_STRINGTIE = [
    ("stringtie_expression.py", "Подсчёт экспрессии StringTie"),
    ("extract_fpkm.py", "Извлечение FPKM + log2"),
    ("GTF_results_pvalues.py", "Расчёт p-values"),
    ("pvalues_log2.py", "Объединение p-values и log2")
]

STEPS_DESEQ2 = [
    ("deseq2_analysis.py", "DESeq2 анализ"),
    ("extract_Deseq2.py", "Извлечение результатов DESeq2")
]

class PipelineApp(QWidget):
    def __init__(self):
        super().__init__()
        self.settings = {}
        self.load_settings()
        self.init_ui()

    def load_settings(self):
        if os.path.exists(SETTINGS_FILE):
            with open(SETTINGS_FILE, "r") as f:
                self.settings = json.load(f)
        else:
            self.settings = {
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
                    "use_deseq2": False
                },
                "gene_mapping": {},
                "visualization": {
                    "show_p_values": True
                }
            }
            self.save_settings()

    def save_settings(self):
        with open(SETTINGS_FILE, "w") as f:
            json.dump(self.settings, f, indent=4)

    def init_ui(self):
        self.setWindowTitle("PipeSeq - Запуск пайплайна")
        self.resize(600, 500)

        self.setStyleSheet("""
            QWidget { background-color: #2b2b40; color: #f5f5fa; font-family: Arial; font-size: 14px; }
            QPushButton { background-color: #3b3b5c; border: 1px solid #444; padding: 6px; border-radius: 6px; }
            QPushButton:hover { background-color: #50507a; }
            QTextEdit { background-color: #3b3b5c; color: #f5f5fa; border-radius: 5px; padding: 5px; }
            QLabel { color: #f5f5fa; }
        """)

        layout = QVBoxLayout()

        title = QLabel("PipeSeq - Запуск пайплайна")
        title.setFont(QFont("Arial", 20, QFont.Weight.Bold))
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)



        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        layout.addWidget(self.log_output)

        self.run_pipeline_btn = QPushButton("Запустить пайплайн")
        self.run_pipeline_btn.clicked.connect(self.run_pipeline)
        layout.addWidget(self.run_pipeline_btn)

        self.run_visualization_btn = QPushButton("Визуализация (Heatmap)")
        self.run_visualization_btn.clicked.connect(self.run_visualization)
        layout.addWidget(self.run_visualization_btn)

        self.run_replace_btn = QPushButton("Заменить Base Names (GUI)")
        self.run_replace_btn.clicked.connect(self.run_replace_base_names)
        layout.addWidget(self.run_replace_btn)

        self.cleanup_btn = QPushButton("Очистить всё (кроме генома)")
        self.cleanup_btn.clicked.connect(self.cleanup_all_data)
        layout.addWidget(self.cleanup_btn)

        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        layout.addWidget(self.progress_bar)

        self.setLayout(layout)

    def log(self, message):
        self.log_output.append(message)
        with open(PIPELINE_LOG, "a", encoding="utf-8") as log_file:
            log_file.write(message + "\n")

    def show_error_dialog(self, title, message):
        """
        Выводит окно с описанием ошибки и тремя вариантами:
        «Запустить снова», «Пропустить шаг» и «Завершить процессы».
        Возвращает: "retry", "skip" или "cancel"
        """
        msg_box = QMessageBox(self)
        msg_box.setIcon(QMessageBox.Icon.Critical)
        msg_box.setWindowTitle(title)
        msg_box.setText(message)

        retry_btn = msg_box.addButton("Запустить снова", QMessageBox.ButtonRole.AcceptRole)
        skip_btn = msg_box.addButton("Пропустить шаг", QMessageBox.ButtonRole.DestructiveRole)
        cancel_btn = msg_box.addButton("Завершить процессы", QMessageBox.ButtonRole.RejectRole)

        msg_box.exec()

        if msg_box.clickedButton() == retry_btn:
            return "retry"
        elif msg_box.clickedButton() == skip_btn:
            return "skip"
        else:
            return "cancel"

    def run_pipeline(self):
        self.progress_bar.setValue(0)
        self.log_output.clear()
        self.log("Запуск пайплайна...\n")


        if self.settings["options"].get("fix_genome"):
            self.log("Исправление геномного файла...")
            while True:
                try:
                    fix_script = os.path.join(script_dir, "fix.gtf.py")
                    subprocess.run(["python", fix_script], check=True)
                    self.log("Геном исправлен\n")
                    break
                except subprocess.CalledProcessError as e:
                    self.log(f"Ошибка исправления генома: {e}")
                    choice = self.show_error_dialog("Ошибка исправления генома", f"Ошибка исправления генома: {e}\n\nХотите повторить попытку, пропустить этап или завершить процессы?")
                    if choice == "retry":
                        continue
                    elif choice == "skip":
                        self.log("Этап исправления генома пропущен.\n")
                        break
                    else:
                        self.log("Пайплайн прерван пользователем.")
                        QMessageBox.information(self, "Прервано", "Пайплайн был завершён.")
                        return


        steps = []
        if self.settings["options"].get("use_stringtie", True):
            steps += STEPS_STRINGTIE.copy()
        if self.settings["options"].get("use_deseq2", False):
            for step in STEPS_DESEQ2:
                if step not in steps:
                    steps.append(step)
        if not steps:
            self.log("Не выбран ни один метод анализа.")
            QMessageBox.warning(self, "Ошибка", "В настройках не выбран ни один метод анализа (StringTie или DESeq2).")
            return
        if not any("GTF_results_pvalues.py" in s for s, _ in steps):
            steps.insert(-1, ("GTF_results_pvalues.py", "Расчёт p-values"))

        self.progress_bar.setMaximum(len(steps))
        self.progress_bar.setValue(0)


        for script, description in steps:
            full_path = os.path.join(script_dir, script)
            while True:
                if not os.path.exists(full_path):
                    self.log(f"Скрипт {script} не найден!")
                    choice = self.show_error_dialog("Ошибка", f"Скрипт {script} для шага '{description}' не найден.\n\nХотите запустить этот шаг снова, пропустить его или завершить процессы?")
                    if choice == "retry":
                        continue
                    elif choice == "skip":
                        break
                    else:
                        self.log("Пайплайн прерван пользователем.")
                        QMessageBox.information(self, "Прервано", "Пайплайн был завершён.")
                        return
                self.log(f"Запуск {description} [{script}]...")
                try:
                    subprocess.run(["python", full_path], check=True)
                    self.log(f"{description} завершено.\n")
                    self.progress_bar.setValue(self.progress_bar.value() + 1)
                    break
                except subprocess.CalledProcessError as e:
                    self.log(f"Ошибка в {description}: {e}")
                    choice = self.show_error_dialog("Ошибка", f"Ошибка при выполнении {description}: {e}\n\nХотите повторить шаг, пропустить его или завершить процессы?")
                    if choice == "retry":
                        continue
                    elif choice == "skip":
                        break
                    else:
                        self.log("Пайплайн прерван пользователем.")
                        QMessageBox.information(self, "Прервано", "Пайплайн был завершён.")
                        return


        if self.settings["options"].get("delete_intermediate_files"):
            self.clean_intermediate_files()

        self.log("Пайплайн успешно завершён!\n")
        QMessageBox.information(self, "Готово", "Анализ завершён успешно.")

    def clean_intermediate_files(self):
        self.log("Удаление промежуточных файлов...")
        bam_folder = self.settings["folders"].get("bam_folder", "")
        if bam_folder and os.path.exists(bam_folder):
            for file in os.listdir(bam_folder):
                if file.endswith((".sam", ".bam")) and not file.endswith(".gtf"):
                    path = os.path.join(bam_folder, file)
                    try:
                        os.remove(path)
                        self.log(f"Удалён файл: {path}")
                    except Exception as e:
                        self.log(f"Не удалось удалить {path}: {e}")

    def run_visualization(self):
        script_name = "temp_card_p.py" if self.settings["visualization"].get("show_p_values", True) else "temp_card.py"
        script_path = os.path.join(script_dir, script_name)
        if not os.path.exists(script_path):
            self.log(f"Визуализация: {script_name} не найден!")
            QMessageBox.critical(self, "Ошибка", f"{script_name} отсутствует.")
            return

        self.log(f"Запуск {script_name}...")
        try:
            subprocess.run(["python", script_path], check=True)
            self.log(f"Визуализация завершена.")
        except subprocess.CalledProcessError as e:
            self.log(f"Ошибка в {script_name}: {e}")

    def run_replace_base_names(self):
        script = os.path.join(script_dir, "Replace_Base_Names_Gui.py")
        if not os.path.exists(script):
            self.log(f"Скрипт Replace_Base_Names_Gui.py не найден!")
            QMessageBox.critical(self, "Ошибка", f"Replace_Base_Names_Gui.py отсутствует.")
            return
        self.log("Запуск Replace_Base_Names_Gui.py")
        try:
            subprocess.run(["python", script], check=True)
            self.log(f"Заменены имена Base Name.")
        except subprocess.CalledProcessError as e:
            self.log(f"Ошибка в {script}: {e}")

    def cleanup_all_data(self):
        script_dir_local = os.path.dirname(os.path.abspath(__file__))
        for filename in os.listdir(script_dir_local):
            if filename.endswith("_log.txt"):
                path = os.path.join(script_dir_local, filename)
                try:
                    os.remove(path)
                    self.log(f"Удалён лог: {filename}")
                except Exception as e:
                    self.log(f"Не удалось удалить {filename}: {e}")

        folders_to_clean = ["fastq_folder", "bam_folder", "gtf_folder", "results_folder"]
        for key in folders_to_clean:
            folder = self.settings["folders"].get(key, "")
            if os.path.isdir(folder):
                for filename in os.listdir(folder):
                    path = os.path.join(folder, filename)
                    if key == "bam_folder" and filename.endswith((".gtf", ".fa", ".fasta", ".ht2")):
                        continue  
                    try:
                        if os.path.isfile(path):
                            os.remove(path)
                        elif os.path.isdir(path):
                            import shutil
                            shutil.rmtree(path)
                    except Exception as e:
                        self.log(f"Не удалось удалить {path}: {e}")
                self.log(f"Очищена папка: {key}")
            else:
                self.log(f"Папка {key} не найдена или не указана.")

if __name__ == "__main__":
    from PyQt6.QtCore import QTimer
    app = QApplication(sys.argv)
    window = PipelineApp()
    window.show()

    QTimer.singleShot(0, window.run_pipeline)
    sys.exit(app.exec())
