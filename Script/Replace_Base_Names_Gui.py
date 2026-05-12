import sys
import os
import pandas as pd
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QPushButton, QLabel,
    QFileDialog, QLineEdit, QFormLayout, QMessageBox
)
from PyQt6.QtCore import Qt

class BaseNameReplacer(QWidget):
    def __init__(self):
        super().__init__()
        self.base_names = {}
        self.name_inputs = {}
        self.file_path = None
        self.setWindowTitle('Замена Base Name')
        self.setMinimumWidth(500)
        self.form_layout = QFormLayout()
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        self.label = QLabel('Выберите файл для замены Base Name')
        layout.addWidget(self.label)

        self.btn_open = QPushButton('Открыть файл')
        self.btn_open.clicked.connect(self.openFile)
        layout.addWidget(self.btn_open)

        layout.addLayout(self.form_layout)

        self.btn_replace = QPushButton('Заменить и перезаписать файл')
        self.btn_replace.clicked.connect(self.replaceNames)
        self.btn_replace.setEnabled(False)
        layout.addWidget(self.btn_replace)

        self.setLayout(layout)

    def openFile(self):
        file_name, _ = QFileDialog.getOpenFileName(
            self,
            'Выбери файл',
            '',
            'Text Files (*.txt);;All Files (*)'
        )
        if file_name:
            self.file_path = file_name
            self.label.setText(f'Выбран файл: {file_name}')
            self.loadBaseNames()

    def loadBaseNames(self):
        try:
            df = pd.read_csv(self.file_path, sep='\t')
        except Exception as e:
            QMessageBox.critical(self, "Ошибка", f"Не удалось прочитать файл:\n{e}")
            return

        if 'Base Name' not in df.columns:
            QMessageBox.critical(self, "Ошибка", "В файле нет колонки 'Base Name'")
            return


        while self.form_layout.count():
            child = self.form_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        self.name_inputs.clear()

        unique_names = sorted(df['Base Name'].dropna().unique())
        for name in unique_names:
            input_field = QLineEdit()
            input_field.setPlaceholderText("Новое имя (оставь пустым, чтобы не менять)")
            self.name_inputs[name] = input_field
            self.form_layout.addRow(f"{name} ➜", input_field)

        self.btn_replace.setEnabled(True)

    def replaceNames(self):
        if not self.file_path:
            QMessageBox.warning(self, "Ошибка", "Сначала выбери файл.")
            return

        try:
            df = pd.read_csv(self.file_path, sep='\t')
        except Exception as e:
            QMessageBox.critical(self, "Ошибка", f"Не удалось прочитать файл:\n{e}")
            return

        changes_made = False
        for old_name, input_field in self.name_inputs.items():
            new_name = input_field.text().strip()
            if new_name and new_name != old_name:
                df['Base Name'] = df['Base Name'].replace(old_name, new_name)
                changes_made = True

        if not changes_made:
            QMessageBox.information(self, "Без изменений", "Никакие имена не были изменены.")
            return

        try:
            df.to_csv(self.file_path, sep='\t', index=False)
            QMessageBox.information(self, "Готово", f"Файл успешно перезаписан:\n{self.file_path}")
            self.label.setText(f'Файл перезаписан: {self.file_path}')
        except Exception as e:
            QMessageBox.critical(self, "Ошибка", f"Не удалось сохранить файл:\n{e}")

def launch_gui():
    try:
        app = QApplication.instance()
        if not app:
            app = QApplication(sys.argv)

        window = BaseNameReplacer()
        window.show()


        sys.exit(app.exec())
    except Exception as e:
        import traceback
        with open("Replace_Base_Names_error.log", "w", encoding="utf-8") as f:
            f.write("Ошибка запуска GUI:\n")
            f.write(traceback.format_exc())

if __name__ == '__main__':
    launch_gui()
