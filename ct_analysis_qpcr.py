import sys
import os
import re
import numpy as np
import pandas as pd
import json
from scipy.stats import shapiro, ttest_ind, mannwhitneyu
from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QTableWidget, QTableWidgetItem, QMessageBox, QLabel
)
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QInputDialog

class CtAnalysisApp(QWidget):

    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("RT-qPCR Ct Analysis")
        self.resize(800, 600)
        main_layout = QVBoxLayout()

        instr = QLabel(
            "Instructions:\n"
            "- Double-click on the title to rename.\n"
            "- Prefix Control- for control.\n"
            "- Suffixes _1, _2 … for tech/bio repeats.\n"
            "- Ctrl+C/Ctrl+V — copy/paste from Excel."
        )
        instr.setWordWrap(True)
        main_layout.addWidget(instr, alignment=Qt.AlignmentFlag.AlignLeft)

        self.table = QTableWidget(2, 2)
        self.table.setHorizontalHeaderLabels(["Gene_1", "Gene_2"])
        self.table.setVerticalHeaderLabels  (["Sample_1", "Sample_2"])

        self.table.horizontalHeader().setSectionsClickable(True)
        self.table.verticalHeader  ().setSectionsClickable(True)
        self.table.horizontalHeader().sectionDoubleClicked.connect(self.edit_horizontal_header)
        self.table.verticalHeader  ().sectionDoubleClicked.connect(self.edit_vertical_header)
        self.table.horizontalHeader().setSectionsClickable(True)
        self.table.verticalHeader  ().setSectionsClickable(True)
        self.table.horizontalHeader().sectionDoubleClicked.connect(self.edit_horizontal_header)
        self.table.verticalHeader  ().sectionDoubleClicked.connect(self.edit_vertical_header)

        main_layout.addWidget(self.table)


        self.add_row_btn = QPushButton("Add sample")
        self.add_row_btn.clicked.connect(self.add_sample)
        self.add_col_btn = QPushButton("Add gene")
        self.add_col_btn.clicked.connect(self.add_gene)
        self.compute_btn = QPushButton("Calculate ΔCt and statistics")
        self.compute_btn.clicked.connect(self.compute_statistics)
        self.del_row_btn = QPushButton("Remove sample")
        self.del_row_btn.clicked.connect(self.remove_sample)
        self.del_col_btn = QPushButton("Remove gene")
        self.del_col_btn.clicked.connect(self.remove_gene)

        btn_layout = QHBoxLayout()
        btn_layout.addWidget(self.add_row_btn)
        btn_layout.addWidget(self.del_row_btn)    
        btn_layout.addWidget(self.add_col_btn)
        btn_layout.addWidget(self.del_col_btn)    
        btn_layout.addWidget(self.compute_btn)
        main_layout.addLayout(btn_layout)

        self.setLayout(main_layout)


    def keyPressEvent(self, event):
        if event.modifiers() == Qt.KeyboardModifier.ControlModifier:
            if event.key() == Qt.Key.Key_C:
                self.copy_selection()
            elif event.key() == Qt.Key.Key_V:
                self.paste_selection()
        super().keyPressEvent(event)


    def add_sample(self):
        row = self.table.rowCount()
        self.table.insertRow(row)
        self.table.setVerticalHeaderItem(row, QTableWidgetItem(f"Sample_{row+1}"))

    def add_gene(self):
        col = self.table.columnCount()
        self.table.insertColumn(col)
        self.table.setHorizontalHeaderItem(col, QTableWidgetItem(f"Gene_{col+1}"))

    def remove_sample(self):
        row = self.table.currentRow()
        if row < 0 or row >= self.table.rowCount():
            row = self.table.rowCount() - 1
        if row >= 0:
            self.table.removeRow(row)
        else:
            QMessageBox.warning(self, "Removing a sample", "There are no rows available to delete.")

    def remove_gene(self):

        col = self.table.currentColumn()
        if col < 0 or col >= self.table.columnCount():
            col = self.table.columnCount() - 1
        if col >= 0:
            self.table.removeColumn(col)
        else:
            QMessageBox.warning(self, "Removing a gene", "There are no columns available to delete.")


    def parse_headers(self):

        gene_headers = [self.table.horizontalHeaderItem(j).text() 
                        for j in range(self.table.columnCount())]
        genes = [re.sub(r'_(\d+)$', '', h) for h in gene_headers]


        samples = []
        for i in range(self.table.rowCount()):
            s = self.table.verticalHeaderItem(i).text().strip()

            m = re.match(r'^(\d+)(Control-)?(.+?)_(\d+)$', s)
            if not m:
                QMessageBox.critical(self, "Error", 
                    f"Incorrect sample format: «{s}»\n"
                    "Should be: 1Control-Condition, time_1 or 1Condition, time_1")
                return None, None
            exp_num, ctrl_pref, cond, bio_rep = m.groups()
            samples.append({
                'label':    s,
                'exp_num':  int(exp_num),
                'is_ctrl':  bool(ctrl_pref),
                'cond':     cond.strip(),   
                'bio_rep':  int(bio_rep),
            })
        return genes, samples
    
    def edit_horizontal_header(self, index: int):
        
        old = [self.table.horizontalHeaderItem(j).text() for j in range(self.table.columnCount())]
        text, ok = QInputDialog.getText(
            self, 
            "Rename genes",
            "Insert via Tab N names starting from this column:",
            text=old[index]
        )
        if not (ok and text):
            return
        parts = text.split('\t')
        for k, name in enumerate(parts):
            col = index + k
            if col < self.table.columnCount():
                self.table.setHorizontalHeaderItem(col, QTableWidgetItem(name))

    def edit_vertical_header(self, index: int):
        
        old = [self.table.verticalHeaderItem(i).text() for i in range(self.table.rowCount())]
        text, ok = QInputDialog.getText(
            self,
            "Rename samples",
            "Insert via Tab or Enter M names starting from this line:",
            text=old[index]
        )
        if not (ok and text):
            return
       
        parts = re.split(r'[\t\n]+', text)
        for k, name in enumerate(parts):
            row = index + k
            if row < self.table.rowCount():
                self.table.setVerticalHeaderItem(row, QTableWidgetItem(name))

    def copy_selection(self):
        selection = self.table.selectedRanges()[0]
        data = []
        for row in range(selection.topRow(), selection.bottomRow() + 1):
            row_data = []
            for col in range(selection.leftColumn(), selection.rightColumn() + 1):
                item = self.table.item(row, col)
                row_data.append(item.text() if item else "")
            data.append('\t'.join(row_data))
        QApplication.clipboard().setText('\n'.join(data))

    def paste_selection(self):
       
        text = QApplication.clipboard().text()

        rows = [r.split('\t') for r in text.splitlines()]

        sel = self.table.selectedRanges()
        if sel:

            r0, c0 = sel[0].topRow(), sel[0].leftColumn()
            rmax, cmax = sel[0].bottomRow(), sel[0].rightColumn()
        else:

            r0, c0 = self.table.currentRow(), self.table.currentColumn()
            rmax, cmax = self.table.rowCount()-1, self.table.columnCount()-1

        for i, row in enumerate(rows):
            for j, val in enumerate(row):
                r = r0 + i
                c = c0 + j
                if r > rmax or c > cmax:
                    continue

                self.table.setItem(r, c, QTableWidgetItem(val))

    def compute_statistics(self):

        with open("settings.json", "r", encoding="utf-8") as f:
            settings = json.load(f)
        res_dir = settings["folders"]["results_folder"]
        os.makedirs(res_dir, exist_ok=True)
        out_path = os.path.join(res_dir, "ct_analysis_results.txt")


        genes, samples = self.parse_headers()
        if genes is None:
            return

        raw_ct  = {}
        mean_ct = {}
        for i, samp in enumerate(samples):
            lab = samp['label']
            raw_ct[lab]  = {}
            mean_ct[lab] = {}
            for g in genes:
                vals = []
                for j in range(self.table.columnCount()):
                    hdr = self.table.horizontalHeaderItem(j).text()
                    if re.sub(r'_(\d+)$', '', hdr) == g:
                        it = self.table.item(i, j)
                        if it and it.text().strip():
                            try:
                                vals.append(float(it.text().replace(',', '.')))
                            except ValueError:
                                pass
                raw_ct[lab][g]  = vals
                mean_ct[lab][g] = np.nanmean(vals) if vals else np.nan


        ref_genes = [g for g in genes if g.startswith("Control-")]
        if not ref_genes:
            QMessageBox.critical(self, "Error", "Reference gene with prefix Control- was not found")
            return
        ref = ref_genes[0]

        results = []

        uniq = sorted({(s['exp_num'], s['cond']) 
                       for s in samples if not s['is_ctrl']})
        for exp_num, cond in uniq:

            ctrls = [s['label'] for s in samples 
                     if s['exp_num']==exp_num and s['is_ctrl']]

            exps  = [s['label'] for s in samples 
                     if s['exp_num']==exp_num and not s['is_ctrl'] and s['cond']==cond]
            if not ctrls or not exps:
                continue

            for gene in genes:
                if gene == ref:
                    continue  

                d_ctrl = []
                for lab in ctrls:

                    d_ctrl += [x - y for x, y in zip(raw_ct[lab][gene], raw_ct[lab][ref])]
                d_exp = []
                for lab in exps:
                    d_exp += [x - y for x, y in zip(raw_ct[lab][gene], raw_ct[lab][ref])]

                mean_ctrl = np.nanmean(d_ctrl)
                mean_exp  = np.nanmean(d_exp)
                log2fc    = mean_exp - mean_ctrl


                if len(d_ctrl)>=3 and len(d_exp)>=3:
                    p_nc = shapiro(d_ctrl)[1]
                    p_ne = shapiro(d_exp)[1]
                    if p_nc>0.05 and p_ne>0.05:
                        pval = ttest_ind(d_exp, d_ctrl, equal_var=False)[1]
                    else:
                        pval = mannwhitneyu(d_exp, d_ctrl, alternative='two-sided')[1]
                else:
                    pval = ttest_ind(d_exp, d_ctrl, equal_var=False)[1]

                results.append({
                    "Base Name":         cond,
                    "Gene ID":           "CHLRE_01g025050v5",
                    "Gene Name":         gene,
                    "p-value":           pval,
                    "log2(Exp/Control)": log2fc
                })
        
        df = pd.DataFrame(results)        
        df.to_csv(out_path, sep="\t", index=False)
        self.fix_output_format(out_path)
        QMessageBox.information(self, "Ready", f"The results are recorded in {out_path}")

    def fix_output_format(self, out_path):

        df = pd.read_csv(out_path, sep="\t")

        df["Gene Name"] = df["Gene Name"].astype(str)
        df["Gene Name"] = df["Gene Name"].str.replace('"', '', regex=False)
        df["Gene Name"] = df["Gene Name"].str.replace('\n', '', regex=False)
        df["Gene Name"] = df["Gene Name"].str.strip()



        if "Gene Name" in df.columns:
            df = df.rename(columns={"Gene Name": "GATA Name"})


        df = df.drop_duplicates(subset=["Base Name", "GATA Name"])

        df["GATA_num"] = df["GATA Name"].str.extract(r'(\d+)').astype(float)
        df = df.sort_values(by=["Base Name", "GATA_num"])
        df = df.drop(columns=["GATA_num"])

        cols = ["Base Name", "GATA Name", "p-value", "log2(Exp/Control)"]
        df = df[cols]
        df.to_csv(out_path, sep="\t", index=False)
        
if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = CtAnalysisApp()
    win.show()
    sys.exit(app.exec())