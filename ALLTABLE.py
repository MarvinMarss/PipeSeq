# ALLTABLE.py
import os
import sys
import re
import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QTableWidget, QTableWidgetItem, QFileDialog,
    QMessageBox, QCheckBox, QLabel, QHeaderView, QInputDialog,
    QComboBox, QDoubleSpinBox,
    QDialog, QListWidget, QListWidgetItem
)
from PyQt6.QtCore import Qt

# ---------- Визуальные настройки ----------
sns.set_style("whitegrid")
sns.set_context("talk", font_scale=0.9)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SETTINGS_PATH = os.path.join(SCRIPT_DIR, "settings.json")

# SciPy опционален
try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False

# Требуемые семантики столбцов
REQ_SEMANTICS = {
    "base_name":  [r"^base\s*_?\s*name$", r"^basename$"],
    "gene":       [r"^gata\s*_?\s*name$", r"^gene$", r"^symbol$"],
    "pvalue":     [r"^p-?value$", r"^p[_\s]?val(ue)?$"],
    "log2fc":     [r"^log2\s*\(\s*exp\s*/\s*control\s*\)$", r"^log2fc$", r"^log2\s*fold\s*change$"],
}

# ---------- Вспомогательные функции ----------

def load_results_folder():
    if os.path.exists(SETTINGS_PATH):
        try:
            with open(SETTINGS_PATH, "r", encoding="utf-8") as f:
                st = json.load(f)
            rf = st.get("folders", {}).get("results_folder", "")
            if rf and os.path.isdir(rf):
                return rf
        except Exception:
            pass
    return os.getcwd()

def find_column(df: pd.DataFrame, patterns) -> str | None:
    cols = list(df.columns)
    norm = {c: re.sub(r"\s+", " ", str(c)).strip() for c in cols}
    for pat in patterns:
        rx = re.compile(pat, re.IGNORECASE)
        for c, cn in norm.items():
            if rx.match(cn):
                return c
    return None

def resolve_columns(df: pd.DataFrame):
    col_base = find_column(df, REQ_SEMANTICS["base_name"])
    col_gene = find_column(df, REQ_SEMANTICS["gene"])
    col_pval = find_column(df, REQ_SEMANTICS["pvalue"])
    col_l2fc = find_column(df, REQ_SEMANTICS["log2fc"])

    missing = []
    if col_base is None: missing.append("Base Name")
    if col_gene is None: missing.append("GATA Name")
    if col_pval is None: missing.append("p-value")
    if col_l2fc is None: missing.append("log2(Exp/Control)")

    if missing:
        raise ValueError(f"Отсутствуют требуемые столбцы: {', '.join(missing)}")

    return col_base, col_gene, col_pval, col_l2fc

_time_rx = re.compile(r"(?P<val>\d+(?:[.,]\d+)?)\s*(?:h|hour|hours)\b", re.IGNORECASE)

def _extract_time_hours(s: str) -> str | None:
    m = _time_rx.search(s)
    if not m:
        m2 = re.search(r"(\d+(?:[.,]\d+)?)\s*(?:hour|hours)\b", s, re.IGNORECASE)
        if not m2:
            return None
        val = m2.group(1).replace(",", ".")
    else:
        val = m.group("val").replace(",", ".")
    try:
        x = float(val)
        return f"{int(x)}h" if x.is_integer() else f"{x}h"
    except Exception:
        return None

def _detect_regime(s: str) -> str | None:
    text = s.lower()
    if re.search(r"high\s*[- ]?\s*light", text): return "high light"
    if re.search(r"\blight\b", text):             return "light"
    if re.search(r"\bdark(ness)?\b", text):       return "dark"
    return None

def infer_condition_label(base_name: str, auto_normalize: bool) -> str:
    s = str(base_name).strip()
    s_norm = re.sub(r"\s+", " ", s)
    if not auto_normalize:
        return s_norm
    regime = _detect_regime(s_norm)
    t = _extract_time_hours(s_norm)
    if regime and t:   return f"{t} {regime}"
    if regime and not t: return regime
    if t and not regime: return t
    return s_norm

def collapse_gata_4(index_series: pd.Index) -> pd.Index:
    collapsed = []
    for g in index_series:
        g = str(g)
        g2 = re.sub(r"^(GATA-4)(?:_t\d+)?$", r"\1", g)
        collapsed.append(g2)
    return pd.Index(collapsed)

# ---------- Построение матриц ----------

def build_matrices(datasets, auto_normalize: bool):
    """
    Возвращает:
      combined_z : DataFrame (гены × (условие\\nметод)), Z-score по строкам (для теплокарты, ВСЕ эксперименты).
      M_log2     : DataFrame ((ген, ConditionNorm) × методы), log2FC (для совместимости; может не использоваться далее).
      M_p        : DataFrame ((ген, ConditionNorm) × методы), p-value.
      LONG       : DataFrame со строками на уровень наблюдения:
                   колонки: [GATA, ConditionNorm, BaseRaw, Method, log2FC, pval]
                   (используется в диалоге сопоставления по полным названиям и для реконструкции групп).
    """
    wide_blocks = []
    long_blocks = []

    for label, path in datasets:
        df = pd.read_csv(path, sep="\t")
        col_base, col_gene, col_pval, col_l2fc = resolve_columns(df)

        df = df.copy()
        base_raw = df[col_base].astype(str).map(lambda x: re.sub(r"\s+", " ", x.strip()))
        df["__Condition__"] = df[col_base].astype(str).map(lambda x: infer_condition_label(x, auto_normalize))
        df[col_gene] = collapse_gata_4(df[col_gene])

        # wide (теплокарта)
        pivot_val = df.pivot_table(index=col_gene, columns="__Condition__", values=col_l2fc, aggfunc="mean")
        pivot_p   = df.pivot_table(index=col_gene, columns="__Condition__", values=col_pval, aggfunc="min")
        pivot_p   = pivot_p.reindex(columns=pivot_val.columns)
        pivot_val.columns = [f"{c}\n{label}" for c in pivot_val.columns]
        wide_blocks.append(pivot_val)

        # long (для корреляций/сопоставления)
        tmp = pd.DataFrame({
            "GATA": df[col_gene].values,
            "ConditionNorm": df["__Condition__"].values,
            "BaseRaw": base_raw.values,
            "Method": [label]*len(df),
            "log2FC": df[col_l2fc].values,
            "pval":   df[col_pval].values
        })
        long_blocks.append(tmp)

    combined = pd.concat(wide_blocks, axis=1)
    mean = combined.mean(axis=1)
    std  = combined.std(axis=1).replace(0, np.nan)
    combined_z = combined.sub(mean, axis=0).div(std, axis=0)

    LONG = pd.concat(long_blocks, ignore_index=True)

    # совместимые «старые» своды по нормализованным условиям (могут не пригодиться далее)
    M_log2 = LONG.pivot_table(index=["GATA", "ConditionNorm"], columns="Method", values="log2FC", aggfunc="mean").sort_index()
    M_p    = LONG.pivot_table(index=["GATA", "ConditionNorm"], columns="Method", values="pval",   aggfunc="min").reindex(index=M_log2.index, columns=M_log2.columns)

    return combined_z, M_log2, M_p, LONG

# ---------- Корреляции (автовыбор метода) ----------

def _iqr_outlier_fraction(x: np.ndarray) -> float:
    x = x[~np.isnan(x)]
    if x.size < 4: return 0.0
    q1, q3 = np.percentile(x, [25, 75])
    iqr = q3 - q1
    if iqr == 0: return 0.0
    lo = q1 - 1.5 * iqr
    hi = q3 + 1.5 * iqr
    return float(np.mean((x < lo) | (x > hi)))

def _normality_ok(x: np.ndarray) -> bool:
    x = x[~np.isnan(x)]
    if x.size < 8: return False
    if SCIPY_AVAILABLE:
        try:
            p = stats.shapiro(x).pvalue if 8 <= x.size <= 5000 else stats.normaltest(x).pvalue
            return bool(p > 0.05)
        except Exception:
            pass
    uniq = np.unique(np.round(x, 3)).size
    return (uniq / max(1, x.size)) > 0.5

def _pearson(x, y):
    r = np.corrcoef(x, y)[0, 1]
    p = np.nan
    if SCIPY_AVAILABLE:
        try:
            r, p = stats.pearsonr(x, y)
        except Exception:
            p = np.nan
    return float(r), (float(p) if np.isfinite(p) else np.nan)

def _spearman(x, y):
    if SCIPY_AVAILABLE:
        try:
            res = stats.spearmanr(x, y, nan_policy='omit')
            return float(res.correlation), float(res.pvalue)
        except Exception:
            pass
    sx = pd.Series(x).rank().to_numpy()
    sy = pd.Series(y).rank().to_numpy()
    r = np.corrcoef(sx, sy)[0, 1]
    return float(r), np.nan

def _kendall(x, y):
    if SCIPY_AVAILABLE:
        try:
            res = stats.kendalltau(x, y, nan_policy='omit')
            r = float(res.correlation) if res.correlation is not None else np.nan
            p = float(res.pvalue) if res.pvalue is not None else np.nan
            return r, p
        except Exception:
            return np.nan, np.nan
    return np.nan, np.nan

def _auto_method(x, y):
    n = len(x)
    if n < 3:  return "kendall" if SCIPY_AVAILABLE and n >= 2 else "spearman"
    if n < 10: return "kendall" if SCIPY_AVAILABLE else "spearman"
    out_frac = max(_iqr_outlier_fraction(x), _iqr_outlier_fraction(y))
    norm_ok = _normality_ok(x) and _normality_ok(y)
    return "pearson" if (norm_ok and out_frac < 0.10) else "spearman"

def compute_pairwise_correlation_matrices(M_log2: pd.DataFrame,
                                          M_p: pd.DataFrame,
                                          alpha: float,
                                          pair_rule: str,
                                          corr_mode: str):
    methods = list(M_log2.columns)
    m = len(methods)
    R  = pd.DataFrame(np.nan, index=methods, columns=methods)
    P  = pd.DataFrame(np.nan, index=methods, columns=methods)
    N  = pd.DataFrame(0,     index=methods, columns=methods, dtype=int)
    Mth = pd.DataFrame("NA", index=methods, columns=methods, dtype=object)

    for i in range(m):
        R.iloc[i, i] = 1.0
        P.iloc[i, i] = 0.0
        N.iloc[i, i] = M_log2[methods[i]].notna().sum()
        Mth.iloc[i, i] = "—"

    for i in range(m):
        for j in range(i+1, m):
            xi = M_log2.iloc[:, i].to_numpy()
            xj = M_log2.iloc[:, j].to_numpy()
            pi = M_p.iloc[:, i].to_numpy()
            pj = M_p.iloc[:, j].to_numpy()

            mask = np.isfinite(xi) & np.isfinite(xj)
            if pair_rule == 'both':
                mask &= np.isfinite(pi) & np.isfinite(pj) & (pi <= alpha) & (pj <= alpha)
            elif pair_rule == 'any':
                any_ok = ((np.isfinite(pi) & (pi <= alpha)) | (np.isfinite(pj) & (pj <= alpha)))
                mask &= any_ok
            elif pair_rule == 'none':
                pass
            else:
                raise ValueError("Unknown pair_rule")

            x = xi[mask]; y = xj[mask]
            n = x.size
            N.iloc[i, j] = n; N.iloc[j, i] = n

            if n < 2 or np.all(x == x[0]) or np.all(y == y[0]):
                r = np.nan; p = np.nan; code = "NA"
            else:
                meth = corr_mode if corr_mode != 'auto' else _auto_method(x, y)
                def _compute(method_name):
                    if method_name == 'pearson':  return _pearson(x, y), 'P'
                    if method_name == 'spearman': return _spearman(x, y), 'S'
                    if method_name == 'kendall':  return _kendall(x, y), 'K'
                    return ((np.nan, np.nan), 'NA')
                (r, p), code = _compute(meth)
                if not np.isfinite(r):
                    for fb in (['spearman','kendall','pearson'] if meth=='pearson'
                               else ['pearson','kendall','spearman'] if meth=='spearman'
                               else ['spearman','pearson','kendall']):
                        (r, p), code = _compute(fb)
                        if np.isfinite(r): break

            R.iloc[i, j] = r; R.iloc[j, i] = r
            P.iloc[i, j] = p; P.iloc[j, i] = p
            Mth.iloc[i, j] = code; Mth.iloc[j, i] = code

    return R, P, N, Mth

# ---------- Визуализация ----------

def plot_heatmap(combined_z: pd.DataFrame, out_path_png: str):
    plt.figure(figsize=(max(12, combined_z.shape[1]*0.4), max(6, combined_z.shape[0]*0.4)))
    vmax = np.nanmax(np.abs(combined_z.values))
    vmax = 1.0 if not np.isfinite(vmax) else max(1.0, float(vmax))
    ax = sns.heatmap(
        combined_z, cmap="vlag", center=0, vmin=-vmax, vmax=vmax,
        linewidths=0.5, linecolor="gray", cbar_kws={"label": "Z-score", "shrink": 0.7}
    )
    ax.set_title("Combined heatmap (Condition × Method), Z-score per gene", fontsize=14, pad=12)
    ax.set_ylabel("Genes", fontsize=12)
    ax.set_xlabel("")
    plt.xticks(rotation=60, ha="right"); plt.yticks(rotation=0)
    plt.tight_layout(); plt.savefig(out_path_png, dpi=300); plt.close()

def plot_corr_advanced(R: pd.DataFrame, P: pd.DataFrame, N: pd.DataFrame, Mth: pd.DataFrame, out_path_png: str):
    annot = R.copy()
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            if i == j:
                annot.iloc[i, j] = "1.00"
            else:
                r = R.iloc[i, j]
                annot.iloc[i, j] = f"{r:.2f}" if np.isfinite(r) else "NA"

    plt.figure(figsize=(max(6, 1.6*R.shape[1]), max(5, 1.3*R.shape[0])))
    ax = sns.heatmap(R, annot=annot, fmt="", cmap="vlag", center=0,
                     square=True, cbar_kws={"label": "r"}, linewidths=0.5, linecolor="gray")
    ax.set_title("Correlation between methods (r)", fontsize=14, pad=12)  # необязательно менять заголовок
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(out_path_png, dpi=300)
    plt.close()

def plot_corr_methods(Mth: pd.DataFrame, N: pd.DataFrame, out_path_png: str):
    annot = Mth.copy().astype(str)
    for i in range(Mth.shape[0]):
        for j in range(Mth.shape[1]):
            annot.iloc[i, j] = "—" if i == j else f"{Mth.iloc[i, j]}\n(n={int(N.iloc[i, j])})"
    map_code = {"P":1.0, "S":0.5, "K":0.0, "—":np.nan, "NA":np.nan}
    Mnum = Mth.replace(map_code)
    plt.figure(figsize=(max(6, 1.6*Mth.shape[1]), max(5, 1.3*Mth.shape[0])))
    ax = sns.heatmap(Mnum, annot=annot, fmt="", cmap="coolwarm", center=0.5,
                     square=True, cbar=False, linewidths=0.5, linecolor="gray")
    ax.set_title("Chosen correlation method per pair (P/S/K) and sample size n", fontsize=14, pad=12)
    plt.xticks(rotation=45, ha="right"); plt.yticks(rotation=0)
    plt.tight_layout(); plt.savefig(out_path_png, dpi=300); plt.close()

def _find_rt_qpcr_label(methods: list[str]) -> str | None:
    rx = re.compile(r"rt[\s\-_]*q?pcr", re.IGNORECASE)
    for m in methods:
        if rx.search(m): return m
    return None

def _paired_xy(M_log2: pd.DataFrame, M_p: pd.DataFrame, a: float, rule: str, m1: str, m2: str):
    xi = M_log2[m1].to_numpy(); xj = M_log2[m2].to_numpy()
    pi = M_p[m1].to_numpy();   pj = M_p[m2].to_numpy()
    mask = np.isfinite(xi) & np.isfinite(xj)
    if rule == 'both':
        mask &= np.isfinite(pi) & np.isfinite(pj) & (pi <= a) & (pj <= a)
    elif rule == 'any':
        mask &= ((np.isfinite(pi) & (pi <= a)) | (np.isfinite(pj) & (pj <= a)))
    elif rule == 'none':
        pass
    return xi[mask], xj[mask]

def plot_panel_c_scatter(M_log2: pd.DataFrame, M_p: pd.DataFrame,
                         alpha: float, pair_rule: str,
                         corr_mode: str, out_path_png: str):
    methods = list(M_log2.columns)
    target = _find_rt_qpcr_label(methods) or (methods[0] if len(methods) >= 2 else None)
    if target is None: return
    pairs = [(target, m) for m in methods if m != target]

    nplots = len(pairs)
    if nplots == 0: return
    ncol = 3; nrow = int(np.ceil(nplots / ncol))
    plt.figure(figsize=(6*ncol, 5*nrow))

    idx = 1
    for (m1, m2) in pairs:
        x, y = _paired_xy(M_log2, M_p, alpha, pair_rule, m1, m2)
        plt.subplot(nrow, ncol, idx)
        if x.size >= 2:
            plt.scatter(x, y, s=28, alpha=0.8)
            try:
                coef = np.polyfit(x, y, 1)
                xs = np.linspace(np.min(x), np.max(x), 100)
                ys = coef[0]*xs + coef[1]
                plt.plot(xs, ys, linewidth=2)
            except Exception:
                pass
            method_used = corr_mode if corr_mode != 'auto' else _auto_method(x, y)
            if method_used == 'pearson':   r, p = _pearson(x, y);  code = "P"
            elif method_used == 'spearman': r, p = _spearman(x, y); code = "S"
            else:                           r, p = _kendall(x, y);  code = "K"
            n = x.size
            plt.title(f"{m1} vs {m2}\nr={r:.2f} (p={p:.3g}) | {code}, n={n}")
        else:
            plt.title(f"{m1} vs {m2}\nНедостаточно пар")
        plt.xlabel(m1); plt.ylabel(m2); idx += 1

    plt.tight_layout(); plt.savefig(out_path_png, dpi=300); plt.close()

def plot_panel_d_bland_altman(M_log2: pd.DataFrame, M_p: pd.DataFrame,
                              alpha: float, pair_rule: str,
                              out_path_png: str):
    methods = list(M_log2.columns)
    target = _find_rt_qpcr_label(methods) or (methods[0] if len(methods) >= 2 else None)
    if target is None: return
    pairs = [(target, m) for m in methods if m != target]

    nplots = len(pairs)
    if nplots == 0: return
    ncol = 3; nrow = int(np.ceil(nplots / ncol))
    plt.figure(figsize=(6*ncol, 5*nrow))

    idx = 1
    for (m1, m2) in pairs:
        x, y = _paired_xy(M_log2, M_p, alpha, pair_rule, m1, m2)
        plt.subplot(nrow, ncol, idx)
        if x.size >= 2:
            mean = (x + y) / 2.0
            diff = (y - x)
            md = np.mean(diff)
            sd = np.std(diff, ddof=1) if diff.size > 1 else 0.0
            loa1 = md - 1.96*sd; loa2 = md + 1.96*sd
            plt.scatter(mean, diff, s=28, alpha=0.8)
            plt.axhline(md, linestyle='--'); plt.axhline(loa1, linestyle=':'); plt.axhline(loa2, linestyle=':')
            plt.title(f"Bland–Altman: {m1} vs {m2}\nΔ= {md:.2f} ± 1.96·SD")
            plt.xlabel("Mean"); plt.ylabel("Difference (m2 - m1)")
        else:
            plt.title(f"Bland–Altman: {m1} vs {m2}\nНедостаточно пар")
        idx += 1

    plt.tight_layout(); plt.savefig(out_path_png, dpi=300); plt.close()

# ---------- Диалоги выбора ----------

class ConditionSelectDialog(QDialog):
    """Множественный выбор условий/групп для корреляций (после сопоставления)."""
    def __init__(self, conditions: list[str], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Выбор групп экспериментов для корреляций")
        self.resize(520, 520)

        v = QVBoxLayout(self)
        self.listw = QListWidget()
        self.listw.setSelectionMode(QListWidget.SelectionMode.MultiSelection)

        uniq = sorted({str(c) for c in conditions})
        for c in uniq:
            it = QListWidgetItem(c)
            it.setSelected(True)
            self.listw.addItem(it)

        v.addWidget(QLabel("Отметьте группы, которые войдут в расчёт корреляций и панелей (c, d):"))
        v.addWidget(self.listw)

        h = QHBoxLayout()
        btn_all = QPushButton("Выбрать все")
        btn_none = QPushButton("Снять все")
        btn_ok = QPushButton("OK")
        btn_cancel = QPushButton("Отмена")
        h.addWidget(btn_all); h.addWidget(btn_none)
        h.addStretch(1)
        h.addWidget(btn_ok); h.addWidget(btn_cancel)
        v.addLayout(h)

        btn_all.clicked.connect(self._select_all)
        btn_none.clicked.connect(self._select_none)
        btn_ok.clicked.connect(self.accept)
        btn_cancel.clicked.connect(self.reject)

    def _select_all(self):
        for i in range(self.listw.count()):
            self.listw.item(i).setSelected(True)

    def _select_none(self):
        for i in range(self.listw.count()):
            self.listw.item(i).setSelected(False)

    def selected_conditions(self) -> list[str]:
        return [self.listw.item(i).text()
                for i in range(self.listw.count())
                if self.listw.item(i).isSelected()]

class ConditionPairDialog(QDialog):
    """
    Сопоставление экспериментов (по ПОЛНЫМ названиям BaseRaw + метод).
    Пользователь формирует группы эквивалентных экспериментов; группы используются для корреляций.
    """
    def __init__(self, records: list[tuple[str, str, str]], parent=None):
        """
        records: список кортежей (Method, BaseRaw, ConditionNorm)
        """
        super().__init__(parent)
        self.setWindowTitle("Сопоставление экспериментов (полные названия)")
        self.resize(920, 620)

        self.records = records  # список
        # внутренние id для записей
        self.ids = [f"{m}||{b}" for (m, b, c) in records]
        self.id2rec = {i: rec for i, rec in zip(self.ids, records)}

        v = QVBoxLayout(self)
        v.addWidget(QLabel("Выберите элементы (минимум 2) и нажмите «Сформировать группу».\n"
                           "Отображаются ПОЛНЫЕ исходные названия Base Name (без сокращения) с указанием метода."))

        # список всех элементов (полные имена)
        self.list_all = QListWidget()
        self.list_all.setSelectionMode(QListWidget.SelectionMode.MultiSelection)
        for i, (m, b, c) in zip(self.ids, records):
            # полное имя (метод + исходный BaseRaw); для справки добавим нормализованное в скобках
            text = f"[{m}]  {b}    (норм.: {c})"
            it = QListWidgetItem(text)
            it.setData(Qt.ItemDataRole.UserRole, i)
            self.list_all.addItem(it)
        v.addWidget(self.list_all)

        # панель управления группами
        h = QHBoxLayout()
        self.btn_make = QPushButton("Сформировать группу из выделенных")
        self.btn_remove = QPushButton("Удалить выбранную группу")
        self.btn_clear = QPushButton("Очистить все группы")
        h.addWidget(self.btn_make); h.addWidget(self.btn_remove); h.addWidget(self.btn_clear)
        v.addLayout(h)

        v.addWidget(QLabel("Группы (эквивалентные эксперименты):"))

        # список групп
        self.list_groups = QListWidget()
        v.addWidget(self.list_groups)

        # кнопки OK/Cancel
        h2 = QHBoxLayout()
        btn_ok = QPushButton("OK"); btn_cancel = QPushButton("Отмена")
        h2.addStretch(1); h2.addWidget(btn_ok); h2.addWidget(btn_cancel)
        v.addLayout(h2)

        # структуры групп: список словарей {'label': ..., 'members': set(ids)}
        self.groups = []

        # события
        self.btn_make.clicked.connect(self._make_group)
        self.btn_remove.clicked.connect(self._remove_group)
        self.btn_clear.clicked.connect(self._clear_groups)
        btn_ok.clicked.connect(self.accept)
        btn_cancel.clicked.connect(self.reject)

    def _make_group(self):
        sel = [self.list_all.item(i) for i in range(self.list_all.count()) if self.list_all.item(i).isSelected()]
        if len(sel) < 2:
            QMessageBox.warning(self, "Недостаточно", "Выберите минимум два элемента для формирования группы.")
            return
        members = set([it.data(Qt.ItemDataRole.UserRole) for it in sel])

        # проверка: элементы могут уже быть в группах — в таком случае объединим группы (union)
        to_merge = []
        for g in self.groups:
            if len(members & g['members']) > 0:
                to_merge.append(g)
        merged_members = set(members)
        merged_label_parts = []
        for g in to_merge:
            merged_members |= g['members']
            merged_label_parts.append(g['label'])
        # удаляем объединяемые группы
        for g in to_merge:
            self.groups.remove(g)

        # имя группы: спросим у пользователя; дефолт — конкатенация первых двух BaseRaw
        # соберём небольшую подпись из уникальных BaseRaw
        raws = []
        for mid in list(merged_members)[:3]:
            m, b, c = self.id2rec[mid]
            raws.append(b)
        default_label = " ≡ ".join(raws)
        text, ok = QInputDialog.getText(self, "Название группы", "Введите имя группы:",
                                        text=default_label[:120] if default_label else "Group")
        if not ok or not text.strip():
            text = default_label if default_label else "Group"
        label = text.strip()

        self.groups.append({'label': label, 'members': set(merged_members)})
        self._refresh_groups_view()

    def _remove_group(self):
        it = self.list_groups.currentItem()
        if not it:
            return
        idx = it.data(Qt.ItemDataRole.UserRole)
        if idx is None:
            return
        self.groups.pop(idx)
        self._refresh_groups_view()

    def _clear_groups(self):
        self.groups.clear()
        self._refresh_groups_view()

    def _refresh_groups_view(self):
        self.list_groups.clear()
        for k, g in enumerate(self.groups):
            # соберём краткое описание состава
            members_txt = []
            for mid in list(g['members'])[:4]:
                m, b, c = self.id2rec[mid]
                members_txt.append(f"[{m}] {b}")
            if len(g['members']) > 4:
                members_txt.append(f"... (+{len(g['members'])-4})")
            text = f"{g['label']}  ←  " + " | ".join(members_txt)
            it = QListWidgetItem(text)
            it.setData(Qt.ItemDataRole.UserRole, k)
            self.list_groups.addItem(it)

    def mapping(self) -> dict:
        """
        Возвращает словарь сопоставления: (Method||BaseRaw) -> group_label
        Только для элементов, попавших в группы.
        """
        mp = {}
        for g in self.groups:
            for mid in g['members']:
                mp[mid] = g['label']
        return mp

# ---------- GUI ----------

class AllTableApp(QWidget):
    def __init__(self):
        super().__init__()
        self.results_folder = load_results_folder()
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("ALLTABLE — объединение данных, сопоставление экспериментов и корреляции")
        self.resize(1140, 680)

        main = QVBoxLayout()

        lbl = QLabel(
            "Колонки = методы/датасеты.\n"
            "• Одинарный клик по строке «Файл» — выбрать файл.\n"
            "• Двойной клик по строке «Имя» — переименовать метку метода.\n"
            "• «Автонормализация названий условий» влияет только на теплокарту и изначальный свод.\n"
            "• Сопоставление экспериментов выполняется по ПОЛНЫМ Base Name; только сопоставленные группы\n"
            "  включаются в расчёт корреляций и панелей (c,d). Z-score строится по всем экспериментам."
        )
        lbl.setWordWrap(True)
        main.addWidget(lbl)

        self.table = QTableWidget(2, 2)
        self.table.setHorizontalHeaderLabels(["Dataset_1", "Dataset_2"])
        self.table.setVerticalHeaderLabels(["Имя", "Файл"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)

        self.table.setItem(0, 0, QTableWidgetItem("Method_1"))
        self.table.setItem(1, 0, QTableWidgetItem(""))
        self.table.setItem(0, 1, QTableWidgetItem("Method_2"))
        self.table.setItem(1, 1, QTableWidgetItem(""))

        self.table.cellClicked.connect(self.on_cell_clicked)
        self.table.cellDoubleClicked.connect(self.on_cell_double_clicked)
        main.addWidget(self.table)

        # Панель управления
        h = QHBoxLayout()
        btn_add = QPushButton("Добавить колонку"); btn_add.clicked.connect(self.add_column)
        btn_del = QPushButton("Удалить выбранную колонку"); btn_del.clicked.connect(self.remove_selected_column)

        self.chk_auto = QCheckBox("Автонормализация названий условий"); self.chk_auto.setChecked(True)

        # alpha
        alpha_box = QHBoxLayout()
        alpha_label = QLabel("α (p-value):")
        self.spn_alpha = QDoubleSpinBox(); self.spn_alpha.setRange(1e-6, 0.5); self.spn_alpha.setSingleStep(0.005)
        self.spn_alpha.setDecimals(6); self.spn_alpha.setValue(0.05)
        alpha_box.addWidget(alpha_label); alpha_box.addWidget(self.spn_alpha)

        # правило отбора
        rule_box = QHBoxLayout()
        rule_label = QLabel("Отбор пар для корреляции:")
        self.cmb_rule = QComboBox()
        self.cmb_rule.addItems(["Оба значимы (p ≤ α)", "Хотя бы один значим (p ≤ α)", "Без фильтра по p-value"])
        rule_box.addWidget(rule_label); rule_box.addWidget(self.cmb_rule)

        # метод корреляции
        meth_box = QHBoxLayout()
        meth_label = QLabel("Correlation method:")
        self.cmb_meth = QComboBox()
        self.cmb_meth.addItems(["Auto", "Pearson", "Spearman", "Kendall"])
        meth_box.addWidget(meth_label); meth_box.addWidget(self.cmb_meth)

        h.addWidget(btn_add); h.addWidget(btn_del); h.addStretch(1)
        h.addLayout(alpha_box); h.addLayout(rule_box); h.addLayout(meth_box); h.addWidget(self.chk_auto)
        main.addLayout(h)

        btn_run = QPushButton("Начать создание"); btn_run.clicked.connect(self.run_build)
        main.addWidget(btn_run, alignment=Qt.AlignmentFlag.AlignRight)

        self.setLayout(main)

    # --- Табличные действия ---
    def add_column(self):
        c = self.table.columnCount()
        self.table.insertColumn(c)
        self.table.setHorizontalHeaderItem(c, QTableWidgetItem(f"Dataset_{c+1}"))
        self.table.setItem(0, c, QTableWidgetItem(f"Method_{c+1}"))
        self.table.setItem(1, c, QTableWidgetItem(""))

    def remove_selected_column(self):
        col = self.table.currentColumn()
        if col < 0:
            QMessageBox.warning(self, "Удаление", "Сначала выделите колонку."); return
        if self.table.columnCount() <= 1:
            QMessageBox.warning(self, "Удаление", "Должна остаться хотя бы одна колонка."); return
        self.table.removeColumn(col)

    def on_cell_clicked(self, row, col):
        if row == 1:
            start_dir = self.results_folder
            file_name, _ = QFileDialog.getOpenFileName(
                self, "Выберите файл с данными", start_dir,
                "Text/TSV files (*.txt *.tsv);;All files (*.*)"
            )
            if file_name:
                self.table.setItem(row, col, QTableWidgetItem(file_name))

    def on_cell_double_clicked(self, row, col):
        if row == 0:
            cur = self.table.item(row, col).text() if self.table.item(row, col) else ""
            text, ok = QInputDialog.getText(self, "Имя метода", "Введите название метода/набора:", text=cur)
            if ok and text.strip():
                self.table.setItem(row, col, QTableWidgetItem(text.strip()))

    # --- Запуск анализа ---
    def run_build(self):
        try:
            datasets = self.collect_datasets()
            if len(datasets) == 0:
                QMessageBox.warning(self, "Нет данных", "Добавьте хотя бы один файл."); return

            auto_normalize = self.chk_auto.isChecked()
            alpha = float(self.spn_alpha.value())

            rule_map = {"Оба значимы (p ≤ α)": "both", "Хотя бы один значим (p ≤ α)": "any", "Без фильтра по p-value": "none"}
            pair_rule = rule_map[self.cmb_rule.currentText()]

            meth_map = {"Auto": "auto", "Pearson": "pearson", "Spearman": "spearman", "Kendall": "kendall"}
            corr_mode = meth_map[self.cmb_meth.currentText()]

            # 1) Построение теплокарты по всем данным
            combined_z, M_log2_norm, M_p_norm, LONG = build_matrices(datasets, auto_normalize)
            out_dir = self.results_folder if os.path.isdir(self.results_folder) else os.getcwd()
            heatmap_png = os.path.join(out_dir, "ALLTABLE_heatmap.png")
            plot_heatmap(combined_z, heatmap_png)

            # 2) Сопоставление экспериментов по ПОЛНЫМ названиям
            #    records: (Method, BaseRaw, ConditionNorm) — без дублей
            rec_df = LONG[["Method", "BaseRaw", "ConditionNorm"]].drop_duplicates().copy()
            records = list(rec_df.itertuples(index=False, name=None))  # (Method, BaseRaw, ConditionNorm)

            if len(records) < 2:
                QMessageBox.warning(self, "Недостаточно экспериментов",
                                    "Найдено менее двух уникальных экспериментов (по полным названиям)."); return

            dlg_map = ConditionPairDialog(records, self)
            if dlg_map.exec() != QDialog.DialogCode.Accepted:
                return
            mapping = dlg_map.mapping()  # ключ: "Method||BaseRaw" -> group_label
            if len(mapping) == 0:
                QMessageBox.warning(self, "Пустое сопоставление",
                                    "Не создано ни одной группы сопоставления. Корреляции не будут посчитаны."); return

            # 3) Построение сводов на основе групп сопоставления
            #    Оставляем ТОЛЬКО те наблюдения, которые попали в группы.
            LONG = LONG.copy()
            LONG["__MID__"] = LONG["Method"] + "||" + LONG["BaseRaw"]
            LONG["Group"] = LONG["__MID__"].map(mapping)
            LONG_sel = LONG[LONG["Group"].notna()].copy()

            if LONG_sel.empty:
                QMessageBox.warning(self, "Нет данных для корреляций",
                                    "После сопоставления не осталось данных."); return

            # Индекс для корреляций: (GATA, Group)
            M_log2_grp = LONG_sel.pivot_table(index=["GATA", "Group"], columns="Method", values="log2FC", aggfunc="mean").sort_index()
            M_p_grp    = LONG_sel.pivot_table(index=["GATA", "Group"], columns="Method", values="pval",   aggfunc="min").reindex(index=M_log2_grp.index, columns=M_log2_grp.columns)

            # 4) Дополнительный выбор групп (какие именно группы участвуют)
            groups_all = list(M_log2_grp.index.get_level_values("Group").unique())
            with open(os.path.join(out_dir, "groups_all.txt"), "w", encoding="utf-8") as f_all:
                f_all.write("\n".join(sorted(map(str, groups_all))))

            dlg_sel = ConditionSelectDialog(groups_all, self)
            if dlg_sel.exec() != QDialog.DialogCode.Accepted:
                return
            selected_groups = dlg_sel.selected_conditions()
            if len(selected_groups) == 0:
                QMessageBox.warning(self, "Пустой выбор", "Нужно выбрать хотя бы одну группу."); return
            with open(os.path.join(out_dir, "groups_selected.txt"), "w", encoding="utf-8") as f_sel:
                f_sel.write("\n".join(sorted(map(str, selected_groups))))

            mask_idx = M_log2_grp.index.get_level_values("Group").isin(set(selected_groups))
            M_log2_sel = M_log2_grp.loc[mask_idx]
            M_p_sel    = M_p_grp.loc[mask_idx]

            if M_log2_sel.shape[0] == 0:
                QMessageBox.warning(self, "Нет данных для корреляций",
                                    "После выбора групп не осталось пар для расчёта."); return

            # 5) Корреляции (расширенные) и панели (c,d) по выбранным ГРУППАМ
            R, P, N, Mth = compute_pairwise_correlation_matrices(
                M_log2_sel, M_p_sel, alpha=alpha, pair_rule=pair_rule, corr_mode=corr_mode
            )
            corr_png = os.path.join(out_dir, "ALLTABLE_corr.png")
            plot_corr_advanced(R, P, N, Mth, corr_png)

            corr_meth_png = os.path.join(out_dir, "ALLTABLE_corr_methods.png")
            plot_corr_methods(Mth, N, corr_meth_png)

            panel_c_png = os.path.join(out_dir, "ALLTABLE_panel_c_scatter.png")
            plot_panel_c_scatter(M_log2_sel, M_p_sel, alpha, pair_rule, corr_mode, panel_c_png)

            panel_d_png = os.path.join(out_dir, "ALLTABLE_panel_d_bland_altman.png")
            plot_panel_d_bland_altman(M_log2_sel, M_p_sel, alpha, pair_rule, panel_d_png)

            # 6) Детализированная выгрузка
            details = []
            methods = list(M_log2_sel.columns)
            for i in range(len(methods)):
                for j in range(i+1, len(methods)):
                    mi, mj = methods[i], methods[j]
                    details.append({
                        "method_i": mi, "method_j": mj,
                        "r": R.loc[mi, mj], "p": P.loc[mi, mj], "n": int(N.loc[mi, mj]),
                        "chosen": Mth.loc[mi, mj],
                        "pair_rule": pair_rule, "corr_mode_request": corr_mode
                    })
            det_df = pd.DataFrame(details)
            det_csv = os.path.join(out_dir, "ALLTABLE_corr_details.csv")
            det_df.to_csv(det_csv, index=False)

            QMessageBox.information(
                self, "Готово",
                "Файлы сохранены:\n"
                f"{heatmap_png}\n"
                f"{corr_png}\n"
                f"{corr_meth_png}\n"
                f"{panel_c_png}\n"
                f"{panel_d_png}\n"
                f"{det_csv}\n"
                f"{os.path.join(out_dir, 'groups_all.txt')}\n"
                f"{os.path.join(out_dir, 'groups_selected.txt')}"
            )
        except Exception as e:
            QMessageBox.critical(self, "Ошибка", f"{type(e).__name__}: {e}")

    def collect_datasets(self):
        pairs = []
        for j in range(self.table.columnCount()):
            label_item = self.table.item(0, j)
            path_item = self.table.item(1, j)
            label = (label_item.text().strip() if label_item else f"Method_{j+1}") or f"Method_{j+1}"
            path = (path_item.text().strip() if path_item else "")
            if not path: continue
            if not os.path.exists(path):
                raise FileNotFoundError(f"Файл не найден: {path}")
            pairs.append((label, path))
        return pairs

# ---------- Точка входа ----------
if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = AllTableApp()
    win.show()
    sys.exit(app.exec())
