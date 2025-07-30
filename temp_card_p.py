import os
import json
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
from tkinter import Tk, Label, Entry, Button, StringVar, BooleanVar, Checkbutton, filedialog, Listbox
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm


script_dir = os.path.dirname(os.path.abspath(__file__))
settings_path = os.path.join(script_dir, "settings.json")
config_path = os.path.join(script_dir, "heatmap_config.json")


with open(settings_path, "r") as f:
    settings = json.load(f)
results_folder = settings["folders"]["results_folder"]


root_file = Tk()
root_file.withdraw()
selected_file = filedialog.askopenfilename(
    initialdir=results_folder,
    title="Выберите файл для тепловой карты",
    filetypes=(("Text files", "*.txt"), ("All files", "*.*"))
)
root_file.destroy()
if not selected_file:
    raise FileNotFoundError("Файл не выбран")
print(f"Используем файл: {selected_file}")


df = pd.read_csv(selected_file, sep="\t")
required_columns = ["Base Name", "GATA Name", "log2(Exp/Control)", "p-value"]
missing_columns = [col for col in required_columns if col not in df.columns]
if missing_columns:
    raise ValueError(f"Отсутствуют необходимые столбцы: {missing_columns}")

possible_conditions = sorted(df["Base Name"].unique())
possible_genes = sorted(
    df["GATA Name"].unique(),
    key=lambda x: (int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else float('inf'), x)
)


def load_or_ask_config(possible_conditions, possible_genes):

    config = {
        "title": "Heatmap of Log2 Fold Change",
        "xlabel": "GATA factors",
        "ylabel": "Conditions",
        "color_low": "blue",
        "color_mid": "green",
        "color_high": "red",
        "show_p_values": True,
        "custom_condition_order": possible_conditions,
        "custom_gene_order": possible_genes,
        "sort_conditions_by_max_exp": False,
        "sort_genes_by_max_exp": False,
        "remove_color_for_high_p": False
    }
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            file_config = json.load(f)
            config.update(file_config)

    root = Tk()
    root.title("Настройки тепловой карты")


    title_var = StringVar(value=config.get("title", ""))
    xlabel_var = StringVar(value=config.get("xlabel", ""))
    ylabel_var = StringVar(value=config.get("ylabel", ""))
    color_low_var = StringVar(value=config.get("color_low", ""))
    color_mid_var = StringVar(value=config.get("color_mid", ""))
    color_high_var = StringVar(value=config.get("color_high", ""))
    show_p_values_var = BooleanVar(value=config.get("show_p_values", True))
    sort_conditions_var = BooleanVar(value=config.get("sort_conditions_by_max_exp", False))
    sort_genes_var = BooleanVar(value=config.get("sort_genes_by_max_exp", False))
    remove_color_var = BooleanVar(value=config.get("remove_color_for_high_p", False))


    Label(root, text="Заголовок:").grid(row=0, column=0, sticky="w")
    Entry(root, textvariable=title_var).grid(row=0, column=1, columnspan=2, sticky="we")

    Label(root, text="Подпись X (вверху):").grid(row=1, column=0, sticky="w")
    Entry(root, textvariable=xlabel_var).grid(row=1, column=1, columnspan=2, sticky="we")

    Label(root, text="Подпись Y (слева):").grid(row=2, column=0, sticky="w")
    Entry(root, textvariable=ylabel_var).grid(row=2, column=1, columnspan=2, sticky="we")

    Label(root, text="Цвет min (нижний):").grid(row=3, column=0, sticky="w")
    Entry(root, textvariable=color_low_var).grid(row=3, column=1, columnspan=2, sticky="we")

    Label(root, text="Цвет 0 (середина):").grid(row=4, column=0, sticky="w")
    Entry(root, textvariable=color_mid_var).grid(row=4, column=1, columnspan=2, sticky="we")

    Label(root, text="Цвет max (верхний):").grid(row=5, column=0, sticky="w")
    Entry(root, textvariable=color_high_var).grid(row=5, column=1, columnspan=2, sticky="we")

    Checkbutton(root, text="Показывать p-value", variable=show_p_values_var).grid(row=6, column=0, columnspan=3, sticky="w")


    summary_label = Label(root, text="", justify="left", wraplength=400)
    summary_label.grid(row=7, column=0, columnspan=3, sticky="w")


    Checkbutton(root, text="Сортировать условия по максимальной экспрессии (log2)", variable=sort_conditions_var).grid(row=8, column=0, columnspan=3, sticky="w")
    Checkbutton(root, text="Сортировать гены по максимальной экспрессии (Control всегда первым)", variable=sort_genes_var).grid(row=9, column=0, columnspan=3, sticky="w")
    Checkbutton(root, text="Убрать цвет у ячеек с p-value > 0.05 и уменьшить шрифт в 2 раза", variable=remove_color_var).grid(row=10, column=0, columnspan=3, sticky="w")

    def update_summary(*args):
        summary = (
            f"Текущие настройки:\n"
            f"Заголовок: {title_var.get()}\n"
            f"Подпись X: {xlabel_var.get()}\n"
            f"Подпись Y: {ylabel_var.get()}\n"
            f"Цвет min: {color_low_var.get()}\n"
            f"Цвет 0: {color_mid_var.get()}\n"
            f"Цвет max: {color_high_var.get()}\n"
            f"Показывать p-value: {show_p_values_var.get()}\n"
            f"Сортировать условия по max экспрессии: {sort_conditions_var.get()}\n"
            f"Сортировать гены по max экспрессии (Control первым): {sort_genes_var.get()}\n"
            f"Убрать цвет у p-value > 0.05 и уменьшить шрифт: {remove_color_var.get()}\n"
            f"(Порядок условий и генов будет сохранён при сохранении)"
        )
        summary_label.config(text=summary)

    title_var.trace("w", update_summary)
    xlabel_var.trace("w", update_summary)
    ylabel_var.trace("w", update_summary)
    color_low_var.trace("w", update_summary)
    color_mid_var.trace("w", update_summary)
    color_high_var.trace("w", update_summary)
    show_p_values_var.trace("w", update_summary)
    sort_conditions_var.trace("w", update_summary)
    sort_genes_var.trace("w", update_summary)
    remove_color_var.trace("w", update_summary)
    update_summary()


    Label(root, text="Порядок условий:").grid(row=11, column=0, sticky="w")
    condition_listbox = Listbox(root, selectmode="single", exportselection=False, height=6)
    condition_listbox.grid(row=11, column=1)
    for item in possible_conditions:
        condition_listbox.insert("end", item)

    def move_up_cond():
        selection = condition_listbox.curselection()
        if not selection or selection[0] == 0:
            return
        i = selection[0]
        text = condition_listbox.get(i)
        condition_listbox.delete(i)
        condition_listbox.insert(i - 1, text)
        condition_listbox.select_set(i - 1)
        update_summary()

    def move_down_cond():
        selection = condition_listbox.curselection()
        if not selection or selection[0] == condition_listbox.size() - 1:
            return
        i = selection[0]
        text = condition_listbox.get(i)
        condition_listbox.delete(i)
        condition_listbox.insert(i + 1, text)
        condition_listbox.select_set(i + 1)
        update_summary()

    Button(root, text="↑", command=move_up_cond).grid(row=11, column=2, sticky="n")
    Button(root, text="↓", command=move_down_cond).grid(row=11, column=2, sticky="s")


    Label(root, text="Порядок генов:").grid(row=12, column=0, sticky="w")
    gene_listbox = Listbox(root, selectmode="single", exportselection=False, height=6)
    gene_listbox.grid(row=12, column=1)
    for item in possible_genes:
        gene_listbox.insert("end", item)

    def move_up_gene():
        selection = gene_listbox.curselection()
        if not selection or selection[0] == 0:
            return
        i = selection[0]
        text = gene_listbox.get(i)
        gene_listbox.delete(i)
        gene_listbox.insert(i - 1, text)
        gene_listbox.select_set(i - 1)
        update_summary()

    def move_down_gene():
        selection = gene_listbox.curselection()
        if not selection or selection[0] == gene_listbox.size() - 1:
            return
        i = selection[0]
        text = gene_listbox.get(i)
        gene_listbox.delete(i)
        gene_listbox.insert(i + 1, text)
        gene_listbox.select_set(i + 1)
        update_summary()

    Button(root, text="↑", command=move_up_gene).grid(row=12, column=2, sticky="n")
    Button(root, text="↓", command=move_down_gene).grid(row=12, column=2, sticky="s")


    def save_and_close():
        config["title"] = title_var.get()
        config["xlabel"] = xlabel_var.get()
        config["ylabel"] = ylabel_var.get()
        config["color_low"] = color_low_var.get()
        config["color_mid"] = color_mid_var.get()
        config["color_high"] = color_high_var.get()
        config["show_p_values"] = show_p_values_var.get()
        config["custom_condition_order"] = [condition_listbox.get(i) for i in range(condition_listbox.size())]
        config["custom_gene_order"] = [gene_listbox.get(i) for i in range(gene_listbox.size())]
        config["sort_conditions_by_max_exp"] = sort_conditions_var.get()
        config["sort_genes_by_max_exp"] = sort_genes_var.get()
        config["remove_color_for_high_p"] = remove_color_var.get()
        with open(config_path, "w") as f:
            json.dump(config, f, indent=2)
        root.destroy()

    Button(root, text="Сохранить и продолжить", command=save_and_close).grid(row=13, column=0, columnspan=3, pady=10)

    root.mainloop()
    return config

cfg = load_or_ask_config(possible_conditions, possible_genes)


df["Base Name"] = df["Base Name"].str.replace(r"(_rep\d+|_\d+)$", "", regex=True)


if cfg.get("sort_conditions_by_max_exp", False):
    condition_max = df.groupby("Base Name")["log2(Exp/Control)"].max()
    base_order = list(condition_max.sort_values(ascending=False).index)
else:
    base_order = cfg.get("custom_condition_order") or sorted(df["Base Name"].unique())

if cfg.get("sort_genes_by_max_exp", False):
    gene_max = df.groupby("GATA Name")["log2(Exp/Control)"].max()
    control_gene = None
    for gene in gene_max.index:
        if gene.startswith("Control"):
            control_gene = gene
            break
    if control_gene:
        non_control_genes = [g for g in gene_max.index if g != control_gene]
        non_control_sorted = list(gene_max[non_control_genes].sort_values(ascending=False).index)
        gene_order = [control_gene] + non_control_sorted
    else:
        gene_order = list(gene_max.sort_values(ascending=False).index)
else:
    gene_order = cfg.get("custom_gene_order") or sorted(df["GATA Name"].unique())

df["Base Name"] = pd.Categorical(df["Base Name"], categories=base_order, ordered=True)
df["GATA Name"] = pd.Categorical(df["GATA Name"], categories=gene_order, ordered=True)

heatmap_data = df.pivot(index="Base Name", columns="GATA Name", values="log2(Exp/Control)")
annotations = df.pivot(index="Base Name", columns="GATA Name", values="p-value")

vmax = max(abs(heatmap_data.min().min()), abs(heatmap_data.max().max()))
vmin = -vmax
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", [cfg["color_low"], cfg["color_mid"], cfg["color_high"]])
norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

if cfg["show_p_values"]:
    def format_pvalue(log2_val, p_val):

        return f"{log2_val:.1f}\n({p_val:.2f})" if p_val >= 0.05 else f"*{log2_val:.1f}*\n({p_val:.2f})"
    annot_text = np.vectorize(format_pvalue)(heatmap_data.values, annotations.values)
else:
    annot_text = np.vectorize(lambda x: f"{x:.1f}")(heatmap_data.values)

plt.figure(figsize=(max(12, len(gene_order) * 0.8), max(8, len(base_order) * 0.5)))

ax = sns.heatmap(
    heatmap_data,
    cmap=custom_cmap,
    norm=norm,
    annot=annot_text,
    fmt="",
    linewidths=0.5,
    cbar=True,
    linecolor='black',
    square=True
)


if cfg["show_p_values"]:
    for i in range(heatmap_data.shape[0]):
        for j in range(heatmap_data.shape[1]):
            if annotations.iloc[i, j] < 0.05:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='black', lw=1))


if cfg.get("remove_color_for_high_p", False):
    for i in range(heatmap_data.shape[0]):
        for j in range(heatmap_data.shape[1]):
            if annotations.iloc[i, j] > 0.05:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, fill=True, color="white", lw=0))
    for text in ax.texts:
        x, y = text.get_position()
        col = int(x)
        row = int(y)
        if row < heatmap_data.shape[0] and col < heatmap_data.shape[1]:
            if annotations.iloc[row, col] > 0.05:
                text.set_fontsize(text.get_fontsize() / 2)


cbar = ax.collections[0].colorbar
cbar.ax.yaxis.set_ticks_position("left")
scale_labels = {
    6: "+64x", 5: "+32x", 4: "+16x", 3: "+8x", 2: "+4x", 1: "+2x", 
    0: "no change", 
    -1: "÷2", -2: "÷4", -3: "÷8", -4: "÷16", -5: "÷32", -6: "÷64"
}
for value, label in scale_labels.items():
    if vmin <= value <= vmax:
        pos = (value - vmin) / (vmax - vmin)
        cbar.ax.text(1.2, pos, label, ha='left', va='center', fontsize=9, transform=cbar.ax.transAxes)

plt.title(cfg["title"], fontsize=14, fontweight="bold")
plt.xlabel(cfg["xlabel"], fontsize=12)
plt.ylabel(cfg["ylabel"], fontsize=12)
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()

def sanitize_filename(name):
    return re.sub(r'[\\/*?:"<>|]', "_", name)

filename = sanitize_filename(cfg["title"]).strip()
if not filename:
    filename = "heatmap"

output_path = os.path.join(results_folder, f"{filename}.png")
plt.savefig(output_path, dpi=300)
plt.show()

print(f"Тепловая карта сохранена: {output_path}")
